% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

clear variables
close all
clc

%% Sensitivity analysis for the Methylomicrobium buryatense model

% Test of the function deFBA fot the M. buryatense model
load model_Mburyatense.mat
model.tf = 20;
model.h = 5;
model.N = round(model.tf/model.h);
model.nobj = 1;
model.ProtQuotaProdID = ['Cell_wall','DNA','RNA','Lipid','Protein','Carbohydrate'];
model.quotaInitial = [0.09127,0.04,0.097,0.09,0,0.0378,zeros(1,145),0]';
model.objectiveWeights = [zeros(size(model.quotaWeights));model.proteinWeights]';

model.Y0(5) = 30;
model.Y0(7) = 30;
model.Y0(8) = 30;

model.continuousMet = zeros(model.sizeYmet,1);
model.continuousMet(2) = 1; % Methane
model.continuousMet(4) = 1; % Oxygen

kLaCH4 = 35; % Initial guess
kLaO = kLaCH4*1.37; % According to literature --> see thesis Koen Michiels

model.massTransfer = zeros(model.sizeYmet,1);
model.massTransfer(4) = kLaO;
model.massTransfer(2) = kLaCH4;

model.saturation = zeros(model.sizeYmet,1); % Saturation concentrations of gaseous substartes and products at 30 °C and 1 bar
model.saturation(4) = 0.036*0.99567/32*1000; % g O2/ kg water, 0.99567 kg/L water density --> total g/L, / M --> total mmol/L
model.saturation(2) = 0.018*0.99567/16*1000; % g CH4/ kg water, 0.99567 kg/L water density --> total g/L, / M --> total mmol/L

model.volume = 1; % L, the liquid volume in the bioreactor

model.Y0(2) = model.saturation(2)*model.volume;
model.Y0(4) = model.saturation(4)*model.volume;

kcat = []; % For the M. buryatense model, the units are 1/h
for i = 1:size(model.Kcat_f,2)
    ind = find(model.Kcat_f(:,i));
    if i ~= size(model.Kcat_f,2)
        kcat = [kcat; model.Kcat_f(ind(1),i)];
    else
        kcat = [kcat; model.Kcat_f(ind,i)];
    end
end

L = model.Kcat_f;
model = rmfield(model,'Kcat_f');

start = tic;
y = deFBA(model,kcat,'b',1);
telapsed = toc(start);

%% Implementation of the Morris screening method

% Define lower and upper bounds for the catalytic constants
% For a first version, these are taken to be 50% less and more than the
% current value

lb = kcat.*0.5;
ub = kcat.*1.5;

% Define parameters of the Morris method
M = 15;
r = 10; % Number of trajectories
p = 4; % Number of levels between lower and upper bounds

[mu,sigma] = Morris(@ (kcat) deFBA(model,kcat,'b',1),lb,ub,M,p,r); % y12 for lactate production

%% Plot the results of the Morris method

plot(mu,sigma,'bo','markerSize',6)
xlabel('\mu^*','fontSize',16)
ylabel('\sigma','fontSize',16)

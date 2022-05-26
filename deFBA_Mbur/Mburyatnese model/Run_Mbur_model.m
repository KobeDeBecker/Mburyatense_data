% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

clear variables
close all
clc

load('model_Mburyatense4.mat')

%% Run deFBA
tf = 40;
N = 40;
model.tf = tf;
model.N = N;
model.h = model.tf/model.N;
model.epsilon = 10^(-4);

model.Y0(6) = 30; % Nitrate as theoretically present in  the reactor medium
model.Y0(8) = 7; % Phosphate as theoretically present in  the reactor medium
model.Y0(9) = 2.5; % Sulphate as theoretically present in  the reactor medium
model.Y0(3) = model.Y0(3)*model.volume;
model.Y0(5) = model.Y0(5)*model.volume;

r = run_deFBA(model,'cplex');
if isempty(r)
    dummy = toVector(emptySolutionStruct(model),model);
    r = zeros(size(dummy));
end
s = toStruct(r,model);

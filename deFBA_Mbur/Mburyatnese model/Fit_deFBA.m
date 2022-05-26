% This script constructs the fitting of the deFBA model to the biomass and
% DO curves in order to determine the most sensitive katalytic constants as
% determiend during the sensitivity analysis. In addition, the kLa-values
% will also be estimated

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

clear variables
close all
clc

%% Load data
% Biomass
OD2 = load('OD2.mat');
OD2 = OD2.OD;
OD3 = load('OD3.mat');
OD3 = OD3.OD;
OD4 = load('OD4.mat');
OD4 = OD4.ODAU;

indx2 = isnan(OD2(:,2));
indx3 = isnan(OD3(:,2));
indx4 = isnan(OD4(:,2));
OD2(indx2,:) = []; % Remove NAN measurements
OD3(indx3,:) = []; % Remove NAN measurements
OD4(indx4,:) = []; % Remove NAN measurements
OD2(:,2) = OD2(:,2)-min(OD2(:,2)); % Correction for blanc measurement
OD3(:,2) = OD3(:,2)-min(OD3(:,2)); % Correction for blanc measurement
OD4(:,2) = OD4(:,2)-min(OD4(:,2)); % Correction for blanc measurement

% Select the correct time points
lambda2 = 6.2098; %3.1159;
t2 = 29.8438; % 60.27; Extra correction due to stationary phase
indx_OD2 = find(OD2(:,1)>t2,1);
OD2a = OD2(1:indx_OD2-1,:);
indx_OD2 = find(OD2a(:,1)<lambda2);
OD2a(indx_OD2,:) = [];

lambda3 = 6.8730; %5.47831;
t3 = 24.4547;
indx_OD3 = find(OD3(:,1)>t3,1);
OD3a = OD3(1:indx_OD3-1,:);
indx_OD3 = find(OD3a(:,1)<lambda3);
OD3a(indx_OD3,:) = [];

lambda4 = 0.1509; %0.0123;
t4 = 15.1039;
indx_OD4 = find(OD4(:,1)>t4,1);
OD4a = OD4(1:indx_OD4-1,:);
indx_OD4 = find(OD4a(:,1)<lambda4);
OD4a(indx_OD4,:) = [];

% Convert OD to biomass
a = 0.6639; % Coefficient as determined from calibration
X1 = [OD2a(:,1)-lambda2 OD2a(:,2)*a*3];
X2 = [OD3a(:,1)-lambda3 OD3a(:,2)*a*3];
X3 = [OD4a(:,1)-lambda4 OD4a(:,2)*a*3];

% Remove zero values (not compatible with parameter estimation)
indx1 = X1(:,2)==0;
X1(indx1,:) = [];
indx2 = X2(:,2)==0;
X2(indx2,:) = [];
indx3= X3(:,2)==0;
X3(indx3,:) = [];

% Dissolved oxygen
DO2 = load('DO2.mat');
DO2 = DO2.DO;
DO3 = load('DO3.mat');
DO3 = DO3.DO;
DO4 = load('DO4.mat');
DO4 = DO4.DO;

indx2 = isnan(DO2(:,2));
indx3 = isnan(DO3(:,2));
indx4 = isnan(DO4(:,2));
DO2(indx2,:) = []; % Remove NAN measurements
DO3(indx3,:) = []; % Remove NAN measurements
DO4(indx4,:) = []; % Remove NAN measurements

% Time point selection
indx_DO2 = find(DO2(:,1)>t2,1);
DO2a = DO2(1:indx_DO2-1,:);
indx_DO2 = find(DO2a(:,1)<lambda2);
DO2a(indx_DO2,:) = [];

indx_DO3 = find(DO3(:,1)>t3,1);
DO3a = DO3(1:indx_DO3-1,:);
indx_DO3 = find(DO3a(:,1)<lambda3);
DO3a(indx_DO3,:) = [];

indx_DO4 = find(DO4(:,1)>t4,1);
DO4a = DO4(1:indx_DO4-1,:);
indx_DO4 = find(DO4a(:,1)<lambda4);
DO4a(indx_DO4,:) = [];

max2 = max(DO2a(:,2));
max3 = max(DO3a(:,2));
max4 = max(DO4a(:,2));
DO2a(:,2) = DO2a(:,2)/max2*100; % Set maximum value at 100%
DO3a(:,2) = DO3a(:,2)/max3*100; % Set maximum value at 100%
DO4a(:,2) = DO4a(:,2)/max4*100; % Set maximum value at 100%

% Final data storage in vectors for parameter estimation
DO1 = [DO2a(:,1)-lambda2 DO2a(:,2)];
DO2 = [DO3a(:,1)-lambda3 DO3a(:,2)];
DO3 = [DO4a(:,1)-lambda4 DO4a(:,2)];

% Manual correction on data of DO2!
% t = 42.66;
% indx = DO2(:,1) >= 42.66;
% DO2(indx,:) = [];


%% Parameter estimation - X
model = load('model_Mburyatense3.mat');
model = model.model;
model.volume = 3;

model.massTransfer(5) = 10^10;
model.massTransfer(3) = 0.7299*model.massTransfer(5); % Harsh assumption...

tf = 40;
N = 12;

% Select parameter starting points: [kcat_pMMO,kcat_acpSmal,kcat_G3P_deh,kcat_argsucc,kcat_pglu_deh,kcat_PlsX,kcat_PlsY,kcat_PlsC,w_glc]
x01 = [5000,500,500,1474,1000,500,500,500,0.18,0.012,0.023,0.023]; % Initial guesses
%x0 = ones(size(x01));
lb = zeros(size(x01));
ub = 10^6*3600.*ones(size(x01));
ub(9) = 100;
ub(10:12) = 10^(-1);

% param_num = [1,9,10,11,12];
% for i = 1:length(param_num) %1:length(x0)
%     ne = 10;
%     obj(i,:) = zeros(1,ne);
%     if param_num(i) == 9
%         param(i,:) = linspace(0,20,ne);
%     elseif param_num(i) == 8 || param_num(i) == 7 || param_num(i) == 6 || param_num(i) == 4 || param_num(i) == 5
%         param(i,:) = linspace(1,2000,ne);
%     elseif param_num(i) == 1
%         param(i,:) = linspace(1000,7000,ne);
%     elseif param_num(i) == 2 || param_num(i) == 3
%         param(i,:) = linspace(1,2000,ne);
%     elseif param_num(i) == 10 || param_num(i) == 11 || param_num(i) == 12
%         param(i,:) = linspace(0.001,0.1,ne);
%     else
%         param(i,:) = linspace(100,10^6,ne);
%     end
%     for j = 1:length(param(i,:))
%         p = x0;
%         p(param_num(i)) = param(i,j);
%         obj(i,j) = obj_deFBAfit_X(p,X1,X2,X3,model,x01,tf,N);
%     end
%     figure(param_num(i))
%     plot(param(i,:),obj(i,:),'b-','LineWidth',2)
%     xlabel('p','FontSize',16)
%     ylabel('Obj','FontSize',16)
% end

% kLa_CH4 cannot be estimated since dissolved methane concentrations are not available. Empirical literature relations will be used.
options = optimoptions('lsqnonlin','PlotFcn','optimplotresnorm');
[param_deFBA,SSE,resid,~,~,~,J] = lsqnonlin(@ (p) obj_deFBAfit_X(p,X1,X2,X3,model,x01,tf,N),x01,lb,ub,options);
% param_deFBA = ga(@ (p) obj_deFBAfit(p,X1,X2,X3,DO1,DO2,DO3,w,model),9,[],[],[],[],lb,[]);
% options = optimoptions('particleswarm','PlotFcn','pswplotbestf');
% param_deFBA = particleswarm(@ (p) obj_deFBAfit(p,X1,X2,X3,DO1,DO2,DO3,w,model),10,lb,ub,options);
% options = optimoptions('simulannealbnd','PlotFcn','saplotbestf');
% param_deFBA = simulannealbnd(@ (p) obj_deFBAfit(p,X1,X2,X3,DO1,DO2,DO3,w,model),x0,lb,ub);
%% Comparison graphs
% load('parameters_Mbur3.mat')

tf = 40;
N = 12;
[Xm1,~] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(10)],X1(:,1),DO1(:,1),model,tf,N);
[Xm2,~] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(11)],X2(:,1),DO2(:,1),model,tf,N);
[Xm3,~] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(12)],X3(:,1),DO3(:,1),model,tf,N);

% Plot X
figure
p1 = scatter(X1(:,1),X1(:,2),'b','Filled');
hold on
p2 = scatter(X2(:,1),X2(:,2),'g','Filled');
p3 = scatter(X3(:,1),X3(:,2),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(X1(:,1),Xm1,'-b','LineWidth',2)
plot(X2(:,1),Xm2,'-g','LineWidth',2)
plot(X3(:,1),Xm3,'-r','LineWidth',2)
xlabel('t [h]','FontSize',16)
ylabel('X [g_{DW}]','FontSize',16)

lgd = legend({'Measurements 1','Measurements 2','Measurements 3','Fit 1','Fit 2','Fit 3'},'Location','NorthWest');
box on

% Plot X on log-scale
figure
plot(X1(:,1),log10(X1(:,2)),'ob','LineWidth',2)
hold on
plot(X2(:,1),log10(X2(:,2)),'og','LineWidth',2)
plot(X3(:,1),log10(X3(:,2)),'or','LineWidth',2)
plot(X1(:,1),log10(Xm1),'-b','LineWidth',2)
plot(X2(:,1),log10(Xm2),'-g','LineWidth',2)
plot(X3(:,1),log10(Xm3),'-r','LineWidth',2)
xlabel('t [h]','FontSize',16)
ylabel('X [g_{DW}]','FontSize',16)

%% Statistical analysis
df = length(X1(:,1)) + length(X2(:,1)) + length(X3(:,1)) - length(x01);
MSE = SSE/df;
P = MSE./(J'*J);

%% Save parameter results
save 'parameters_Mbur3.mat'

%% Determine kLa based on DO measurements
load('parameters_Mbur3.mat')
kLa0 = 5;
lb_kLa = 0;
ub_kLa = 10^3;

tf = 40;
N = 12;

options = optimoptions('lsqnonlin','PlotFcn','optimplotresnorm');
[kLa,SSE_kLa,resid_kLa,~,~,~,J_kLa] = lsqnonlin(@ (p) obj_deFBAfit_DO(p,param_deFBA,DO1,DO2,DO3,model,tf,N),kLa0,lb_kLa,ub_kLa,options);

%% Comparison graphs DO
load('parameters_Mbur3.mat')

tf = 40;
N = 12;
[Xm1_kLa,DOm1,OC1] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(10),kLa],X1(:,1),DO1(:,1),model,tf,N);
[Xm2_kLa,DOm2,OC2] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(11),kLa],X2(:,1),DO2(:,1),model,tf,N);
[Xm3_kLa,DOm3,OC3] = deFBA_paramfit([param_deFBA(1:9),param_deFBA(12),kLa],X3(:,1),DO3(:,1),model,tf,N);

% Plot X
figure
p1 = scatter(X1(:,1),X1(:,2),'b','Filled');
hold on
p2 = scatter(X2(:,1),X2(:,2),'g','Filled');
p3 = scatter(X3(:,1),X3(:,2),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(X1(:,1),Xm1_kLa,'-b','LineWidth',2)
plot(X2(:,1),Xm2_kLa,'-g','LineWidth',2)
plot(X3(:,1),Xm3_kLa,'-r','LineWidth',2)
xlabel('t [h]','FontSize',16)
ylabel('X [g_{DW}]','FontSize',16)

lgd = legend({'Measurements 1','Measurements 2','Measurements 3','Fit 1','Fit 2','Fit 3'},'Location','NorthWest');
box on

% Plot DO
figure
p1 = scatter(DO1(:,1),DO1(:,2),'b','Filled');
hold on
p2 = scatter(DO2(:,1),DO2(:,2),'g','Filled');
p3 = scatter(DO3(:,1),DO3(:,2),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(DO1(:,1),DOm1,'-b','LineWidth',2)
plot(DO2(:,1),DOm2,'-g','LineWidth',2)
plot(DO3(:,1),DOm3,'-r','LineWidth',2)
xlabel('t [h]','FontSize',16)
ylabel('DO [-]','FontSize',16)

lgd = legend({'Measurements 1','Measurements 2','Measurements 3','Fit 1','Fit 2','Fit 3'},'Location','SouthWest');
box on

r1 = load('r_O2_CH4_2.mat');
r2 = load('r_O2_CH4_3.mat');
r3 = load('r_O2_CH4_4.mat');
r1 = r1.r_O2_CH4;
r2 = r2.r_O2_CH4;
r3 = r3.r_O2_CH4;

indx_r1 = find(r1(:,1)<lambda2);
r1(indx_r1,:) = [];
r1(:,1) = r1(:,1)-lambda2;
indx_r2 = find(r2(:,1)<lambda3);
r3(indx_r2,:) = [];
r2(:,1) = r2(:,1)-lambda3;
indx_r3 = find(r3(:,1)<lambda4);
r3(indx_r3,:) = [];
r3(:,1) = r3(:,1)-lambda4;

% Plot O2/CH4-ratio
figure
p1 = scatter(r1(:,1),r1(:,2),'b','Filled');
hold on
p2 = scatter(r2(:,1),r2(:,2),'g','Filled');
p3 = scatter(r3(:,1),r3(:,2),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(tf/N/2:tf/N:tf,OC1,'-k','LineWidth',2)
plot([0 tf],[1.06 1.06],'--k','LineWidth',2)
xlim([0 25])
ylim([0 3])
xlabel('t [h]','FontSize',16)
ylabel('O_2/CH_4 ratio [-]','FontSize',16)

lgd = legend({'Measurements 1','Measurements 2','Measurements 3','Simulation','Experimental mean value'},'Location','NorthEast');
box on

%% Statistical analysis
df_kLa = length(DO1(:,1)) + length(DO2(:,1)) + length(DO3(:,1)) - 1;
MSE_kLa = SSE_kLa/df_kLa;
P_kLa = MSE_kLa./(J_kLa'*J_kLa);

%% Save final results
save 'parameters_Mbur3.mat'

% This script constructs the fitting of the biomasss curve in order to
% determine lag-time, maximal growth rate and start of the stationary
% phase. In addition, the presence of a linear growth rate will be
% statistically validated.

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

clear variables
close all
clc

%% Load all three data sets

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
t2 = 60.27;
indx_OD2 = find(OD2(:,1)>t2,1);
OD2a = OD2(1:indx_OD2-1,:);

t3 = 50.99;
indx_OD3 = find(OD3(:,1)>t3,1);
OD3a = OD3(1:indx_OD3-1,:);

t4 = 42.47;
indx_OD4 = find(OD4(:,1)>t4,1);
OD4a = OD4(1:indx_OD4-1,:);

% Convert OD to biomass
a = 0.6639; % Coefficient as determined from calibration
X1 = [OD2a(:,1) OD2a(:,2)*a];
X2 = [OD3a(:,1) OD3a(:,2)*a];
X3 = [OD4a(:,1) OD4a(:,2)*a];

% Remove zero values (not compatible with parameter estimation)
indx1 = find(X1(:,2)==0);
X1(indx1,:) = [];
indx2 = find(X2(:,2)==0);
X2(indx2,:) = [];
indx3= find(X3(:,2)==0);
X3(indx3,:) = [];

%% Fit the curves with only one exponential growth phase
x0 = [-2;-2;-2;13;13;5;0.1;30;30;30]; % Select parameter starting points: [lnX01;lnX02;lnX03;lambda1;lambda2;lambda3;mu;ts1;ts2;ts3]
[param_1growth,SSE_R,~,~,~,~,J_R] = lsqnonlin(@ (p) obj_1growthphase(p,X1,X2,X3),x0,[-Inf;-Inf;-Inf;0;0;0;0;0;0;0],[+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf]);

t = 0:0.1:60;
p1 = scatter(X2(:,1),log(X2(:,2)),'g','Filled');
hold on
p2 = scatter(X1(:,1),log(X1(:,2)),'b','Filled');
p3 = scatter(X3(:,1),log(X3(:,2)),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(t,model_1growthphase([param_1growth(1);param_1growth(4);param_1growth(7);param_1growth(8)],t),'-b','LineWidth',2)
plot(t,model_1growthphase([param_1growth(2);param_1growth(5);param_1growth(7);param_1growth(9)],t),'-g','LineWidth',2)
plot(t,model_1growthphase([param_1growth(3);param_1growth(6);param_1growth(7);param_1growth(10)],t),'-r','LineWidth',2)

xlabel('t [h]','FontSize',16)
ylabel('ln(X) [/]','FontSize',16)

%% Fit curves with an additional linear growth phase
x0 = [-2;-2;-2;13;13;5;0.2;0.1;10;10;10;30;30;30]; % Select parameter starting points: [lnX01;lnX02;lnX03;lambda1;lambda2;lambda3;mu_exp;mu_lin;tl1;tl2;tl3;ts1;ts2;ts3]
[param_2growth,SSE_F,~,~,~,~,J_F] = lsqnonlin(@ (p) obj_2growthphase(p,X1,X2,X3),x0,[-Inf;-Inf;-Inf;0;0;0;0;0;0;0;0;0;0;0],[+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf,+Inf]);

figure
t = 0:0.1:60;
p1 = scatter(X2(:,1),log(X2(:,2)),'g','Filled');
hold on
p2 = scatter(X1(:,1),log(X1(:,2)),'b','Filled');
p3 = scatter(X3(:,1),log(X3(:,2)),'r','Filled');
p1.MarkerFaceAlpha = 0.25;
p2.MarkerFaceAlpha = 0.25;
p3.MarkerFaceAlpha = 0.25;
plot(t,model_2growthphase([param_2growth(1);param_2growth(4);param_2growth(7);param_2growth(8);param_2growth(9);param_2growth(12)],t),'-b','LineWidth',2)
plot(t,model_2growthphase([param_2growth(2);param_2growth(5);param_2growth(7);param_2growth(8);param_2growth(10);param_2growth(13)],t),'-g','LineWidth',2)
plot(t,model_2growthphase([param_2growth(3);param_2growth(6);param_2growth(7);param_2growth(8);param_2growth(11);param_2growth(14)],t),'-r','LineWidth',2)

xlabel('t [h]','FontSize',16)
ylabel('ln(X) [/]','FontSize',16)

%% Statistical analysis (F-test + variances)
% F-test
df_F = (size(X1,1)+size(X2,1)+size(X3,2)) - length(param_2growth); % Degrees of freedom of full model
df_R = (size(X1,1)+size(X2,1)+size(X3,2)) - length(param_1growth); % Degress of freedom of reduced model
alpha = 0.05; % Corresponding to 95% confidence
F_star = (SSE_R-SSE_F)/(df_R-df_F)*(df_F/SSE_F);
F = finv(1-alpha,df_R-df_F,df_F);
if F_star > F % Accept more complex model
    fprintf('The value for F^* equals %d. The value for F with a confidence of %d equals %d. \n The complex model is accepted. \n',F_star,(1-alpha)*100,F)
else % Reject more complex model
    fprintf('The value for F^* equals %d. The value for F with a confidence of %d equals %d. \n The complex model is rejected. \n',F_star,(1-alpha)*100,F)
end

% Variances
MSE_R = SSE_R/df_R;
MSE_F = SSE_F/df_F;
P_R = MSE_R./(J_R'*J_R);
P_F = MSE_F./(J_F'*J_F); 
% Infinite values for stationary switch time of X1 and X2 due to short
% measurements

%% Plot results in an exponential graph
figure
t = 0:0.1:60;
p1 = scatter(X2(:,1),X2(:,2),'b','Filled');
hold on
p2 = scatter(X1(:,1),X1(:,2),'g','Filled');
p3 = scatter(X3(:,1),X3(:,2),'r','Filled');
p1.MarkerFaceAlpha = 0.05;
p2.MarkerFaceAlpha = 0.05;
p3.MarkerFaceAlpha = 0.05;

plot(t,exp(model_2growthphase([param_2growth(1);param_2growth(4);param_2growth(7);param_2growth(8);param_2growth(9);param_2growth(12)],t)),'-b','LineWidth',2)
plot(t,exp(model_2growthphase([param_2growth(2);param_2growth(5);param_2growth(7);param_2growth(8);param_2growth(10);param_2growth(13)],t)),'-g','LineWidth',2)
plot(t,exp(model_2growthphase([param_2growth(3);param_2growth(6);param_2growth(7);param_2growth(8);param_2growth(11);param_2growth(14)],t)),'-r','LineWidth',2)

plot(t,exp(model_1growthphase([param_1growth(1);param_1growth(4);param_1growth(7);param_1growth(8)],t)),'--b','LineWidth',2)
plot(t,exp(model_1growthphase([param_1growth(2);param_1growth(5);param_1growth(7);param_1growth(9)],t)),'--g','LineWidth',2)
plot(t,exp(model_1growthphase([param_1growth(3);param_1growth(6);param_1growth(7);param_1growth(10)],t)),'--r','LineWidth',2)

xlabel('t [h]','FontSize',16)
ylabel('C_X [g_{DW/L}]','FontSize',16)
xlim([0 60])

lgd = legend({'Measurements 1','Measurements 2','Measurements 3','Two-phase growth 1','Two-phase growth 2','Two-phase growth 3','Exponential phase 1','Exponential phase 2','Exponential phase 3'},'Location','NorthWest');
box on
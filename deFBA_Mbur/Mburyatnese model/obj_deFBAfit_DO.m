function obj = obj_deFBAfit_DO(kLa,p,DO1,DO2,DO3,model,tf,N)
% This function descibes the objective function to be minimized by solving
% a non-linear least squares problem in order to determine the kLa values
% for M. buryatense fermentations, given the adjusted model - version 3.0.

% Written by Kobe De Becker - 14/12/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

% Determine parameter vector
p1 = [p(1:10),kLa];
p2 = [p(1:9),p(11),kLa];
p3 = [p(1:9),p(12),kLa];

% Simulate model results
[~,DO1_m] = deFBA_paramfit(p1,DO1(:,1),DO1(:,1),model,tf,N);
[~,DO2_m] = deFBA_paramfit(p2,DO2(:,1),DO2(:,1),model,tf,N);
[~,DO3_m] = deFBA_paramfit(p3,DO3(:,1),DO3(:,1),model,tf,N);

% Subtract from experimental data such that objective minimizes J = ||DO -
% DO_m||^2
F1 = DO1(:,2) - DO1_m;
F2 = DO2(:,2) - DO2_m;
F3 = DO3(:,2) - DO3_m;

% Construct objective error vector
obj = [F1;F2;F3]; 

end
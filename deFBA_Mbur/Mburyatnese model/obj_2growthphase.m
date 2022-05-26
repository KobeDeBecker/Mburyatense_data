% This function descibes the objective function to be minimized by solving
% a non-linear least squares problem in case one (exponetial) growth
% phase is present followed by a linear growth phase. 

% Written by Kobe De Becker - 25/01/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function obj = obj_2growthphase(p,X1,X2,X3)

% Allocate parameter values
p1 = [p(1);p(4);p(7);p(8);p(9);p(12)];
p2 = [p(2);p(5);p(7);p(8);p(10);p(13)];
p3 = [p(3);p(6);p(7);p(8);p(11);p(14)];

% Simulate model data
lnX1_m = model_2growthphase(p1,X1(:,1));
lnX2_m = model_2growthphase(p2,X2(:,1));
lnX3_m = model_2growthphase(p3,X3(:,1));

% Define weighting factors to minimze the effect of early measurements
% (more noisy)
W1 = ones(1,length(X1(:,1)));
indx = find(X1(:,1)>10,1);
W1(1:indx) = 0.1;

W2 = ones(1,length(X2(:,1)));
indx = find(X2(:,1)>10,1);
W2(1:indx) = 0.1;

W3 = ones(1,length(X3(:,1)));
indx = find(X3(:,1)>10,1);
W3(1:indx) = 0.1;

% Subtract from experimental data
F1 = W1'.*(log(X1(:,2)) - lnX1_m);
F2 = W2'.*(log(X2(:,2)) - lnX2_m);
F3 = W3'.*(log(X3(:,2)) - lnX3_m);

% Construct objective error vector
obj = [F1;F2;F3];
end
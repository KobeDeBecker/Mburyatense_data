% This function descibes the objective function to be minimized by solving
% a non-linear least squares problem in order to determine the most
% sensitive catalytic constants for M. buryatense
% fermentations, given the adjusted model - version 3.0.

% Written by Kobe De Becker - 29/01/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function obj = obj_deFBAfit_X(p,X1,X2,X3,model,x01,tf,N)

% Allocate initial biomass values
p1 = [p(1:9),p(10)];
p2 = [p(1:9),p(11)];
p3 = [p(1:9),p(12)];

% Simulate model data
[X1_m,~] = deFBA_paramfit(p1,X1(:,1),X1(:,1),model,tf,N);
[X2_m,~] = deFBA_paramfit(p2,X2(:,1),X2(:,1),model,tf,N);
[X3_m,~] = deFBA_paramfit(p3,X3(:,1),X3(:,1),model,tf,N);

% Take logarithm to work in the same order of magnitude
X1_y = log10(X1(:,2));
X2_y = log10(X2(:,2));
X3_y = log10(X3(:,2));
X1_m = log10(X1_m);
X2_m = log10(X2_m);
X3_m = log10(X3_m);

% Define weighting matrices
W1 = ones(1,length(X1_y));
W2 = ones(1,length(X2_y));
W3 = ones(1,length(X3_y));

ind1 = X1(:,1)<=10;
ind2 = X1(:,1)<=10;
ind3 = X1(:,1)<=10;
W1(ind1) = 0.1;
W2(ind2) = 0.1;
W3(ind3) = 0.1;

% Subtract from experimental data such that objective minimizes J = W*||X -
% X_m||^2
F1 = W1'.*(X1_y - X1_m);
F2 = W2'.*(X2_y - X2_m);
F3 = W3'.*(X3_y - X3_m);

% Construct objective error vector
obj = [F1;F2;F3];
%obj = norm([F1;F2;F3])^2;
end
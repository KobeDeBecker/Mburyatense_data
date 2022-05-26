function [mu, sigma] = Morris(model,lb,ub,M,p,r)
% This function executes the Morris sensitvity analysis.
% INPUTS:
%- model: the model with as inputs the paraeters for which a sensitvity
% analysis is required
% - lb: lower bound on parameter values
% - ub: upper bound on parameter values
% - M: number of Morris trajectories generated
% - p: number of grid points
% - r: number of chosen Morris trajectories
% OUTPUTS:
% - mu: a vector containing all averages of the elementary effects of a
% parameter
% - sigma: a vector containing all standard deviations of the elementary 
% effects of a parameter

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

%% Generate required constants

% Define number of parameters
k = size(lb,1);

% Create the stepsize of the elementary effects
delta = p/(2*(p-1));

% Create base values
base = 0;
v = 1;
while v/(p-1) < 1-delta
    base = [base, v/(p-1)];
    v = v+1;
end

% Create a matrix B with the following properties:
% - size k+1 by k witk k the number of parameters
% - all elements are either 0 or 1
% - for every column index j, j=1, …, k, there are two rows of B that
% differ only in the jth entry

B = zeros(k+1,k);
for j = 2:size(B,1)
    B(j,1:j-1) = ones(1,size(1:j-1,2));
end

%% Step 1: Generate M trajectories

X = zeros(M*(k+1),k);
for i = 1:M
    % Generate a random point xi in the hypercube on the p-level grid
    xi = zeros(1,k);
    for z = 1:k
        ind = randperm(size(base,2));
        xi(z) = base(ind(1));
    end
    
    % Generate a diogonal matrix D with on its diagonal only 1 or -1 with
    % equal probability
    d = zeros(1,k);
    for z = 1:k
        element = -1+2*rand;
        if element < 0
            element = -1;
        elseif element >= 0
            element = 1;
        end
        d(z) = element;
    end
    D = diag(d);
    
    % Generate a random permutation matrix P
    P = zeros(k,k);
    indx = randperm(k);
    for t = 1:length(indx)
        P(t,indx(t)) = 1;
    end
    
    % Calculate a sampling matrix Bi
    J = ones(k+1,k);
    Bi = ones(k+1,1)*xi + ((delta/2).*((2.*B-J)*D+J))*P;
    
    % Store the sampling matrix in X with size M*(k+1) by k
    X((i-1)*(k+1)+1:i*(k+1),:) = Bi;
end

%% Step 2: Choose r trajectories with the highest dispersion in x-space

% Calculate all pairwise distances between the trajectories
distances = zeros(M,M);
for i=1:M
    for j=i:M
        if i == j
            d = 0;
        else
            d = determinedistance(X((i-1)*(k+1)+1:i*(k+1),:),X((j-1)*(k+1)+1:j*(k+1),:));
        end
        distances(i,j) = d;
        distances(j,i) = d;
    end
end

% Choose the ones with the highest "spread Ds" (r trajectories are selected)
comb = combnk(1:M,r);
Ds = zeros(size(comb,1),1);
for i = 1:size(comb,1)
    comb1 = combnk(comb(i,:),2);
    for j = 1:size(comb1,1)
        Ds1 = distances(comb1(j,1),comb1(j,2));
        Ds(i) = Ds(i) + Ds1^2;
    end
    Ds(i) = sqrt(Ds(i));
end
[~,ind] = sort(Ds,'descend');
indmax = ind(1);
runs = comb(indmax,:);

disp('Generated and chosen Morris trajectories')

%% Step 3: Calculate elementary effects

EEtot = zeros(k,length(runs));
for t = 1:length(runs)
    Xi = X((runs(t)-1)*(k+1)+1:runs(t)*(k+1),:);
    EE = determineelementaryeffects(model,Xi,lb,ub);
    EEtot(:,t) = EE;
end

%% Step 4: Calculate mu and sigma

% Note that the absolute value of the elementary effects is taken while
% calculating mu!
mu = zeros(size(EEtot,1),1);
sigma = zeros(size(EEtot,1),1);
for i = 1:size(EEtot,1)
    mu(i) = mean(abs(EEtot(i,:)));
    sigma(i) = std(EEtot(i,:));
end

end
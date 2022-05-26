function EE = determineelementaryeffects(model,r,lb,ub)
% This function determines the elementary effects for one run.
% INPUTS:
% - model: the model to generate the outputs
% - r: the Morris trajectory (or run)
% - lb: lower bounds for parameter values
% - ub: upper bounds for parameter values
% OUTPUTS:
% - EE: a vector containing all the elemenatry effects of all the inputs of
% the model.

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

% Transform trajectories to the paramater coordinates
for i=1:size(r,1)
    for j=1:size(r,2)
        r(i,j) = (r(i,j))*(ub(j)-lb(j))+lb(j);
    end
end

% Solve all model outputs
output = zeros(size(r,1),1);
for s = 1:size(r,1)
    output1 = model(r(s,:));
    output(s,1) = output1;
end

% Find the normalized elementary effects
EE = zeros(size(lb,1),1);
for i = 2:size(r,1)
    ind = find(r(i,:)~=r(i-1,:));
    EE(ind) = (output(i,1)-output(i-1,1))/(r(i,ind)-r(i-1,ind))*((ub(ind)+lb(ind))/2)/(mean(output));
end

end
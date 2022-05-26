% This function generates the biomass amount in g_DW for parameter values
% given in p, at the time points given in tp. The model assumes one
% exponential growth phase

% Written by Kobe De Becker - 25/01/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function lnX = model_1growthphase(p,tp)

lnX0 = p(1); % Initial biomass amount
lambda = p(2); % lag time
mu = p(3); % Exponential growth rate
ts = p(4); % Start of stationary phase

lnX = zeros(length(tp),1);
for i = 1:length(tp)
    if tp(i) <= lambda % Lag pahse
        lnX(i) = lnX0;
    elseif tp(i) > lambda && tp(i) <= ts % Exponential growth phase
        lnX(i) = lnX0 + mu*(tp(i)-lambda); 
    else
        lnX(i) = lnX0 + mu*(ts-lambda);
    end
end
end
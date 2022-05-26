% This function generates the biomass amount in g_DW for parameter values
% given in p, at the time points given in tp. The model assumes one
% exponential growth phase folowed by a linear growth rate;

% Written by Kobe De Becker - 25/01/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function lnX = model_2growthphase(p,tp)

lnX0 = p(1); % Initial biomass amount
lambda = p(2); % Lag time
mu_exp = p(3); % Exponential growth rate
mu_lin = p(4); % Linear growth rate
tl = p(5); % Start of linear growth pahse
ts = p(6); % Start of stationary phase

lnX = zeros(length(tp),1);
for i = 1:length(tp)
    if tp(i) <= lambda % Lag pahse
        lnX(i) = lnX0;
    elseif tp(i) > lambda && tp(i) <= tl % Exponential growth phase
        lnX(i) = lnX0 + mu_exp*(tp(i)-lambda); 
    elseif tp(i) > tl && tp(i) <= ts
        lnX(i) = log(exp(lnX0)*exp(mu_exp*(tl-lambda)) + mu_lin*(tp(i)-tl));
    else
        lnX(i) = log(exp(lnX0)*exp(mu_exp*(tl-lambda)) + mu_lin*(ts-tl));
    end
end
end
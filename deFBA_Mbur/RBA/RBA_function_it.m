function [P,S,mumax] = RBA_function_it(B,model,mu0,delta0)
% This function calculates the biomass composition distribution given a
% total amount of biomass and a compatible Resource Allocation Model.
% INPUTS: - B: Total biomass amount
%         - model: Model compatible with RBA code
% OUTPUT: - P: Biomass composition distribution
%         - S: Amount of storage compounds
%         - mu: Maximal expected growth rate

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

mu = mu0;
delta = delta0;
counter = 1;
while delta >= 10^(-4)
    [Aeq,beq,Aineq,bineq,lb,ub] = formulateConstraintsRBA(model,mu,B);

    cons = emptySolutionStructRBA(model);
    c = zeros(1,size(toVectorRBA(cons,model),2));
    
    % create cplex object
    lpProb = Cplex;
    % set parameters
    lpProb.Param.simplex.tolerances.optimality.Cur = 1e-9;
    lpProb.Param.simplex.tolerances.feasibility.Cur = 1e-9;
    lpProb.Param.feasopt.tolerance.Cur = 1e-9;
    lpProb.Param.sifting.display.Cur = 0;
    lpProb.Param.emphasis.numerical.Cur = 1;
    lpProb.Param.simplex.display.Cur = 0;
    lpProb.Param.tune.display.Cur = 3;
    lpProb.Param.barrier.convergetol.Cur = 1e-12;
    lpProb.Param.barrier.display.Cur = 0;
    lpProb.Param.threads.Cur = 12;
    lpProb.Param.workmem.Cur = 2048;     
    lpProb.Param.read.scale.Cur = 0;
    lpProb.DisplayFunc =[];

    % add bounds and objective
    lpProb.addCols(c',[],lb',ub');
    lpProb.Model.sense = 'maximize';
    lpProb.addRows([beq,-inf*ones(1,size(Aineq,1))]',[Aeq;Aineq],[beq,bineq]');

    lpProb.solve();
    if lpProb.Solution.status ~= 1 && counter == 1 % Initial mu does not result in a feasible solution
        mu = mu/10;
    elseif lpProb.Solution.status == 1 % Solution found
        result = lpProb.Solution.x;
        s = toStructRBA(result,model);
        P = s.p;
        S = s.s;
        mumax = mu;
        mu = mu + delta;
        counter = 2;
    elseif lpProb.Solution.status ~= 1 && counter ~= 1 % Infeasible or other problem
        delta = delta/10;
    end
end
if isempty(find(result)) % Check if it is not the zero flux solution
            disp('Trivial RBA solution obtained')
%         else
%             disp('Non-trivial RBA solution obtained')
end

end
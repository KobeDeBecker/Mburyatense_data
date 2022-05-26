function [P,S] = RBA_function(B,model,mu)
% This function calculates the biomass composition distribution given a
% total amount of biomass and a compatible Resource Allocation model.
% INPUTS: - B: Total biomass amount
%         - model: Model compatible with RBA code
%         - mu: Maximal expected growth rate
% OUTPUT: - P: Biomass composition distribution

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

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
if lpProb.Solution.status == 1 % Solution found
    result = lpProb.Solution.x;
    if isempty(find(result)) % Check if it is not the zero flux solution
        disp('Trivial RBA solution obtained')
    else
        disp('Non-trivial RBA solution obtained')
    end
    s = toStructRBA(result,model);
    P = s.p;
    S = s.s;
    
elseif lpProb.Solution.status ~= 1 % Infeasible or other problem
    disp('Infeasible/unbounded RBA problem, consider using other initial points.')
    if lpProb.Solution.status == 3 % Infeasible problem
        a = lpProb.refineConflict;
        indc = a.Conflict.colind;
        indr = a.Conflict.rowind;
    end
    P = [];
end

end
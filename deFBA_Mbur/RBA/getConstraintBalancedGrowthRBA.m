function [Aeq,beq] = getConstraintBalancedGrowthRBA(model,mu)

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

cons = emptySolutionStructRBA(model);
Aeq = sparse(model.sizePmet,size(toVectorRBA(cons,model),2));
beq = zeros(1,model.sizePmet);

 for i = 1:model.sizePmet
     idx = getIndexVariableRBA(model,'p',i);
     Aeq(i,idx) = mu;
     Aeq(i,1:model.noRxn) = -model.S(model.sizeXmet+model.sizeYmet+i,:);
 end

end
function [Aeq,beq] = getConstraintTotalBiomassRBA(model,B)
% getConstraintTotalBiomassRBA returns, given a deFBA/RBA model structure, 
% the equality constraint matrix and its corresponding right hand side 
% for imposing that the sum of all biomass components should be equal to a
% predefined value B.
%
% INPUT:
% model             deFBA model structure with the fields:
%   rxns                        cell array of all reaction IDs
%   rxnNames                    cell array of all metabolite names
%   rev                         0-1 array describing if the reactions are reversible (1) or irreversible (0)
%   spontaneousRxn              0-1 array describing if the reactions are spontaneous (1) or enzyme-catalyzed (0)
%   noRxn                       number or all reactions
%   mets                        cell array of all species IDs
%   metNames                    cell array of all species names
%   S                           stoichiometric matrix
%   epsilon                     scaling factor for improved numerics (see Waldherr et al 2015)
%   sizeXmet                    number of internal metabolites
%   sizeYmet                    number of external and storage metabolites
%   sizeQuotaMet                number of quota components
%   sizePmet                    number of quota and enzyme species
%   sizeXrxn                    number of reactions producing internal metabolites
%   sizeYrxn                    number of reactions involving storage or external metabolites
%   sizeQuotaRxn                number of reactions producing quota metabolites
%   sizePrxn                    number of reactions producing quota or enzymes
%   noStorage                   number of storage metabolites
%   storageWeight               array with molecular weight of storage metabolites in kDa
%   enz                         cell array with species IDs of enzymes
%   proteinWeights              array storing the molecular weights of the enzymes in kDa
%   rxnEnzRules                 length(rxns)xlength(enz) 0-1 matrix describing if enzyme j catalysises reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   Kcat_f                      length(rxns)xlength(enz) matrix storing the kcat of enzyme j for the forward direction of reaction i (entry at row i, column j is kcat)
%   initialBiomass              initial values for the quota and enzymes amounts
%   maintenanceID               ID of the maintenance reaction
%   maintenanceValue            maintenance dependence on total biomass
% B                     Predefined total biomass amount in g.
%OUTPUT
% Aeq         matrix containing total biomass equality constraint for deFBA/RBA: rows = constraints, columns = LP variables
% beq         vector containing the right hand side of the total biomass equality constraint in Aeq, s.t. Aeq*x=beq, where x is the solution vector

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

cons = emptySolutionStructRBA(model);    
Aeq = sparse(1,size(toVectorRBA(cons,model),2));
beq = B;
ind = getIndexVariableRBA(model,'p',1:model.sizeQuotaMet);
Aeq(1,ind) = 1;
ind = getIndexVariableRBA(model,'p',model.sizeQuotaMet+1:model.sizePmet);
Aeq(1,ind) = model.proteinWeights;
ind = getIndexVariableRBA(model,'s',1);
Aeq(1,ind) = model.storageWeight;
end
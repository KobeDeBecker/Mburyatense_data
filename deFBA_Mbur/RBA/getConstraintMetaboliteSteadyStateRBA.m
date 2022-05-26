function [Aeq,beq] = getConstraintMetaboliteSteadyStateRBA(model)
% getConstraintMetaboliteSteadyStateRBA returns, given a deFBA/RBA model structure, 
% the equality constraint matrix and its corresponding right hand side 
% for imposing that internal metabolites are at steady
% state
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
%OUTPUT
% Aeq         matrix containing metabolite steady state equality constraint for deFBA: rows = constraints, columns = LP variables
% beq         vector containing the right hand side of the metabolite steady state equality constraint in Aeq, s.t. Aeq*x=beq, where x is the solution vector

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

    sizeXmet = size(model.S,1)-model.sizeYmet-model.sizePmet;
    sizeXrxn = size(model.S,2)-model.sizeYrxn-model.sizePrxn;
    assert(sizeXmet == model.sizeXmet)
    assert(sizeXrxn == model.sizeXrxn) 
    
    cons = emptySolutionStructRBA(model);
    Aeq = sparse(model.sizeXmet,size(toVectorRBA(cons,model),2));
    for i = 1:sizeXmet
        % compute constraint coefficient for vx and vy
        idx = getIndexVariableRBA(model,'v',1:(sizeXrxn+model.sizeYrxn));
        Aeq(i,idx) = model.S(i,1:(sizeXrxn+model.sizeYrxn));
        % compute constraint coefficient for vp for quota rxns
        idx = getIndexVariableRBA(model,'v',(sizeXrxn+model.sizeYrxn+1):(sizeXrxn+model.sizeYrxn+model.sizeQuotaRxn));
        Aeq(i,idx) = model.S(i,(sizeXrxn+model.sizeYrxn+1):(sizeXrxn+model.sizeYrxn+model.sizeQuotaRxn));
        % for enzyme production rxns
        idx = getIndexVariableRBA(model,'v',(sizeXrxn+model.sizeYrxn+model.sizeQuotaRxn+1):(sizeXrxn+model.sizeYrxn+model.sizePrxn));
        Aeq(i,idx) = model.S(i,(sizeXrxn+model.sizeYrxn+model.sizeQuotaRxn+1):end); %model.epsilon*...
    end
    beq = zeros(1,model.sizeXmet);
end

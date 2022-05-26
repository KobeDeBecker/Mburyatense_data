function [Aineq,bineq] = getConstraintBiomassCompositionIncludingEnzymesRBA(model)
% getConstraintBiomassCompositionRBA returns, given a deFBA/RBA model structure, 
% the inequality constraint matrix and its corresponding right hand side 
% for imposing that each quota metabolite and possibly enzyme has to make up
% a certain percentage of the total biomass
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
% Aineq         matrix containing biomass composition inequality constraint for deFBA: rows = constraints, columns = LP variables
% bineq         vector containing the right hand side of the biomass composition inequality constraint in Aineq, s.t. Aineq*x<=bineq, where x is the solution vector

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

    cons = emptySolutionStructRBA(model);    
    Aineq = sparse(length(model.quotaInitial),size(toVectorRBA(cons,model),2));
    bineq = zeros(1,length(model.quotaInitial));
    
    if (size(model.proteinWeights,1)==1)
        model.proteinWeights = model.proteinWeights';
    end
    
    if (size(model.quotaWeights,1)==1)
        model.quotaWeights = model.quotaWeights';
    end
    
    weights = [ones(model.sizeQuotaMet,1);model.proteinWeights];% [model.quotaWeights;model.proteinWeights];
    assert(length(weights)==length(model.quotaInitial));
    
    for j=1:length(model.quotaInitial)
            if (model.quotaInitial(j)~=0)
                    idx = getIndexVariableRBA(model,'p',j);
                    Aineq(j,idx) = weights(j)*(model.quotaInitial(j)-1);
            
                    qidx = [1:1:j-1, j+1:1:length(model.quotaInitial)];
                    idx = getIndexVariableRBA(model,'p',qidx);
                    Aineq(j,idx) = weights(qidx)*model.quotaInitial(j);
                    
                    idx = getIndexVariableRBA(model,'s',1:model.noStorage);
                    Aineq(j,idx) = model.storageWeight*model.quotaInitial(j);

%                     idx = getIndexVariableRBA(model,'p',model.sizeQuotaMet+1:1:model.sizePmet);
%                     Aineq(j,idx) = Aineq(j,idx)*model.epsilon;

%                     idx = getIndexVariable(model,'y',i,1:1:model.noStorage);
%                     if ~isempty(idx)
%                         Aineq(t,idx) = model.storageWeight*model.quotaInitial(j);
%                     end
            end
    end

end
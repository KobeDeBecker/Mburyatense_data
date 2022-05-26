function [lb,ub] = getBoundsRBA(model)
% getBounds returns the lower and upper variable bounds needed to run RBA given a deFBA/RBA model structure
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
% lb            vector containing all variable lower bounds
% ub            vector containing all variable upper bounds

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

    cons = emptySolutionStructRBA(model);
    lb = -inf*ones(1,size(toVectorRBA(cons,model),2));
    ub = inf*ones(1,size(toVectorRBA(cons,model),2));
    assert(length(toVectorRBA(emptySolutionStructRBA(model),model)) == length(lb))
    assert(length(lb) == length(ub))
    
    % positive amounts p and s
    lb(getIndexVariableRBA(model,'p',1:1:model.sizePmet)) = 0;
    lb(getIndexVariableRBA(model,'s',1:1:model.noStorage)) = 0;
    
    % helper variables for reversible reactions should be positive
    lb(getIndexVariableRBA(model,'vRev',1:1:sum(model.rev))) = 0;
    
    % reaction irreversibility
    for k=1:model.noRxn
        if (model.rev(k)==0)
            lb(getIndexVariableRBA(model,'v',k)) = 0;
        end
    end
    
%     % substrates not present can not be taken up
%     for i = 1:length(model.Y0)
%         if model.Y0(i) == 0
%             indx = find(model.S(model.sizeXmet+i,:));
%             for j=1:length(indx)
%                 if model.rev(j) == 1
%                     if model.S(model.sizeXmet+i,:) > 0
%                         lb(j) = 0;
%                         ub(j) = Inf;
%                     else
%                         lb(j) = -Inf;
%                         ub(j) = 0;
%                     end
%                 else
%                     lb(j) = 0;
%                     ub(j) = 0;
%                 end
%             end
%         end
%     end
    
    % impose additional reaction bounds if they exist
%     if isfield(model,'rxnExtraBounds')
%         % impose customize lower bounds for uptake  
%         for i=1:length(model.rxnExtraBounds)
%             idxTime = find(model.rxnExtraBounds{i}.timePoint<=model.N);
%             if ~isempty(idxTime)
%                 idx = getIndexVariable(model,'v',model.rxnExtraBounds{i}.timePoint(idxTime),findRxnIDs(model,model.rxnExtraBounds{i}.rxnID));
%                 lb(idx) = model.rxnExtraBounds{i}.lb;
%                 ub(idx) = model.rxnExtraBounds{i}.ub;
%             end
%         end
%     end
%         
%     % impose switches in nutrients if they exist
%     if isfield(model,'switches')
%         for i=1:length(model.switches)
%             if model.switches(i).n<=model.N
%                 lb(getIndexVariable(model,'y',model.switches(i).n,model.switches(i).idx)) = model.switches(i).amount;
%                 ub(getIndexVariable(model,'y',model.switches(i).n,model.switches(i).idx)) = model.switches(i).amount;
%             end
%         end
%     end
end
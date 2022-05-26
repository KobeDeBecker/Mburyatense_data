function [Aineq,bineq] = getConstraintEnzymeCapacityRBA(model)
% getConstraintEnzymeCapacityRBA returns, given a deFBA/RBA model structure, 
% the inequality constraint matrix and its corresponding right hand side 
% for imposing that each reaction flux is bound by its
% catalyzing enzyme and the kcat
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
% Aineq         matrix containing enzyme capacity inequality constraint for deFBA: rows = constraints, columns = LP variables
% bineq         vector containing the right hand side of the enzyme capacity inequality constraint in Aineq, s.t. Aineq*x<=bineq, where x is the solution vector


% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

    sizeXYQ = size(model.S,2)-model.sizePrxn + model.sizeQuotaRxn;
    assert(sizeXYQ == model.sizeXrxn + model.sizeYrxn + model.sizeQuotaRxn)
    revRxn = model.rxns(logical(model.rev));
    
    cons = emptySolutionStructRBA(model);
    Aineq = sparse(length(model.enz),size(toVectorRBA(cons,model),2));
    bineq = zeros(1,length(model.enz));
    Aineq2 = sparse(2*length(revRxn),size(toVectorRBA(cons,model),2));
    bineq2 = zeros(1,2*length(revRxn));
    d = 1;
    
    for  i=1:length(model.enz)
        % find all rxns calalysed by this enzyme
        idxRxns = find(model.rxnEnzRules(:,i));

        % get enzyme indices
        idxp = getIndexVariableRBA(model,'p',i+model.sizeQuotaMet); % i-1 in original code
        Aineq(i,idxp) = -1;

        for j=1:length(idxRxns)
            %if we have a reversible reaction then we have to
            %use the helpers
            % assumption: kcatforward = kcatreverse
            if model.rev(idxRxns(j))==1
                idxHelp = find(ismember(revRxn,model.rxns(idxRxns(j))));
                idxHelp = getIndexVariableRBA(model,'vRev',idxHelp);
                Aineq(i,idxHelp) = 1./(model.Kcat_f(idxRxns(j),i)); % ... *model.epsilon

                % make sure helpers bound absolute value of v
                Aineq2(d,idxHelp) = -1;
                Aineq2(d,getIndexVariableRBA(model,'v',idxRxns(j))) = -1;
                d = d+1;
                Aineq2(d,idxHelp) = -1;
                Aineq2(d,getIndexVariableRBA(model,'v',idxRxns(j))) = 1;
                d = d+1;
            else
                idxv = getIndexVariableRBA(model,'v',idxRxns(j));
                Aineq(i,idxv) = 1./(model.Kcat_f(idxRxns(j),i)); % ... *model.epsilon
            end
        end              
    end
            
    Aineq = [Aineq;Aineq2];
    bineq = [bineq,bineq2];
end

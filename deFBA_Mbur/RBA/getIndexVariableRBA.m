function idx = getIndexVariableRBA(model,field,k)
% getIndexVariable computes index for a given variable in the RBA solution 
% vector given its name and entry
%     
% solutions have the following variables
% v    -  dimensions model.noRxn x 1
% p    -  dimensions model.sizePmet x 1
% vRev -  dimensions sum(model.rev) x 1
%
% For all variables slicing is supported
%
% The order of variables is
% v, p, vRev
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
% field             String containing the type of variable (options: 'v', 'p', 'vRev')
% k                 First dimension index (entry)
% 
%OUTPUT
% idx               The variable index in the solution vector.

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 
     
    % we have one specific growth rate variable
    % we have a flux variable for each reaction
    sizeV = model.noRxn;
    % we have a macromolecule amount at for each macromolecule
    sizeP = model.sizePmet;
    % we have a reversible flux value for each reversible reaction
    sizeVrev = sum(model.rev);
    
    % for each flux variable of a reversible reaction we have in addition a helper variable vRev that helps us get the absolute value

    % the order of variables is
    % v, p, vRev, s
    
    start = 0;
    
    if strcmp(field,'v')
        idx = start + k;
    end
    
    if strcmp(field,'p')
        idx = start + sizeV + k;
    end 
      
    if strcmp(field,'vRev')
        idx = start + sizeV + sizeP + k;
    end
    
    if strcmp(field,'s') % Storage compound
        idx = start + sizeV + sizeP + sizeVrev + k;
    end
end

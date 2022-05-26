function w = toVectorRBA(str, model)
%toVectorRBA transform a struct of strtraint coefficients or solution struct into vector
%
% w = toVectorRBA(str,model)
%
% solutions have the following variables
% v    -  dimensions 1 x model.noRxn
% p    -  dimensions 1 x model.sizePmet
% vRev -  dimensions model.N x sum(model.rev)
%
% the order of variables is
% v, p, vRev
%
%INPUTS
% str        Struct that contains the coefficients for the constraints or the solution entries.
%            The struct should have the fields: v, p, vRev
% model             deFBA/RBA model structure with the fields:
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
%
%OUTPUT
% w           vector version of the coefficients or of the solution

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 
  
    % check structure sizes for consistency
    if (size(str.v,1)~=1 || size(str.v,2)~=model.noRxn ||...
            size(str.p,1)~=1 || size(str.p,2)~=model.sizePmet)
        
        msgID = 'toVector:BadSize';
        msg = 'Structure sizes incompatible with the model';
        baseException = MException(msgID,msg);
        throw(baseException)
    end

    v = reshape(str.v,1,numel(str.v));
    p = reshape(str.p,1,numel(str.p));
    vRev = reshape(str.vRev,1,numel(str.vRev));
    s = reshape(str.s,1,numel(str.s));
    w = sparse([v,p,vRev,s]);    
end

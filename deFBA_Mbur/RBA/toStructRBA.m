%toStructRBA transform a vector of solutions or constraint coefficients into a
%         solution struct. The method is the exact inverse of the toVectorRBA function
%
% str = toStructRBA(w,model)
%
%INPUT:
% w                 Vector to be converted to structure
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
% str         Structure with the fields:
%               v = fluxes at all middles of discretization intervals
%               p = macromolecule amounts
%               vRev = helper variable for bounding of reversible reaction fluxes

%   v(internal reactions) - Unit: mmol/*h
%   v(exchange reactions) - Unit: mmol/*h
%   v(quota reactions) 	  - Unit: 1/*h
%   v(enzyme synthesis)   - Unit: mmol/h
%   p(quota)              - Unit: pg
%   p(enzymes)            - Unit: mmol

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function str = toStructRBA(w,model)
% solutions have the following variables
% v    -  dimensions 1 x model.noRxn
% p    -  dimensions 1 x model.sizePmet
% vRev -  dimensions 1 x sum(model.rev)
%
% the order of variables is
% v, p, vRev

    w=full(w);  
    % we have a flux variable for each reaction
    sizeV = model.noRxn;
    % we have a macromolecule amount for each macromolecule
    sizeP = model.sizePmet;
    % for each flux variable of a reversible reaction we have in addition a helper variable vRev that helps us get the absolute value
    sizeVrev = sum(model.rev);
    
    str.v = w(1:sizeV);
    start = sizeV;
    str.p = w(start+1:start+sizeP);
    start = start + sizeP;
    str.vRev = w(start+1:start+sizeVrev);
    start = start + sizeVrev;
    
    assert(length(w) == start + model.noStorage)
    str.s = w(start+1:end);
end

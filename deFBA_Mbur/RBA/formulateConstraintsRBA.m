function [Aeq,beq,Aineq,bineq,lb,ub] = formulateConstraintsRBA(model,mu,B)
% formulateConstraintsRBA returns the equality and inequality constraint matrices and
% right hand side vectors, and the lower and upper bounds needed to run
% RBA given a deFBA/RBA model structure
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
% OUTPUT:
%  Aeq           matrix containing all equality constraints for deFBA: rows = constraints, columns = LP variables
%  beq           vector containing the right hand side of all equality constraints in Aeq
%  Aineq         matrix containing all inequality constraints for deFBA: rows = constraints, columns = LP variables
%  bineq         vector containing the right hand side of all inequality constraints in Aineq, s.t. Aineq*x<=bineq, where x is the solution vector
%  lb            vector containing all variable lower bounds
%  ub            vector containing all variable upper bounds

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

    [Aeq1,beq1] = getConstraintBalancedGrowthRBA(model,mu);
    %disp('formulated muP = Sv') 
    
    [Aeq2,beq2] = getConstraintMetaboliteSteadyStateRBA(model);
    %disp('formulated constraint steady state')
    
    [Aeq3,beq3] = getConstraintTotalBiomassRBA(model,B);
    %disp('formulated total biomass constraint')
    
    [Aineq1,bineq1] = getConstraintEnzymeCapacityRBA(model);
    %disp('formulated constraint kcat')

    [Aineq2,bineq2] = getConstraintBiomassCompositionIncludingEnzymesRBA(model);
    %disp('formulated biomass composition constraint')
    
    [Aineq3,bineq3] = getConstraintMaintenanceRBA(model);
    %disp('formulated maintenance constraint')
    if isempty(find(Aineq3,1))
        Aineq3 = [];
        bineq3 = [];
    end
    
    [lb,ub] = getBoundsRBA(model);
    %disp('formulated variable bounds')
    
    Aeq = [Aeq1;Aeq2;Aeq3];
    beq = [beq1,beq2,beq3];
     
    Aineq = [Aineq1;Aineq2;Aineq3];
    bineq  = [bineq1,bineq2,bineq3];
end
function Out = deFBA(model,param,output,rib,seq)
% This function solves the deFBA model specified by the input "model" and
% returns the output "Out" as specified by the output vector "output".
% INPUTS
%  - model             deFBA model structure with the fields:
%   rxns                        cell array of all reaction IDs
%   rxnNames                    cell array of all metabolite names
%   rev                         0-1 array describing if the reactions are reversible (1) or irreversible (0)
%   spontaneousRxn              0-1 array describing if the reactions are spontaneous (1) or enzyme-catalyzed (0)
%   noRxn                       number or all reactions
%   grRules                     cell array with strings describing the gene association for each reaction
%   rxnECNumbers                cell array containing the EC number for each reaction
%   mets                        cell array of all species IDs
%   metNames                    cell array of all species names
%   S                           stoichiometric matrix
%   genes                       cell array of all gene names
%   rxnGeneMat                  length(rxns)xlength(genes) 0-1 matrix describing if gene j is involved in the catalysis of reaction i (entry at row i, column j is 1) or not (entry at row i, column j is 0)
%   tf                          final time of simulation (hours)
%   h                           length of each time interval
%   N                           number of discretization points for deFBA timecourses (N = round(tf/h))
%   epsilon                     scaling factor for improved numerics (see Waldherr et al 2015)
%   beta                        dilution factor for enzymes, quota and storage
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
%   quotaInitial                array containing the total biomass percentages that have to be satisfied for each quota compound at each time point
%   initialBiomass              initial values for the quota and enzymes amounts
%   Y0                          initial values for the storage (first noStorage entries) and external metabolites (rest) amounts
%   maintenanceID               ID of the maintenance reaction
%   maintenanceValue            maintenance dependence on total biomass
%   gprComp                     cell array that specifies for each reaction the recipe for building its corresponding enzyme from the individual gene ids ( e.g. 3*gene1 AND 2*gene2 means that the enzyme is made of three copies of gene1 and 2 copies of gene2)
%   continuousMet               vector indicating whether the extracellular metabolites are continuously entering and leaving the reactor
%   volume                      the consntant liquid volume of the reactor
%   saturation                  the saturation concentrations of the continuous extravellular metabolites, otherwise zero
% - param              A vector containing the catalytic constants of all the enzymes in the same order as the enzymes in model.enz. The last parameter is the ribosome constant
% - output             An array indicating the wanted output of the model
% - rib                1 or 0, if rib == 1 the last entries in param are the catalytic constants for the riobosome, if rib == 0 the last entry is the amount of amino acidds added per unit of time.
% - seq                if rib == 0, this array contains the amino acid sequences of the enzymes.
%
% OUTPUTS
% - Out                    The wanted output as specified by output

% Written by Kobe De Becker
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

%% Generating Kcat_f
% Kcat_f is a matrix that indicates which the catalytic constants apply
% for which enzymes and which reactions. It has a same number as rows as
% there aren reactions and the number of columns equals the number of
% enzymes (+ the ribosome)

if rib == 0
% Add kcat for ribosome --> Based on length of amlino acids! kcat  = a/l
% with a = number of aa per unit of time and l the length of the
% aa-sequence. a = 12 (slow growing, a = 17 (fast growing) in E. coli
    pl = length(param);
    a = param(end);
    for i=1:length(model.enz)
        param{pl-1+i} = a/length(seq{i});
    end
    K = model.rxnEnzRules;
    for i = 1:size(K,2)
        ind = find(K(:,i));
        if i ~= size(K,2)
            K(ind,i) = param(i);
        else
            for t = 1:length(ind)
                K(ind(t),i) = param(pl-1+t);
            end
        end
    end
elseif rib == 1
    K = double(model.rxnEnzRules);
    for i = 1:size(K,2)
        ind = find(K(:,i));
        if i ~= size(K,2)
            K(ind,i) = param(i);
        else
            for t = 1:length(ind)
                K(ind(t),i) = param(length(model.enz)+t-1);
            end
        end
    end
else
    disp('Invalid number for rib')
end

model.Kcat_f = K;

%% Solve the deFBA model
r = run_deFBA(model,'cplex');
s = toStruct(r,model);

%% Extract wanted output
if strcmp(output,'b')
    Out = sum(s.p(end,1:model.sizeQuotaMet)) + model.proteinWeights'*s.p(end,model.sizeQuotaMet+1:end)';
elseif strcmp(output(1),'y')
    ind = str2double(output(2:end));
    Out = s.y(end,ind);
elseif strcmp(output(1),'p')
    ind = str2double(output(2:end));
    Out = s.p(end,ind);
else
    disp('Invalid output input')
end
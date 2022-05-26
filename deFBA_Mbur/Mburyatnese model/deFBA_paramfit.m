% This function executes the deFBA DOA approach specifically to obtain
% parameter estimates as decided by the sensitivity anlysis and
% corresponding literature study. The function returns the total biomass
% amount and DO-percentage at predescribed time points. If needed, linear
% interpolation between the discretization points of the deFBA model will
% be used.

% Written by Kobe De Becker - 3/2/2021
% This file is covered by the GNU GENERAL PUBLIC LICENSE in terms of 
% copyright, referencing and distribution. 

function [X,DO,OC] = deFBA_paramfit(p,tpX,tpDO,model,tf,N)

% Define the parametervalues in the model

% p(1) = kcat_pMMO
indx_col = strcmp(model.enz,'Methane monooxygenase');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(1);

% p(2) = kcat_acpSmal
indx_col = strcmp(model.enz,'(acyl-carrier-protein) S-malonyltransferase');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(2);

% p(3) = kcat_G3P_deh
indx_col = strcmp(model.enz,'glycerol-3-phosphate dehydrogenase');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(3);

% p(4) = kcat_argsucc
indx_col = strcmp(model.enz,'argininosuccinate synthase');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(4);

% p(5) = kcat_pglu_deh
indx_col = strcmp(model.enz,'phosphogluconate dehydratase');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(5);

% p(6) = kcat_PlsX
indx_col = strcmp(model.enz,'phosphate acetyltransferase PlsX');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(6);

% p(7) = kcat_PlsY
indx_col = strcmp(model.enz,'PlsY');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(7);

% p(8) = kcat_PlsC
indx_col = strcmp(model.enz,'PlsC');
indx_row = model.Kcat_f(:,indx_col) ~= 0;
model.Kcat_f(indx_row,indx_col) = p(8);

if length(p) > 10
% p(11) = kLa_O2 (and indirectly also kLa_CH4 due to a constant relationship
% between both constants (index 3 = CH4, 5 = O2)
    model.massTransfer(5) = p(11);
    model.massTransfer(3) = 0.7299*p(11); % Harsh assumption...
end

% Resdefine initial biomass via RBA (p(10) is the initial biomass amount)
[P,S] = RBA_function_it(p(10)*10^3,model,0.04,10^(-2));
P = P/10^3; % Actual total biomass amount should be rescaled
S = S/10^3;
model.initialBiomass = P';
model.Y0(1) = S;

% p(9) = w_glc objective weight for glycogen
model.objectiveWeights = [p(9),model.objectiveWeights(2:end)];

% Run deFBA
model.tf = tf;
model.N = N;
model.h = model.tf/model.N;

model.Y0(6) = 30; % Nitrate as theoretically present in  the reactor medium
model.Y0(8) = 7; % Phosphate as theoretically present in  the reactor medium
model.Y0(9) = 2.5; % Sulphate as theoretically present in  the reactor medium
model.Y0(3) = model.Y0(3)*model.volume;
model.Y0(5) = model.Y0(5)*model.volume;

r = run_deFBA(model,'cplex');
if isempty(r)
    dummy = toVector(emptySolutionStruct(model),model);
    r = zeros(size(dummy));
end
s = toStruct(r,model);

% Determine total biomass
biomass0 = p(10);
% biomass0 = biomass0 + model.proteinWeights'*s.p0(model.sizeQuotaMet+1:end)'+model.storageWeight*s.y0(1);
biomass = zeros(model.N,1);
for i = 1:size(s.p,1)
    biomass(i) = sum(s.p(i,1:model.sizeQuotaMet)) + model.proteinWeights'*s.p(i,model.sizeQuotaMet+1:end)'+model.storageWeight*s.y(i,1);
end
biomass = [biomass0;biomass];
X = interp1(0:model.h:model.tf,biomass,tpX);

% Determine dissolved oxygen (DO)
DO_deFBA = [s.y0(5); s.y(:,5)]./model.saturation(5).*100./model.volume;
DO = interp1(0:model.h:model.tf,DO_deFBA,tpDO);

% Determine Oxygen on methane uptake ratio
O2_up = s.v(:,379)+s.v(:,380);
CH4_up = s.v(:,378);
OC = O2_up./CH4_up;

end
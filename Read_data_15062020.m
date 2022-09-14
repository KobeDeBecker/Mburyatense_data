% This is a template script to load bioreactor experiment data with
% Methylomicrobium buryatense 5GB1 into MATLAB. The only required input is
% the name of the three scripts from which these data can be received.
%
% 1. Excel-file with dry-weight, offline OD, and HPLC results
% 2. Csv-file with reactor parameters and sensor measurements (BIOFLO)
% 3. Csv-file with off-gas results (BlueVis)
%
% All data are loaded as a matrix with two columns, in which the first
% column represents the measurement time and the second contains the
% measured value.

clear variables
close all
clc

%% Initialisation
% Add here the names of the respective documents
F1 = 'BiomassHPLC_20200615.xlsx';
F2 = 'CTPC06603.20200616 5GB1.Control.csv'; % Inoculation time seems 2 hours off?!
F3 = 'BlueVis_Export_beginning_from_2020_6_16_10_46.csv'; 

% Add here whether you want to show the plots in MATLAB (1) or if you only
% want to load the data (0)
graph = 1;

%% File 1 (F1)
Data1 = readmatrix(F1);
Data1(:,1:2) = [];

Time1 = Data1(:,1);

% Offline OD [-]
indx = get_index_nan(Data1(:,2));
OD_offline = Data1(:,1:2);
OD_offline(indx,:) = [];
OD_offline(:,2) = OD_offline(:,2) - OD_offline(1,2);

% CDW [g/L]
indx = get_index_nan(Data1(:,3));
CDW = Data1(:,[1,3]);
CDW(indx,:) = [];

% Formic acid [mg/L]
indx = get_index_nan(Data1(:,7));
Form = Data1(:,[1,7]);
Form(indx,:) = [];

% Acetic acid [mg/L]
indx = get_index_nan(Data1(:,8));
Ac = Data1(:,[1,8]);
Ac(indx,:) = [];

% Plot results
if graph % Comment those graphs that are not needed
    figure('Name','Biomass')
    yyaxis left
    plot(OD_offline(:,1),OD_offline(:,2),'bo','MarkerSize',2,'LineWidth',1.5)
    ylabel('OD_{offline} [-]','FontSize',16)
    yyaxis right
    plot(CDW(:,1),CDW(:,2),'ro','MarkerSize',2,'LineWidth',1.5)
    ylabel('CDW [g/L]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    
    figure('Name','Organic acids')
    yyaxis left
    plot(Ac(:,1),Ac(:,2),'bo','MarkerSize',2,'LineWidth',1.5)
    ylabel('Ac [mg/L]','FontSize',16)
    yyaxis right
    plot(Form(:,1),Form(:,2),'ro','MarkerSize',2,'LineWidth',1.5)
    ylabel('Form [mg/L]','FontSize',16)
    xlabel('t [h]','FontSize',16)
end

%% File 2 (F2)
% Always check order and number of rows/columns when loading data! Apparantly,
% this can change during experiments!
Data2 = readmatrix(F2,'OutputType','string');
Data2(1:80,:) = []; % Removing lines on top
Data2(8731:end,:) = []; % removing lines on bottom
% indx = find(~ismissing(Data2(:,3)),1); % The inoculation time is mentioned in the third column
% Inoculation_Time = datestr(datenum(Data2(indx,1)) - hour(Data2(indx,3))/24 - minute(Data2(indx,3))/(60*24) - second(Data2(indx,3))/(60^2*24),'dd-mm-yyyy HH:MM:SS');
Inoculation_Time = '2020-06-16 12:00:00';
formatIn = 'yyyy-mm-dd HH:MM:SS';
indx = find(datenum(Data2(:,1)) > datenum(Inoculation_Time),1);
Data2_Time = datenum(Data2(:,1),formatIn); % The first column contains the time
Data2_Time = (Data2_Time - datenum(Inoculation_Time,formatIn))*24; % Convert units from days to hours
Data2(1:indx-1,:) = []; % Remove all data from before inoculation
Data2_Time(1:indx-1,:) = []; % Remove all data from before inoculation
Data2(:,2:3) = []; % Remove duration and inoculation data from Data2 (always deduct 2 from column number!)
Data2_Val = str2double(Data2(:,2:end)); % Convert string values to numbers with double precision
Data2 = [Data2_Time Data2_Val];

if graph
    % Variables of interest (others are also present in the csv file and
    % can be accessed in a similar way):
    %   DO1.PV --> Column 7 --> 5
    %   ExternalA1.MFM.FPV (off gas flow rate process value) --> Column 9 --> 7
    %   ExternalB1.MFC.FPV (flow rate methane flow) --> Column 10 --> 8
    %   ExternalB1.MFC.FSPM (manual set point for methane flow rate) --> Column 12 --> 10
    %   Unit 1.F1.PV (total gas flow rate (excluding methane)) --> Column 13 --> 11
    %   Unit 1.FAir1.PV --> Column 17 --> 15
    %   Unit 1.FAir1.SP --> Column 18 --> 16
    %   Unit 1.FN21.PV --> Column 31 --> 29
    %   Unit 1.FN21.SP --> Column 32 --> 30
    %   Unit 1.N1.PV (stirrer speed) --> Column 35 --> 33
    %   Unit 1.ODAU1.PV --> Column 37 --> 35
    %   Unit 1.ODCX1.PV --> Column 38 --> 36
    %   Unit 1.pH1.PV --> Column 44 --> 42
    %   Unit 1.VA1.PV (Added acid volume) --> Column 54 --> 52 
    %   Unit 1.VB1.PV (Added base volume) --> Column 56 --> 54
    
    figure('Name','DO')
    plot(Data2(:,1),Data2(:,5),'bo','LineWidth',1.5,'MarkerSize',2)
    ylabel('DO [%]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    xlim([0 80])
    ylim([0 100])
    
    figure('Name','Gas flow rates')
    yyaxis left
    plot(Data2(:,1),Data2(:,7),'bo','MarkerSize',2,'LineWidth',1.5)
    ylim([0 20])
    ylabel('F_{Total} [sL/h]','FontSize',16)
    yyaxis right
    plot(Data2(:,1),Data2(:,8),'ro','MarkerSize',2,'LineWidth',1.5)
    hold on
    plot(Data2(:,1),Data2(:,10),'go','MarkerSize',2,'LineWidth',1.5)
    plot(Data2(:,1),Data2(:,11),'co','MarkerSize',2,'LineWidth',1.5)
    legend({'F_{Total,uit}','F_{CH_4} PV','F_{CH_4} SP','F_{Air} PV'},'FontSize',12,'Location','EastOutside')
    ylim([0 20])
    ylabel('F_{CH_4} [sL/h]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    
    figure('Name','Stirrer speed')
    plot(Data2(:,1),Data2(:,33),'bo','LineWidth',1.5,'MarkerSize',2)
    ylabel('N [rpm]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    ylim([0 1500])
    
    figure('Name','OD')
    % yyaxis left
    plot(Data2(:,1),Data2(:,35),'bo','MarkerSize',2,'LineWidth',1.5)
    ylabel('ODAU [-]','FontSize',16)
%     yyaxis right
%     plot(Data2(:,1),Data2(:,36),'ro','MarkerSize',2,'LineWidth',1.5)
%     ylabel('ODCX [-]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    
    figure('Name','pH')
    yyaxis left
    plot(Data2(:,1),Data2(:,42),'bo','MarkerSize',2,'LineWidth',1.5)
    ylabel('pH [-]','FontSize',16)
    yyaxis right
    plot(Data2(:,1),Data2(:,52),'ro','MarkerSize',2,'LineWidth',1.5)
    hold on
    plot(Data2(:,1),Data2(:,54),'go','MarkerSize',2,'LineWidth',1.5)
    ylabel('V_{Added} [mL]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    legend({'pH','Acid','Base'},'FontSize',12,'Location','NorthEast')
end

%% File 3 (F3)
% Always check order and number of rows/columns when loading data! Apparantly,
% this can change during experiments!
Data3 = readmatrix(F3,'OutputType','string');
Data3(:,39:end) = [];
Data3(:,[1:9,12:3:38,25,26]) = [];
Data3(1:4,:) = [];
indx = 1:size(Data3,2);
even = indx(mod(indx,2) == 0); % Variable values
Data3_Var = str2double(strrep(Data3(:,even),',','.'));
odd = indx(mod(indx,2) ~= 0);
formatIn = 'dd/mm/yyyy HH:MM:SS';
Data3_Time = Data3(:,odd);
Data3_Time = Data3_Time(:);
Data3_Time = datenum(Data3_Time,formatIn);
Data3_Time = Data3_Time - datenum(Inoculation_Time,'yyyy-mm-dd HH:MM:SS'); % Check for inoculation time...??
Data3_Time = reshape(Data3_Time,size(Data3(:,odd))).*24; % There are 24 hours in one day
Data3 = zeros(size(Data3));
Data3(:,even) = Data3_Var;
Data3(:,odd) = Data3_Time;

if graph % Comment those graphs that are not needed
    % Humidity plots (absolute and relative)
    figure('Name','Humidity')
    yyaxis left
    plot(Data3(:,11),Data3(:,12),'-b','LineWidth',2) % Measured by Bluevary
%     hold on
%     plot(Data2(:,13),Data2(:,14),'-k','LineWidth',2) % Measured by BCP-CH4
    ylabel('Absolute humidity [vol%]','FontSize',16)
    yyaxis right
    plot(Data3(:,13),Data3(:,14),'-r','LineWidth',2)
    ylabel('Relative humidity [vol%]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    xlim([0 80])

    % Temperature plot
    figure('Name','Temperature off-gas')
    plot(Data3(:,15),Data3(:,16),'-b','LineWidth',2)
    ylabel('Temperature [°C]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    xlim([0 80])
    
    % Pressure plot
    figure('Name','Pressure off-gas')
    plot(Data3(:,17),Data3(:,18),'-b','LineWidth',2)
    hold on
    plot(Data3(:,3),Data3(:,4),'-r','LineWidth',2)
    ylabel('Pressure [bar]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    legend({'Sensor 1 (Bluevary)','Sensor 2 (BCP-CH4)'},'Location','NorthEast','FontSize',12)
    xlim([0 80])
    
    % Off-gas composition
    figure('Name','Off-gas composition')
    plot(Data3(:,11),Data3(:,12),'-','LineWidth',2)
    hold on
    plot(Data3(:,5),Data3(:,6),'-','LineWidth',2)
    plot(Data3(:,7),Data3(:,8),'-','LineWidth',2)
    plot(Data3(:,1),Data3(:,2),'-','LineWidth',2)
    ylabel('Amount [vol%]','FontSize',16)
    xlabel('t [h]','FontSize',16)
    legend({'H_2O','CO_2','O_2','CH_4'},'Location','EastOutside','FontSize',12)
    xlim([0 80])
end

%% Calculate oxygen transfer rate (OTR), methane transfer rate (MTR)and carbon dioxide transfer rate (CTR)(assuming ideal gas)

F_in = [Data2(:,1), Data2(:,8) + Data2(:,11)];
F_in(isnan(F_in(:,2)),:) = [];

X_CH4_in = [Data2(:,1),Data2(:,8).*0.15./(Data2(:,8) + Data2(:,11))];
X_O2_in = [Data2(:,1),Data2(:,11).*0.21./(Data2(:,8) + Data2(:,11))];
X_CH4_in(isnan(X_CH4_in(:,2)),:) = [];
X_O2_in(isnan(X_O2_in(:,2)),:) = [];
X_H2O_in = 0;
X_CO2_in = 0.0022;

X_CH4_out = [Data3(:,1) Data3(:,2)];
X_O2_out = [Data3(:,7) Data3(:,8)];
X_H2O_out = [Data3(:,11) Data3(:,12)];
X_CO2_out = [Data3(:,5) Data3(:,6)];

X_CH4_out = [F_in(:,1),interp1([0; X_CH4_out(:,1)],[X_CH4_out(1,2); X_CH4_out(:,2)],F_in(:,1))]; 
X_O2_out = [F_in(:,1), interp1([0; X_O2_out(:,1)],[X_O2_out(1,2); X_O2_out(:,2)],F_in(:,1))]; 
X_H2O_out = [F_in(:,1), interp1([0; X_H2O_out(:,1)],[X_H2O_out(1,2); X_H2O_out(:,2)],F_in(:,1))]; 
X_CO2_out = [F_in(:,1), interp1([0; X_CO2_out(:,1)],[X_CO2_out(1,2); X_CO2_out(:,2)],F_in(:,1))]; 

% X_CH4_out(isnan(X_CH4_out(:,2)),:) = [];
% X_O2_out(isnan(X_O2_out(:,2)),:) = [];
% X_H2O_out(isnan(X_H2O_out(:,2)),:) = [];
% X_CO2_out(isnan(X_CO2_out(:,2)),:) = [];

F_out = [F_in(:,1), F_in(:,2).*(1-X_CH4_in(:,2)-X_O2_in(:,2)-X_H2O_in-X_CO2_in)./((100-X_CH4_out(:,2)-X_O2_out(:,2)-X_H2O_out(:,2)-X_CO2_out(:,2))/100)]; % Based on nitrogen balance

V = [Data2(:,1),Data2(:,51)];
V(isnan(V(:,2)),:) = [];
V = [F_in(:,1), interp1(V(:,1),V(:,2),F_in(:,1)).*10^(-3)]; % Convert to [L] instead of [mL]

Ts = 25 + 273.15; % [K] Standard temperature
Ps = 1.01325*10^5; % [Pa] Standard pressure
R = 8.314; % [J/molK] Gas constant
c = Ps/(R*Ts)*10^(-3); % [mol/L] Molar volume 

OTR = [F_in(:,1),(c./V(:,2)).*(F_in(:,2).*X_O2_in(:,2)-F_out(:,2).*X_O2_out(:,2)/100)];
MTR = [F_in(:,1),(c./V(:,2)).*(F_in(:,2).*X_CH4_in(:,2)-F_out(:,2).*X_CH4_out(:,2)/100)];
CTR = [F_in(:,1),(c./V(:,2)).*(-F_in(:,2).*X_CO2_in+F_out(:,2).*X_CO2_out(:,2)/100)];

CO2_produced = [F_in(:,1),c.*(-F_in(:,2).*X_CO2_in+F_out(:,2).*X_CO2_out(:,2)/100)];

if graph
    figure('Name','Transfer rates')
    indx = [];
    for i = 1:size(OTR,1)
        if OTR(i,1) < 5.699
            indx = [indx, i];
        elseif OTR(i,1) >= 15.76 && OTR(i,1) <= 15.77
            indx = [indx, i];
        end
    end
    OTR(indx,:) = [];
%     for i = 1:size(MTR,1)
%         if MTR(i,1) < 5.699
%             indx = [indx, i];
%         elseif MTR(i,1) >= 15.76 && MTR(i,1) <= 15.77
%             indx = [indx, i];
%         end
%     end
    MTR(indx,:) = [];
    CTR = [0 0; CTR];
    MTR = [0 0; MTR];
    OTR = [0 0; OTR];
    plot(OTR(:,1),OTR(:,2),'-b','LineWidth',2)
    hold on
    plot(MTR(:,1),MTR(:,2),'-r','LineWidth',2)
    for i = 1:size(CTR,1)
        if CTR(i,2) < 0
            CTR(i,2) = 0;
        end
    end
    plot(CTR(:,1),CTR(:,2),'-g','LineWidth',2)
    xlabel('t [h]','FontSize',16)
    ylabel('Transfer rates [mol/Lh]','FontSize',16)
    legend({'OTR','MTR','CTR'},'FontSize',12,'Location','NorthEast')
    ylim([0 0.005])
end

% O2/CH4 ratio
r_O2_CH4 = OTR(:,2)./MTR(:,2);
if graph
    figure
    plot(OTR(:,1),r_O2_CH4,'b','LineWidth',2)
    xlabel('t [h]','FontSize',16)
    ylabel('O_2/CH_4 ratio [-]','FontSize',16)
end

% Carbon distribution
F_MTR = [];
for i = 2:size(MTR,1)
    F_MTR_next = trapz(MTR(1:i,1),MTR(1:i,2)*3);
    F_MTR = [F_MTR; MTR(i,1) F_MTR_next]; 
end

F_CTR = [];
indx = isnan(CTR(:,2));
CTR(indx,:) = [];
for i = 2:size(CTR,1)
    F_CTR_next = trapz(CTR(1:i,1),CTR(1:i,2)*3);
    F_CTR = [F_CTR; CTR(i,1) F_CTR_next]; 
end
F_CTR = [F_MTR(:,1) interp1(F_CTR(:,1),F_CTR(:,2),F_MTR(:,1))];

phi = 41.0059*10^(-3); % 1 g_DW contains phi moles carbon atoms (from biomass reaction)
a =  0.6639; % See file Biomass_calibration
calibration = @ (x) a.*x; % Calibration curve from online OD to CDW
OD = [Data2(:,1),Data2(:,35)];
indx = isnan(OD(:,2));
OD(indx,:) = []; % remove NAN fields
OD(:,2) = OD(:,2) - min(OD(:,2)); % Correction for blanc measurement
X = [OD(:,1) calibration(OD(:,2))];
C_Biomass = [X(:,1), phi.*(X(:,2)-min(X(:,2)))*3]; % Correction for X(0)!
C_Biomass = [F_MTR(:,1), interp1(C_Biomass(:,1),C_Biomass(:,2),F_MTR(:,1))]; % Interpolate to obtain the valueas at the same time points as for the other contributions
size_X = size(X);
X_label = mat2cell([[0,0]; X], ones(size_X(1)+1,1), ones(size_X(2),1));;
X_label{1,1} = "Time [h]";
X_label{1,2} = "Calibrated OD [AU]";
writecell(X_label, "OD_15062020.csv", "FileType", "text");

frac_X = C_Biomass(:,2);%./F_MTR(:,2);
frac_CO2 = F_CTR(:,2);%./F_MTR(:,2);
frac_P = F_MTR(:,2) - frac_X - frac_CO2;% 1 - frac_X - frac_CO2;

figure
plot(F_MTR(:,1),frac_X,'-r','LineWidth',2)
hold on
plot(F_MTR(:,1),frac_CO2,'-b','LineWidth',2)
plot(F_MTR(:,1),frac_P,'-g','LineWidth',2)
plot(F_MTR(:,1),F_MTR(:,2),'-k','LineWidth',2)
xlabel('t [h]','FontSize',16)
ylabel('Cummulative moles of carbon [mol]','FontSize',16)
legend({'Biomass','CO_2','Products','Carbon in'},'FontSize',14,'Location','NorthWest')

%% Determine total carbon distribution

% C_in_CH4 = zeros(size(F_in));
% for i = 10:size(F_in,1)
%     C_in_next = trapz([0; F_in(1:i,1)],[0; F_in(1:i,2).*(X_CH4_in(1:i,2)).*c]);
%     C_in_CH4(i,:) = [F_in(i,1), C_in_next];
% end
% C_in_CH4 = [0 0;C_in_CH4]; % Add zero point
% 
% C_out_CH4 = zeros(size(F_out));
% for i = 10:size(F_out,1)
%     C_out_CH4_next = trapz([0; F_out(1:i,1)],[0; F_out(1:i,2).*X_CH4_out(1:i,2)/100.*c]);
%     C_out_CH4(i,:) = [F_out(i,1), C_out_CH4_next];
% end
% C_out_CH4 = [0 0;C_out_CH4]; % Add zero point
% 
% C_out_CO2 = zeros(size(F_out));
% for i = 10:size(F_out,1)
%     C_out_CO2_next = trapz([0; F_out(1:i,1)],[0; F_out(1:i,2).*X_CO2_out(1:i,2)/100.*c]);
%     C_out_CO2(i,:) = [F_out(i,1), C_out_CO2_next];
% end
% C_out_CO2 = [0 0;C_out_CO2]; % Add zero point
% 
% C_in_CO2 = zeros(size(F_in));
% for i = 10:size(F_in,1)
%     C_in_CO2_next = trapz([0; F_in(1:i,1)],[0; F_in(1:i,2).*X_CO2_in.*c]);
%     C_in_CO2(i,:) = [F_in(i,1), C_in_CO2_next];
% end
% C_in_CO2 = [0 0;C_in_CO2]; % Add zero point
% 
% phi = 41.0059*10^(-3); % 1 g_DW contains phi moles carbon atoms (from biomass reaction)
% a =  0.6639; % See file Biomass_calibration
% calibration = @ (x) a.*x; % Calibration curve from online OD to CDW
% OD = [Data2(:,1),Data2(:,35)];
% indx = isnan(OD(:,2));
% OD(indx,:) = []; % remove NAN fields
% OD(:,2) = OD(:,2) - min(OD(:,2)); % Correction for blanc measurement
% X = [OD(:,1) calibration(OD(:,2))];
% C_Biomass = [0, 0; X(:,1), phi.*(X(:,2)-min(X(:,2)))]; % Correction for X(0)!
% C_Biomass = [C_in_CH4(10:end,1), interp1(C_Biomass(:,1),C_Biomass(:,2),C_in_CH4(10:end,1))]; % Interpolate to obtain the valueas at the same time points as for the other contributions
% 
% % To_cells = [C_in_CH4(:,1), C_in_CH4(:,2)-C_out_CH4(:,2)-C_out_CO2(:,2)];
% % indx = To_cells(:,2) <= 0;
% % To_cells(indx,:) = []; % Remove values where there is no significant methane uptake to the liquid phase
% % lambda = To_cells(1,1); % Determine estimate for the lag time
% % indx = find(C_in(:,1) >= lambda,1);
% 
% % frac_CO2 = [C_in_CH4(:,1), (C_out_CO2(:,2)-C_in_CO2(:,2))./(C_in_CH4(:,2)-C_out_CH4(:,2)).*100];
% % frac_X = [C_in_CH4(:,1), C_Biomass(:,2)./(C_in_CH4(:,2)-C_out_CH4(:,2)).*100];
% % frac_Prod = [C_in_CH4(:,1), 100 - frac_CO2(:,2) - frac_X(:,2)];
% 
% C_in_tot = [C_in_CH4(:,1) C_in_CH4(:,2)+C_in_CO2(:,2)];
% C_out_tot = [C_out_CH4(:,1) C_out_CH4(:,2)+C_out_CO2(:,2)];
% frac_C_out = [C_out_CH4(:,1) C_out_tot(:,2)./C_in_tot(:,2)];
% frac_X = [C_out_CH4(10:end,1) C_Biomass(:,2)./C_in_tot(10:end,2)]; % From the total incoming carbon in the reactor
% 
% figure
% plot(frac_C_out(:,1),frac_C_out(:,2),'-b','LineWidth',2)
% hold on
% plot(frac_X(:,1),frac_X(:,2),'-r','LineWidth',2)
% % plot(frac_CO2(:,1),frac_CO2(:,2),'-b','LineWidth',2)
% % hold on
% % plot(frac_X(:,1),frac_X(:,2),'-r','LineWidth',2)
% % plot(frac_Prod(:,1),frac_Prod(:,2),'-g','LineWidth',2)
% % legend({'CO_2','X','P'},'Location','EastOutside')

%% Functions
function indx = get_index_nan(meas)
    indx = [];
    n = size(meas,1);
    for i = 1:n
        if isnan(meas(i))
            indx = [indx; i];
        end
    end
end
close all
clear all

% Aknowledgements
%
% The chl a data and nutrients for south basin are published with the EIDC:
% 
% Maberly, S.C.; Brierley, B.; Carter, H.T.; Clarke, M.A.; De Ville, M.M.; Fletcher, J.M.; James, J.B.; Keenan, P.; 
% Kelly, J.L.; Mackay, E.B.; Parker, J.E.; Patel, M.; Pereira, M.G.; Rhodes, G.; Tanna, B.; Thackeray, S.J.; Vincent, C.J.; Feuchtmayr, H. (2017). 
% Surface temperature, surface oxygen, water clarity, water chemistry and phytoplankton chlorophyll a data from Windermere South Basin, 1945 to 2013. 
% NERC Environmental Information Data Centre. https://doi.org/10.5285/e3c4d368-215d-49b2-8e12-74c99c4c3a9d
% 
% The nutrient concentrations for North Basin are published in the EIDC:
% Maberly, S.C.; Brierley, B.; Carter, H.T.; Clarke, M.A.; De Ville, M.M.; Fletcher, J.M.; James, J.B.; Keenan, P.; Kelly, J.L.; Mackay, E.B.; 
% Parker, J.E.; Patel, M.; Pereira, M.G.; Rhodes, G.; Tanna, B.; Thackeray, S.J.; Vincent, C.J.; Feuchtmayr, H. (2017). Surface temperature, 
% surface oxygen, water clarity, water chemistry and phytoplankton chlorophyll a data from Windermere North Basin, 1945 to 2013. NERC Environmental 
% Information Data Centre. https://doi.org/10.5285/f385b60a-2a6b-432e-aadd-a9690415a0ca

% Sort out directories\paths and switches
C_dir = pwd;addpath([C_dir,'\Sediment-v2.0']);addpath([C_dir,'\Sediment-v2.0\ph-module']);addpath([C_dir,'\IO']);
disp(['Started at:   ', datestr(datetime('now'))]); tic

is_metrics = false; % print metrics in the end
save_initial_conditions = false; % save final concentrations as initial for the next run of sequential years
run_ID = 0;
clim_ID = 0;

%%%% Set start-end date %%%%
m_start = [2009, 1, 1]; %
m_stop = [2009, 12, 31]; %
%%%%%%%%%%% Timeseries of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_of_scenario = 'IO/Scenarios/Wind_2009_Inflow_Chl.csv'%; 
% Initial concentration profiles
init_file = 'Wind_2009_initial_concentrations.txt'
% Mylake results and sediment results save to this file (.mat file)
OUTPUT_file_name = 'IO/Lancaster_JBA_Test_Sim';

%%%%%%%%%%% Load in the parameter values from function %%%%%%%%%%%%%%%%%%%%
[lake_params, sediment_params, tanks_params, algal_params] = load_params();

load([C_dir,'\\IO\Xrec1'],'Xrec1'); % Load parameter sets for acceptable simulations
load([C_dir,'\\IO\Chl_Accept'],'Chl_Accept');

% Run MyLake
[MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
    lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
disp(['Finished at:   ', datestr(datetime('now'))]);

figure; hold on
plot(MyLake_results.basin1.tanks(1).T')
Obs_T = load('SLMSTemp2009.csv');
plot(Obs_T,'ks')
title('Temperature Tank 1')
set(gca,'FontSize',20)


figure; hold on
MyLake_results.basin1.tanks(1).MixStat(12,isnan(MyLake_results.basin1.tanks(1).MixStat(12,:)))=41
plot(42-MyLake_results.basin1.tanks(1).MixStat(12,:),'k-o')
obsmxd = load('MDLIMS_WSB_2009.txt');
plot(41-obsmxd,'r')
title('Mixed depth Tank 1')
legend('Sim','Obs upper','Obs med','Obs lwr')
set(gca,'FontSize',20)

% sum all algae
CHT = zeros(1,365);
CHTD = zeros(1,365);
CHTG = zeros(1,365);
CHTC = zeros(1,365);
% Functional types of algae - loosely Diatoms, Green (other) and
% Cyanobacteria
for ii = 1:6 
CHT = CHT + MyLake_results.basin1.concentrations.tanks(1).Chl{1, ii}(1,:);
if ii < 3
  CHTD = CHTD + MyLake_results.basin1.concentrations.tanks(1).Chl{1, ii}(1,:);  
  end
  if ii > 2 & ii < 5
  CHTG = CHTG + MyLake_results.basin1.concentrations.tanks(1).Chl{1, ii}(1,:);  
  end
  if ii > 4
  CHTC = CHTC + MyLake_results.basin1.concentrations.tanks(1).Chl{1, ii}(1,:);  
  end
end

figure
plot(CHT,'k','linewidth',2)
hold on
plot(CHTD,'r','linewidth',2)
plot(CHTG,'g','linewidth',2)
plot(CHTC,'b','linewidth',2)
obs = load('WSB_CHL_OBS_2009.csv');% Open observed CHL
plot(obs(:,1),obs(:,2),'ks','Markersize',15,'MarkerFaceColor','k')
Obs_CSR = load('CSR.csv');
plot(Obs_CSR(:,1),Obs_CSR(:,5),'ro','Markersize',10,'MarkerFaceColor','r')
plot(Obs_CSR(:,1),Obs_CSR(:,4),'bo','Markersize',10,'MarkerFaceColor','b')
plot(Obs_CSR(:,1),Obs_CSR(:,2)+Obs_CSR(:,3),'go','Markersize',10,'MarkerFaceColor','g')
legend('Sim. TChl','Sim. "Diatoms"','Sim. "Green/other"','Sim. "Cyano"','Obs TChl','Obs D species','Obs CS species','Obs OTHER species')
set(gca,'FontSize',20)


% Plot simulated nutrients
figure; hold on
plot(MyLake_results.basin1.concentrations.tanks(1).Si(:,:)')
title('Si')
figure; hold on
plot(MyLake_results.basin1.concentrations.tanks(1).NO3(:,:)')
title('NO_3')
figure; hold on
plot(MyLake_results.basin1.concentrations.tanks(1).P(:,:)')
title('Phos')






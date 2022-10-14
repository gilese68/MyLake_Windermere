%% Initialise run
close all
clear all

% Sort out directories\paths and switches
C_dir = pwd;addpath([C_dir,'\Sediment-v2.0']);addpath([C_dir,'\Sediment-v2.0\ph-module']);addpath([C_dir,'\IO']);addpath('X:\PROGRAM_FILES\MATLAB\useful_m_files')
disp(['Started at:   ', datestr(datetime('now'))]); tic

is_metrics = false; % print metrics in the end
save_initial_conditions = false; % save final concentrations as initial for the next run of sequential years
run_ID = 0;
clim_ID = 0;

%%%% Set start-end date %%%%
m_start = [2009, 1, 1]; %
m_stop = [2009, 12, 31]; %

%%%% Number of simulations %%%%
N = 1; % Number of acceptable simulations

%%%% Load parameter sets %%%%
load([C_dir,'\\IO\Xrec1'],'Xrec1'); % Load parameter sets for acceptable simulations
load([C_dir,'\\IO\Chl_Accept'],'Chl_Accept');
ParamSet.Tank1 = (Xrec1(1).Tank((Chl_Accept(1,:)),:));
ParamSet.Tank2 = (Xrec1(2).Tank((Chl_Accept(1,:)),:));

%%% Adjust tank 2 for existing FPV array
ParamSet.Tank2(:,1) = ParamSet.Tank2(:,1)*0.81; % Wind shelter coefficient
ParamSet.Tank2(:,8) = ParamSet.Tank2(:,8)*0.812; % Solar shelter coefficient
ParamSet.Tank2(:,9) = ParamSet.Tank2(:,9)*1.016; % Air shelter coefficient

%%%%%%%%%%% Timeseries of inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_of_scenario = 'IO/Scenarios/Wind_2009_Inflow_Chl.csv'%; % P_gradual_increase_no_chl_2015_cutoff_to_0_2200.txt'
% Initial concentration profiles
init_file = 'Wind_2009_initial_concentrations.txt'
% Mylake results and sediment results save to this file (.mat file)
OUTPUT_file_name = 'IO/Lancaster_JBA_Test_Sim'; 

%%%%%%%%%%% Load in the parameter values from function %%%%%%%%%%%%%%%%%%%%
[lake_params, sediment_params, tanks_params, algal_params] = load_params();

%% Experiment 1
% Test the transferability of the model between QEII and Windermere SB
% No floating solar coverage

Sn = 0; % Number of reductions in this scenario

% Manually overwrite paramters
lake_params{5} = {0,0}; % Wind shelter
lake_params{38} = {0.25,0.25}; % Non_PAR
lake_params{39} = {0.25,0.25}; % PAR

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\Windermere_Test\Wind_Ex_1\Wind_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 1; % Which scenario number
Sn = 0;
Sn2 = 1;
Tk = 1;
Type = 1;
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\Windermere_Test\Wind_Ex_1\Experiment_1_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plotting routine
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\Windermere_Test\Wind_Ex_1\Experiment_1_Output'])
ScenarioNum = 1; % Which scenario number

% Plot observed and simulated Chl-a

for ii = 1:1
    plot(Chlorophyll_1.tanks.Reduction.RunNum(ii).Chla_Top_Mean)
    hold on
end
hold on
plot(WNBEIDCS14.Min,'r')
hold on
plot(WNBEIDCS14.Max,'r')
hold off

for ii = 1:75
    plot(Temperature_1.tanks.Reduction.RunNum(ii).Temp(2,:))
    hold on
end
hold on
plot(WNBEIDCS17.Min,'r')
hold on
plot(WNBEIDCS17.Max,'r')

% Plot mixed depth
Mixed_Depth_1.tanks.Sim.Mxd(isnan(Mixed_Depth_1.tanks.Sim.Mxd))=42
plot(Mixed_Depth_1.tanks.Sim.Mxd)
hold on
plot(ObsMxd)'
set(gca, 'YDir','reverse')

% Plot all temps
for ii = 1:42
    plot(Temperature_1.tanks.Reduction.RunNum.Temp(ii,:))
    hold on
end

%% Scenario 1
% Tank 1 sited – 10% coverage increments (solar reduced by Langthwaite, wind by AA’s data)

Sn = 10; % Number of reductions in this scenario
Wr = [1:-((1/Sn)*0.95):0]; % Wind speed reductions (increasing coverage)
Sr = [1:-((1/Sn)*0.94):0]; % Solar radiation reductions (increasing coverage)
Ar = [1:(0.08/Sn):1.08]; % Air temperature warming (increasing coverage)

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*Wr(hh)),ParamSet.Tank2(ii,1)}; % Wind "shelter/enhancement! factor 0.5
        %lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*Sr(hh)),ParamSet.Tank2(ii,8)}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*Ar(hh)),ParamSet.Tank2(ii,9)}; % Air shelter
%         lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
%         lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_1\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 1; % Which scenario number
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,N,Tk,ScenarioNum); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_1\Scenario_1_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plotting routine
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_1\Scenario_1_Output'])
ScenarioNum = 1; % Which scenario number

%% Scenario 2
% Tank 2 sited 10% coverage increments (solar reduced by Langthwaite, wind by AA’s data)

Sn = 10; % Number of reductions in this scenario
Wr = [1:-((1/Sn)*0.95):0]; % Wind speed reductions (increasing coverage)
Sr = [1:-((1/Sn)*0.94):0]; % Solar radiation reductions (increasing coverage)
Ar = [1:(0.08/Sn):1.08]; % Air temperature warming (increasing coverage)

for hh = 11:Sn+1% Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {ParamSet.Tank1(ii,1),(ParamSet.Tank2(ii,1)*Wr(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {ParamSet.Tank1(ii,8),(ParamSet.Tank2(ii,8)*Sr(hh))}; % Solar shelter
        lake_params{37} = {ParamSet.Tank1(ii,9),(ParamSet.Tank2(ii,9)*Ar(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_2\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 2; % Which scenario number
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,N,Tk,ScenarioNum); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_2\Scenario_2_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

%% Scenario 3
% Tank 1 & 2 sited – split the coverage between tanks proportionally to
% their area, 10% coverage increments (solar reduced by Langthwaite, wind by AA’s data) 

Sn = 10; % Number of reductions in this scenario
Wr = [1:-((1/Sn)*0.95):0]; % Wind speed reductions (increasing coverage)
Sr = [1:-((1/Sn)*0.94):0]; % Solar radiation reductions (increasing coverage)
Ar = [1:(0.08/Sn):1.08]; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 0.7; % Proportion occupied by Tank 1
T2a = 0.3; % Proportion occupied by Tank 2

% Modify reductions by tank proportions
WrT1 = 1+((Wr-1)*T1a); % Wind speed reductions (Tank 1)
WrT2 = 1+((Wr-1)*T2a); % Wind speed reductions (Tank 2)
SrT1 = 1+((Sr-1)*T1a); % Solar radiation reductions (Tank 1)
SrT2 = 1+((Sr-1)*T2a); % Solar radiation reductions (Tank 2)
ArT1 = 1+((Ar-1)*T1a); % Air temperature warming (Tank 1)
ArT2 = 1+((Ar-1)*T2a); % Air temperature warming (Tank 2)

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1),(ParamSet.Tank2(ii,1)*WrT2)}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1),(ParamSet.Tank2(ii,8)*SrT2)}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1),(ParamSet.Tank2(ii,9)*ArT2)}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_3\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 3; % Which scenario number
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,N,Tk,ScenarioNum); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_3\Scenario_3_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

%% Scenario 4 (Scenario 5 - inflow volumes)
% Tank 1 and Tank 2 sited (each tank occupies 50% of the reservoir), 
% inflow scaling factor modified for each group of runs, 25% coverage
% increments evenly split 50/50 between each tank

tanks_params{4} = 'IO/test_tanks_bathymetry_2tanks_50.csv';
tanks_params{11} = {'IO/QE2/Tanks/QEII_2018_initial_concentrations_1_50.txt','IO/QE2/Tanks/QEII_2018_initial_concentrations_2_50.txt'};

Sn = 4; % Number of reductions in this scenario
Sn2 = 4; % Number of variable scenarios e.g. inflow factors
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)
InflowFactor = [0.25, 0.5, 1, 1.25]; % Inflow factor

T1b = [0, 0.25, 0.5, 0.75, 1]; % Proportion of Tank 1 covered by array
T2b = [0, 0.25, 0.5, 0.75, 1]; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for gg = 1:length(InflowFactor) % Run for the number of inflow factors

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{14} = {InflowFactor(1,gg),InflowFactor(1,gg)}; % Inflow scaling factor
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d_%d';
        str = sprintf(formatSpec,gg,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_4\QE2_',str],'MyLake_results');
        
    end
end
end

% Post processing
ScenarioNum = 4; % Which scenario number
Type = 2; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
Sn = 4; % Number of reductions in this scenario
Sn2 = [0.25, 0.5, 1, 1.25]; % Number of variable scenarios e.g. inflow factors
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_4\Scenario_4_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_4\Scenario_4_Output'])
ScenarioNum = 4; % Which scenario number
Sn = 4; % Number of reductions in this scenario
Sn2 = [0.25, 0.5, 1, 1.25]; % Number of variable scenarios e.g. inflow factors
[PLP_rec, Mxd_tot] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk)

%%%%%%%%%%% Load in the parameter values from function %%%%%%%%%%%%%%%%%%%%
[lake_params, sediment_params, tanks_params, algal_params] = load_params();

%% Scenario 5 (Scenario 2)
% Array sited on tank 2, 10% coverage intervals which spill over to cover
% tank 1 untill the reservoir is fully covered

Sn = 10; % Number of reductions in this scenario
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

T1b = [0, 0, 0, 0, 0.142857143,	0.285714286, 0.428571429, 0.571428571, 0.714285714, 0.857142857, 1]; % 10% of total/tank area (0.1*T3a)/T1a
T2b = [0, 1/3, 2/3, 1, 1, 1, 1, 1, 1, 1, 1]; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 5; % Which scenario number
Type = 1; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Tk = 2; % Number of tanks
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\Scenario_5_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\Scenario_5_Output'])
ScenarioNum = 5; % Which scenario number
[PLP_rec, Mxd_tot, P1_rec] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk);
PLP_rec_5 = PLP_rec;
P1_rec_5 = P1_rec;
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\PLP_rec_5'],'PLP_rec_5');
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\P1_rec_5'],'P1_rec_5');

%% Scenario 6 (Scenario 3)
% Array sited centrally on reservoir. Once the uncovered area of tank 1 is 
% equal to tank 2 the array is distributed equally between each untill the
% reservoir is fully covered

Sn = 10; % Number of reductions in this scenario
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

% Units of tank covered
T1Cov = [0, 10:10:40, 45:5:70]; % Tank 1 (is 70% of reservoir, so maximum units is 70)
T2Cov = [0, 0, 0, 0, 0, 5:5:30]; % Tank 2 (is 30% of reservoir, so maximum units is 30)

T1b = (T1Cov * ((Sn/(T1a/T3a))))/1000; % Proportion of Tank 1 covered by array
T2b = (T2Cov * ((Sn/(T2a/T3a))))/1000; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for gg = 1:Sn2 % Run for the number of variable scenarios
for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\QE2_',str],'MyLake_results');
        
    end
end
end

% Post processing
ScenarioNum = 6; % Which scenario number
Type = 1; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Tk = 2; % Number of tanks
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\Scenario_6_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\Scenario_6_Output'])
ScenarioNum = 6; % Which scenario number
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
[PLP_rec, Mxd_tot, P1_rec] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk);
PLP_rec_6 = PLP_rec;
P1_rec_6 = P1_rec;
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\PLP_rec_6'],'PLP_rec_6');
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\P1_rec_6'],'P1_rec_6');

%% Scenario 7 (Scenario 4)
% Reservoir volume reduced to 75% - 4m shallower
% Array sited centrally on reservoir. Once the uncovered area of tank 1 is 
% equal to tank 2 the array is distributed equally between each untill the
% reservoir is fully covered

tanks_params{4} = 'IO/test_tanks_bathymetry_2tanks_75vol.csv';
tanks_params{11} = {'IO/QE2/Tanks/QEII_2018_initial_concentrations_1_75vol.txt','IO/QE2/Tanks/QEII_2018_initial_concentrations_2_75vol.txt'};

Sn = 10; % Number of reductions in this scenario
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

% Units of tank covered
T1Cov = [0, 10:10:40, 45:5:70]; % Tank 1 (is 70% of reservoir, so maximum units is 70)
T2Cov = [0, 0, 0, 0, 0, 5:5:30]; % Tank 2 (is 30% of reservoir, so maximum units is 30)

T1b = (T1Cov * ((Sn/(T1a/T3a))))/1000; % Proportion of Tank 1 covered by array
T2b = (T2Cov * ((Sn/(T2a/T3a))))/1000; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_7\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 7; % Which scenario number
Tk = 2; % Number of tanks in scenario
Type = 1; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_7\Scenario_7_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_7\Scenario_7_Output'])
ScenarioNum = 7; % Which scenario number
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
[PLP_rec, Mxd_tot] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk)

%%%%%%%%%%% Load in the parameter values from function %%%%%%%%%%%%%%%%%%%%
[lake_params, sediment_params, tanks_params, algal_params] = load_params();

%% Scenario 8 (Scenario 1)
% Array sited on tank 1, 10% coverage intervals which spill over to cover
% tank 2 untill the reservoir is fully covered

Sn = 10; % Number of reductions in this scenario
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

T1b = [0, 0.142857143, 0.285714286, 0.428571429, 0.571428571, 0.714285714, 0.857142857, 1, 1, 1, 1,]; % 10% of total/tank area (0.1*T3a)/T1a
T2b = [0, 0, 0, 0, 0, 0, 0, 0, 1/3, 2/3, 1]; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\QE2_',str],'MyLake_results');
        
    end
end

% Post processing
ScenarioNum = 8; % Which scenario number
Type = 1;
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Tk = 2;
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\Scenario_8_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\Scenario_8_Output'])
ScenarioNum = 8; % Which scenario number
Type = 1;
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
[PLP_rec, Mxd_tot, P1_rec] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk);
PLP_rec_8 = PLP_rec;
P1_rec_8 = P1_rec;
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\PLP_rec_8'],'PLP_rec_8');
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\P1_rec_8'],'P1_rec_8');

%% Scenario 9 - new loss rates
% Array sited centrally on reservoir. Once the uncovered area of tank 1 is 
% equal to tank 2 the array is distributed equally between each untill the
% reservoir is fully covered

Sn = 10; % Number of reductions in this scenario
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

% Units of tank covered
T1Cov = [0, 10:10:40, 45:5:70]; % Tank 1 (is 70% of reservoir, so maximum units is 70)
T2Cov = [0, 0, 0, 0, 0, 5:5:30]; % Tank 2 (is 30% of reservoir, so maximum units is 30)

T1b = (T1Cov * ((Sn/(T1a/T3a))))/1000; % Proportion of Tank 1 covered by array
T2b = (T2Cov * ((Sn/(T2a/T3a))))/1000; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for gg = 1:Sn2 % Run for the number of variable scenarios
for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_9\QE2_',str],'MyLake_results');
        
    end
end
end

% Post processing
ScenarioNum = 9; % Which scenario number
Type = 1; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_9\Scenario_9_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_9\Scenario_9_Output'])
ScenarioNum = 9; % Which scenario number
Type = 1;
Sn2 = 1; % Number of variable scenarios e.g. inflow factors
Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk)

%% Scenario 10 (All coverages modified inflow)
% Tank 1 and Tank 2 sited (each tank occupies 50% of the reservoir), 
% inflow scaling factor modified for each group of runs, 25% coverage
% increments evenly split 50/50 between each tank

Sn = 10; % Number of reductions in this scenario
Sn2 = 4; % Number of variable scenarios e.g. inflow factors
Wr = 0.95; % Wind speed reductions (increasing coverage)
Sr = 0.94; % Solar radiation reductions (increasing coverage)
Ar = 1.08; % Air temperature warming (increasing coverage)
InflowFactor = [0.25, 0.5, 1, 1.25]; % Inflow factor

% Tank areas/proportions
T1a = 882000; % Tank 1 area
T2a = 378000; % Tank 2 area
T3a = T1a+T2a; % Total area

% Units of tank covered
T1Cov = [0, 10:10:40, 45:5:70]; % Tank 1 (is 70% of reservoir, so maximum units is 70)
T2Cov = [0, 0, 0, 0, 0, 5:5:30]; % Tank 2 (is 30% of reservoir, so maximum units is 30)

T1b = (T1Cov * ((Sn/(T1a/T3a))))/1000; % Proportion of Tank 1 covered by array
T2b = (T2Cov * ((Sn/(T2a/T3a))))/1000; % Proportion of Tank 2 covered by array

% Modify reductions by tank proportions
WrT1b = 1-(Wr*T1b); % Wind speed reductions (Tank 1)
WrT2b = 1-(Wr*T2b); % Wind speed reductions (Tank 2)
SrT1b = 1-(Sr*T1b); % Solar radiation reductions (Tank 1)
SrT2b = 1-(Sr*T2b); % Solar radiation reductions (Tank 2)
ArT1b = 1+((Ar-1)*T1b); % Air temperature warming (Tank 1)
ArT2b = 1+((Ar-1)*T2b); % Air temperature warming (Tank 2)

for gg = 1:length(InflowFactor) % Run for the number of inflow factors

for hh = 1:Sn+1 % Run all variations of scenario
    
    for ii = 1:N % Run all the acceptable simulations
        % Manually overwite parameters to the acceptable simulations
        lake_params{5} = {(ParamSet.Tank1(ii,1)*WrT1b(hh)),(ParamSet.Tank2(ii,1)*WrT2b(hh))}; % Wind "shelter/enhancement! factor 0.5
        lake_params{14} = {InflowFactor(1,gg),InflowFactor(1,gg)}; % Inflow scaling factor
        lake_params{15} = {ParamSet.Tank1(ii,4),ParamSet.Tank2(ii,4)}; % Inflow temperature
        lake_params{36} = {(ParamSet.Tank1(ii,8)*SrT1b(hh)),(ParamSet.Tank2(ii,8)*SrT2b(hh))}; % Solar shelter
        lake_params{37} = {(ParamSet.Tank1(ii,9)*ArT1b(hh)),(ParamSet.Tank2(ii,9)*ArT2b(hh))}; % Air shelter
        lake_params{38} = {ParamSet.Tank1(ii,10),ParamSet.Tank2(ii,10)}; % Non_PAR
        lake_params{39} = {ParamSet.Tank1(ii,11),ParamSet.Tank2(ii,11)}; % PAR
        
        algal_params{1} = {ParamSet.Tank1(ii,15),ParamSet.Tank2(ii,15),ParamSet.Tank1(ii,16),ParamSet.Tank2(ii,16),ParamSet.Tank1(ii,17),ParamSet.Tank2(ii,17)}; % PAR Sat
        algal_params{5} = {ParamSet.Tank1(ii,12),ParamSet.Tank2(ii,12),ParamSet.Tank1(ii,13),ParamSet.Tank2(ii,13),ParamSet.Tank1(ii,14),ParamSet.Tank2(ii,14)}; % Growth
        
        % Run MyLake
        [MyLake_results, Sediment_results]  = fn_MyL_application(m_start, m_stop, sediment_params,...
            lake_params, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, save_initial_conditions, init_file); % runs the model and outputs obs and sim
        disp(['Finished at:   ', datestr(datetime('now'))]);
        
        % Save output
        formatSpec = '%d_%d_%d';
        str = sprintf(formatSpec,gg,hh,ii);
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_10\QE2_',str],'MyLake_results');
        
    end
end
end

% Post processing
ScenarioNum = 10; % Which scenario number
Type = 2; % Type 1 is one variable, type 2 is multiple variables, e.g. varying inflow
Sn = 10; % Number of reductions in this scenario
Sn2 = [0.25, 0.5, 1, 1.25]; % Number of variable scenarios e.g. inflow factors
[Chlorophyll_1, Temperature_1, Mixed_Depth_1] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type); % Function for post-processing
save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_10\Scenario_10_Output'],...
    'Chlorophyll_1', 'Temperature_1', 'Mixed_Depth_1'); % Save to file

%% Plot the scenario and save to file
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_4\Scenario_4_Output'])
ScenarioNum = 4; % Which scenario number
Sn = 4; % Number of reductions in this scenario
Sn2 = [0.25, 0.5, 1, 1.25]; % Number of variable scenarios e.g. inflow factors
[PLP_rec, Mxd_tot] = Results_Plotting(Chlorophyll_1, Temperature_1, Mixed_Depth_1, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk)

%%%%%%%%%%% Load in the parameter values from function %%%%%%%%%%%%%%%%%%%%
[lake_params, sediment_params, tanks_params, algal_params] = load_params();

%% Plot multiple  maximum temperature scenarios on a box plot
% Load data for box plots
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\PLP_rec_5'],'PLP_rec_5');
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\PLP_rec_6'],'PLP_rec_6');
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\PLP_rec_8'],'PLP_rec_8');

% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [18.05 18.05 23.95 23.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [18.05 18.05 23.95 23.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).WaterTemp_Max(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    hold on
end
PLP = PLP_rec_5.tanks(1).WaterTemp_Max(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    hold on
end
PLP = PLP_rec_6.tanks(1).WaterTemp_Max(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Temperature distribution on day of maximum temperature')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Water Temperature (^{o}C)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([18 24])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\MaxTemp')

%% Plot multiple minimum temperature scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [1.55 1.55 4.95 4.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [1.55 1.55 4.95 4.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).WaterTemp_Min(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    hold on
end
PLP = PLP_rec_5.tanks(1).WaterTemp_Min(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    hold on
end
PLP = PLP_rec_6.tanks(1).WaterTemp_Min(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Temperature distribution on day of minimum temperature')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Water Temperature (^{o}C)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([1.5 5])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\MinTemp')

%% Plot multiple max stratification duration scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 34.5 34.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 34.5 34.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Stratification duration')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Stratification duration (days)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 35])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Strat_Max_Duration')

%% Plot multiple cumulative stratification duration scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Cumul(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Cumul(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Cumul(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Cumulative stratification duration')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Cumulative stratification (days)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 80])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Strat_Cumul_Duration')

%% Plot multiple stratification onset scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Start(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Start(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Start(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
hold on
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Stratification onset')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Stratification onset (day of year)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 250])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Strat_Start')

%% Plot multiple stratification overturn scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_End(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_End(5).PLP; % Scenario 2
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_End(5).PLP; % Scenario 3
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
hold on
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Stratification overturn')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Stratification overturn (day of year)')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 250])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Strat_End')

%% Tiled layout - physical parameters
t = tiledlayout(3,2);

% Tile 1 - max temperature
ax1 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [18.05 18.05 23.95 23.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [18.05 18.05 23.95 23.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).WaterTemp_Max(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    hold on
end
PLP = PLP_rec_5.tanks(1).WaterTemp_Max(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    hold on
end
PLP = PLP_rec_6.tanks(1).WaterTemp_Max(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'a) Temperature distribution on day of maximum temperature','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%xlabel('Floating solar array coverage (%)')
ylabel('Water Temperature (^{o}C)')
% legend([b1 b2 b3],'Scenario-Fast','Scenario-Short','Scenario-Central')
xlim([0 340])
ylim([18 24])
set(gca,'FontSize',12)

% Tile 2 - minimum temperature
ax2 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [1.55 1.55 4.95 4.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [1.55 1.55 4.95 4.95];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).WaterTemp_Min(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    hold on
end
PLP = PLP_rec_5.tanks(1).WaterTemp_Min(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    hold on
end
PLP = PLP_rec_6.tanks(1).WaterTemp_Min(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax2,'b) Temperature distribution on day of minimum temperature','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%xlabel('Floating solar array coverage (%)')
ylabel('Water Temperature (^{o}C)')
% legend([b1 b2 b3],'Scenario-Fast','Scenario-Slow','Scenario-Central')
xlim([0 340])
ylim([1.5 5])
set(gca,'FontSize',12)

% Tile 3 - stratification duration
ax3 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 34.5 34.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 34.5 34.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Max_Dur(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax3.FontSize = 12
ax3.TitleHorizontalAlignment = 'left';
subtitle(ax3,'c) Stratification duration','FontSize',10)
ax3.TitleFontSizeMultiplier = 0.83
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%xlabel('Floating solar array coverage (%)')
ylabel('Stratification duration (days)')
% legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 35])
set(gca,'FontSize',12)

% Tile 4 - stratification duration (cumulative)
ax4 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Cumul(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Cumul(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Cumul(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax4.FontSize = 12
ax4.TitleHorizontalAlignment = 'left';
subtitle(ax4,'d) Cumulative stratification duration','FontSize',10)
ax4.TitleFontSizeMultiplier = 0.83;
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%xlabel('Floating solar array coverage (%)')
ylabel('Cumulative stratification (days)')
% legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 80])
set(gca,'FontSize',12)

% Tile 5 - onset
ax5 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_Start(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_Start(5).PLP; % Scenario 2
for jj = 1:Sn+1
   if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_Start(5).PLP; % Scenario 3
for jj = 1:Sn+1
     if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
hold on
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax5.FontSize = 12
ax5.TitleHorizontalAlignment = 'left';
subtitle(ax5,'e) Stratification onset','FontSize',10)
ax5.TitleFontSizeMultiplier = 0.83;
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Stratification onset (day of year)')
% legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 250])
set(gca,'FontSize',12)

% Tile 6 - overturn
ax6 = nexttile
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [2 2 248 248];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Strat_End(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Strat_End(5).PLP; % Scenario 2
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Strat_End(5).PLP; % Scenario 3
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
hold on
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
ax6.FontSize = 12
ax6.TitleHorizontalAlignment = 'left';
subtitle(ax6,'f) Stratification overturn','FontSize',10)
ax6.TitleFontSizeMultiplier = 0.83;
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Stratification overturn (day of year)')
% legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 250])
set(gca,'FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';

clear legend
%legend = legend(labels,'NumColumns',3,'Orientation','horizontal');
legend = legend([b1 b2 b3],'Scenario-Fast','Scenario-Slow','Scenario-Central','NumColumns',3,'Orientation','horizontal')
title(legend,{'Floating Solar Deployment Scenario'},'FontSize',10)
% title(legend,{'Floating Solar Coverage Increment (%)'},'FontSize',10)
legend.Layout.Tile = 'south';

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Physical_Params')

%% Plot multiple polymixis scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.005 0.005 0.995 0.995];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.005 0.005 0.995 0.995];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Polymixis(5).PLP; % Scenario 1 PLP_rec_5.tanks(1).Strat_Max_Dur(5).PLP
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj),1,'*','Color',[0.1059, 0.6196, 0.4667],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    end
    hold on
end
hold on
PLP = PLP_rec_5.tanks(1).Polymixis(5).PLP; % Scenario 2
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+10,1,'*','Color',[0.8510, 0.3725, 0.0078],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    end
    hold on
end
hold on
PLP = PLP_rec_6.tanks(1).Polymixis(5).PLP; % Scenario 3
for jj = 1:Sn+1
    if PLP(jj,3) == 0 % Show marker if no stratification
        plot(G(jj)+20,1,'*','Color',[0.4588, 0.4392, 0.7020],'MarkerSize',10)
    else
        % Median ,pecentile 1,percentile 3,max,min
        SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    end
    hold on
end
hold on
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Polymixis')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Polymixis metric')
legend([b1 b2 b3],'Scenario 1','Scenario 2','Scenario 3')
xlim([0 340])
ylim([0 1])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Polymixis')

%% Plot multiple total chlorophyll-a scenarios on a box plot
% Plotting routine
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 39.5 39.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 39.5 39.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).Chl(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),[0.1059, 0.6196, 0.4667],1);
    hold on
end
PLP = PLP_rec_5.tanks(1).Chl(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,[0.8510, 0.3725, 0.0078],1);
    hold on
end
PLP = PLP_rec_6.tanks(1).Chl(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr_triplet(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,[0.4588, 0.4392, 0.7020],1);
    hold on
end
b1 = bar(350,24,'FaceColor',[0.1059, 0.6196, 0.4667]);
hold on
b2 = bar(351,24,'FaceColor',[0.8510, 0.3725, 0.0078]);
hold on
b3 = bar(352,24,'FaceColor',[0.4588, 0.4392, 0.7020]);
hold on
title('Chlorophyll-a distribution on day of maximum Chlorophyll-a')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Chlorophyll-a distribution (µg L^{-1})')
legend([b1 b2 b3],'Scenario-Fast','Scenario-Slow','Scenario-Central')
xlim([0 340])
ylim([0 40])
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Chl_Total')

%% Plot number of stratified scenarios bar chart
Scr1 = [75 75 75 75 62 0 0 2 7 4 7];
Scr2 = [75 75 75 75 75 75 75 21 0 0 7];
Scr3 = [75 75 75 75 62 2 0 0 0 7 7];
Y = [Scr1; Scr2; Scr3]';
X = [20:30:320];
G = 10:30:310;
x = [5 35 35 5];
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [0.5 0.5 79.5 79.5];
patch(x,y,[0.95 0.95 0.95],'EdgeColor','none')
end
hold on
b = bar(X,Y);
b(1).FaceColor = [0.1059, 0.6196, 0.4667];
b(2).FaceColor = [0.8510, 0.3725, 0.0078];
b(3).FaceColor = [0.4588, 0.4392, 0.7020];
legend([b(1) b(2) b(3)],'Scenario-Fast','Scenario-Slow','Scenario-Central')
xlabel('Floating solar array coverage (%)')
ylabel('Number of scenarios with stratification')
xticks([20:30:320]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
set(gca,'FontSize',12)
saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Num_Strat')

%% Plot differences and comparison plots
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_8\P1_rec_8'],'P1_rec_8'); % Scenario 1
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_5\P1_rec_5'],'P1_rec_5'); % Scenario 2
load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\P1_rec_6'],'P1_rec_6'); % Scenario 3

%% Plot annual surface temperature (using the median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
ylabel('Water Temperature (^{o}C)')
datetick('x','mmm')
xlim ([0 365])
lgd = legend({'cos(x)','cos(2x)','cos(3x)','cos(4x)'},'Location','bestoutside')
title(lgd,'{Floating solar}{coverage increment}')
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled layout annual surface temperature (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
Sn = 10;
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
% a1 = plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51])
% legend([a1], {'Air Temperature'},'FontSize',8,'Location','northwest')
% legend('boxoff')
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario 1','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
%plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51])
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle('Scenario 2','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
%plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51])
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle('Scenario 3','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'x')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Median Water Temperature (^{o}C)','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend = legend(labels,'NumColumns',6,'Orientation','horizontal');
title(legend,{'Floating Solar Coverage Increment (%)'},'FontSize',10)
legend.Layout.Tile = 'south';

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Surf_Temp')

%% Plot annual absolute temperature difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
ylabel('Absolute Water Temperature Difference (^{o}C)')
datetick('x','mmm')
xlim ([0 365])
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled layout annual absolute temperature difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
newcolors = {'#f71735','#090446','#f7cb15','#0b6e4f','#8b6220','#008dd5','#24bc33','#ed33b9','#7d6b91','#80a300','#262a10'};
colororder(newcolors)
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario 1','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax2,'Scenario 2','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax3,'Scenario 3','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'x')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Absolute Water Temperature difference (^{o}C)','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend(labels,'Location','southoutside','NumColumns',6,'Orientation','horizontal');
title(legend,'FPV Coverage Increment (%)','FontSize',10)
legend.Layout.Tile = 'south';

set(gcf, 'Position',  [100, 100, 750, 600]);

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Surf_Temp_Absolute')

%% Plot annual relative temperature difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = ((P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = ((P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
ylabel('Relative Difference (%)')
datetick('x','mmm')
xlim ([0 365])
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = ((P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled annual relative temperature difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
Sn = 10;
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_8.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_8.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end

yyaxis right
a1 = plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51],'LineWidth',2)
legend([a1], {'Air Temperature'},'FontSize',8,'Location','best','Position',[0.778754700156119 0.743781703753327 0.13361168862758 0.0188679240909513])
legend('boxoff')
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario 1','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_5.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_5.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
yyaxis right
plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51],'LineWidth',2)
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle('Scenario 2','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83
ylabel('Observed air temperature (^{o}C)','FontSize',12)

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_6.tanks(1).WaterTemp(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3))./P1_rec_6.tanks(1).WaterTemp(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
yyaxis right
plot(SimDate,AirTemperature,'Color',[0.42 0.47 0.51],'LineWidth',2)
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle('Scenario 3','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'x')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Relative Temperature Difference (%)','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend = legend(labels,'NumColumns',6,'Orientation','horizontal');
title(legend,{'Floating Solar Coverage Increment (%)'},'FontSize',10)
legend.Layout.Tile = 'south';

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Surf_Temp_relative_+Air')

%% Plot annual total chlorophyll-a (using the median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
ylabel('Chlorophyll-a distribution (µg L^{-1})')
datetick('x','mmm')
xlim ([0 365])
% lgd = legend({'cos(x)','cos(2x)','cos(3x)','cos(4x)'},'Location','bestoutside')
% title(lgd,'{Floating solar}{coverage increment}')
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled layout annual total chlorophyll-a (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario 1','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax2,'Scenario 2','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax3,'Scenario 3','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'x')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Chlorophyll-a distribution (µg L^{-1})','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend(labels,'Location','southoutside','NumColumns',6,'Orientation','horizontal');
title(legend,'FPV Coverage Increment (%)','FontSize',10)
legend.Layout.Tile = 'south';

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Chl')

%% Plot annual absolute total chlorophyll-a difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
ylabel('Absolute Chlorophyll-a distribution (µg L^{-1})')
datetick('x','mmm')
xlim ([0 365])
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled layout absolute annual total chlorophyll-a difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
newcolors = {'#f71735','#090446','#f7cb15','#0b6e4f','#8b6220','#008dd5','#24bc33','#ed33b9','#7d6b91','#80a300','#262a10'};
colororder(newcolors)
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario-Fast','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax2,'Scenario-Slow','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3);
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax3,'Scenario-Central','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'xy')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
ax1.YLim = ([-35 5])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Absolute total Chlorophyll-a difference (µg L^{-1})','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend(labels,'Location','southoutside','NumColumns',6,'Orientation','horizontal');
title(legend,'FPV Coverage Increment (%)','FontSize',10)
legend.Layout.Tile = 'south';

set(gcf, 'Position',  [100, 100, 750, 600]);

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Chl_Absolute')

%% Plot annual relative temperature difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
subplot(3,1,1)
for jj = 1:Sn+1
    Curve = ((P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
subplot(3,1,2)
for jj = 1:Sn+1
    Curve = ((P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
ylabel('Relative Chlorophyll-a distribution difference (%)')
datetick('x','mmm')
xlim ([0 365])
hold on
subplot(3,1,3)
for jj = 1:Sn+1
    Curve = ((P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
xlim ([0 365])
xlabel('Date')

%% Tiled layout annual total chlorophyll-a relative difference (median only)
SimDate = 1:length(datenum(m_start):datenum(m_stop));
clear legend
t = tiledlayout(3,1);

% Tile 1
ax1 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_8.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_8.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax1,'Scenario 1','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 2
ax2 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_5.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_5.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax2,'Scenario 2','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

% Tile 3
ax3 = nexttile
for jj = 1:Sn+1
    Curve = ((P1_rec_6.tanks(1).Chl(1).P1(jj).P1(:,3) - P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3))./P1_rec_6.tanks(1).Chl(1).P1(1).P1(:,3))*100;
    plot(SimDate, Curve)
    hold on
end
datetick('x','mmm')
ax=gca;
ax.FontSize = 12
ax.TitleHorizontalAlignment = 'left';
subtitle(ax3,'Scenario 3','FontSize',10)
ax.TitleFontSizeMultiplier = 0.83

xlabel('Date','FontSize',12)
linkaxes([ax1 ax2 ax3],'x')
xticklabels([ax1, ax2],{})
datetick('x','mmm')
ax1.XLim = ([0 365])
t.TileSpacing = 'compact';
t.Padding = 'compact';
ylabel(t,'Relative Chlorophyll-a distribution difference (%)','FontSize',12)

labels = {'0','10','20','30','40','50','60','70','80','90','100'};
legend(labels,'Location','southoutside','NumColumns',6,'Orientation','horizontal');
title(legend,'FPV Coverage Increment (%)','FontSize',10)
legend.Layout.Tile = 'south';

saveas(gcf,'\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions\\QE2_Scenario_1_2_3\\Annual_Chl_Relative')

%% Rename files
for hh = 1:Sn+1 % Run all variations of scenario
    for ii = 1:N % Run all the acceptable simulations
        formatSpec = '%d_%d';
        str = sprintf(formatSpec,hh,ii);
        load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\QE2_1_',str],'MyLake_results');
        save(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_6\QE2_',str],'MyLake_results');
    end
end
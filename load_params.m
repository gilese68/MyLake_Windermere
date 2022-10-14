%% load_params: load parameters for MyLake and Sediment
function [lake_params, sediment_params, tanks_params, algal_params] = load_params()

tanks_params = {
    % TanksPar
    2, 'Nt', % 1 Number of tanks
    'IO/test_advection_matrix.xlsx', 'tank_advection_matrix', % Path to file containing the advection matrix
    'IO/test_diffusion_matrix.xlsx', 'tank_diffusion_matrix', % Path to file containing the diffusion matrix
    'IO/test_tanks_bathymetry.csv', 'tank_areas', % Path to file containing the tanks bathymetry
    1, 'river_inflow_switch', % River inflow switch (0:no river inflow, 1: v2 river inflow module for tank with river inflow, 2: exponential decline for river inflow 3: exponential flipped (flow in from bottom to simulate reservoir bottom-fed inputs))
    1, 'river_inflow_tank_id', % Tank ID with river inflow, required only if river_inflow_switch = 1
    [1,0], 'river_in_coef', % River inflow coefficient per tank, put to zero if using river_inflow_switch e.g. [1,0,0,0,0,0]
    [1,0], 'river_out_coef',% River outflow coefficient per tank e.g. [0,0,0,0,0,1]
    0.75, 'decline_parameter_advection', % Exponential decline parameter for advection
    0.75, 'decline_parameter_diffusion', % Exponential decline parameter for diffusion
    {'IO/Input_Wind/Tanks/Wind_2009_initial_concentrations_1.txt','IO/Input_Wind/Tanks/Wind_2009_initial_concentrations_2.txt'}, 'init_filename_list'
    };

%% PhysPar
C = readcell('IO/lake_params.xlsx');

lake_params = cell(61,2);
lake_params(:,1) = {cell(1,(size(C,2)-1))};

for i = 1:(size(C,2)-1)
    for ii = 1:size(C,1)
        lake_params{ii, 1}{1, i} = C{ii, i+1};
    end
end

% Import parameter names
for i = 1:size(C,1)
    lake_params{i,2} = C{i,1};
end

lake_params = {
    % PhysPar
    {1,1} 'dz',                             % 1  depth resolution (m)    
    {0.0442401,0.0442401}, 'Kz_K1',         % 2     open water diffusion parameter (-)
    {0.000898,0.000898}, 'Kz_K1_ice',       % 3     under ice diffusion parameter (-)
    {7.0e-05,7.0e-05}, 'Kz_N0',             % 4     min. stability frequency (s-2)
    {1,1},  'C_shelter',  %0.5,       % 5     wind shelter parameter (-)
    {54.3,54.3}, 'lat',                     % 6     latitude (decimal degrees)
    {2.9,2.9}, 'lon',                       % 7     longitude (decimal degrees)
    {0.3, 0.3},'alb_melt_ice',              % 8     albedo of melting ice (-)
    {0.77,0.77}, 'alb_melt_snow',           % 9     albedo of melting snow (-)
    {0.45,0.45},  'f_par', %0.45,           % 10    Fraction of PAR in incoming solar radiation (-)
    {5, 5},'lambda_i',                      % 11    PAR light attenuation coefficient for ice (m-1)
    {15,15}, 'lambda_s',                    % 12    PAR light attenuation coefficient for snow (m-1)
    {0.36,0.36}, 'F_sed_sld',               % 13    volume fraction of solids in sediment (= 1-porosity)
    {1,1},  'I_scV', %0.5 0.752,   % 14    scaling factor for inflow volume (-)
    {0,0}, 'I_scT',  %1,           % 15    adjusting delta (i.e not a multiplicative factor a delta T oc) for inflow temperature (-)
%     10, 'I_scC',
    {1,1}, 'I_scPOC',                        % 16    scaling factor for inflow concentration of POC (-)
    {18,0},  'I_scTP', %50,       % 17    scaling factor for inflow concentration of total P (-) TP 1.9
    {1,1}, 'I_scDOP',                        % 18    scaling factor for inflow concentration of diss. organic P (-)
%     10,'I_scChl',
    {1,1}, 'I_scDOC',                        % 19    scaling factor for inflow concentration of DOC  (-)
    {5,5}, 'I_scPOP',              % 20    scaling factor for inflow concentration of POP  (-)
    {0.84622,0.84622}, 'I_scO',             % 21    Scaling factor for inflow concentration of O2 (-)
    {1,1}, 'I_scDIC',                       % 22    Scaling factor for inflow concentration of DIC  (-)
    {6, 0},  'I_scNO3',       % 23    Scaling factor for inflow concentration of NO3 (-)
    {1,1},  'I_scNH4',          % 24    Scaling factor for inflow concentration of NH4 (-)
    {1,1},  'I_scSO4',                      % 25    Scaling factor for inflow concentration of SO4 (-)
    {1,1}, 'I_scFe2',                       % 26    Scaling factor for inflow concentration of Fe2 (-)
    {1,1},  'I_scCa2',                            % 27    Scaling factor for inflow concentration of Ca2 (-)
    {1,1},  'I_scpH',                         % 28    Scaling factor for inflow concentration of pH (-)
    {1,1},  'I_scCH4',                        % 29    Scaling factor for inflow concentration of CH4 (-)
    {600,600},  'I_scFe3',                      % 30    Scaling factor for inflow concentration of Fe3 (-)
    {0.001,0.001},  'I_scAl3',                 % 31    Scaling factor for inflow concentration of Al3 (-)
    {1, 1}, 'I_scFeS',                      % 32    Scaling factor for inflow concentration of FeS (-)
    {2.5,2.5},  'I_scCaCO3',                   % 33    Scaling factor for inflow concentration of CaCO3 (-)
    {1,1},  'I_scCH4g',                      % 34    Scaling factor for inflow concentration of CH4g (-)
    {3,0}, 'I_scSi',               %35,              % 35 TP Scaling factor for inflow concentration of Si (-)
    {1.0 , 1.0},    'C_solar',       % 36    scaling for global radiation and PAR
    {1.0 , 1.0},    'C_air_temp',    % 37    scaling for air temp
    {0.25,0.25},  'swa_b0',  %2.5,   % 1     non-PAR light attenuation coeff. (m-1)
    {0.25,0.25},  'swa_b1',  %1.05,  % 2     PAR light attenuation coeff. (m-1)
    {3.3e-07,3.3e-07}, 'S_res_epi',                 % 3     Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
    {3.3e-08,3.3e-08}, 'S_res_hypo',             % 4     Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
    {0.03,0.03}, 'H_sed',                        % 5     height of active sediment layer (m, wet mass)
    {15,15},'Psat_L',                            % 6     NOTE: NOT USED: Half saturation parameter for Langmuir isotherm
    {30,30}, 'Fm2ax_L',                           % 7     NOTE: NOT USED: Scaling parameter for Langmuir isotherm !!!!!!!!!!!!
    {0.1,0.1}, 'w_s',                             % 8     settling velocity for S (m day-1)
    {1,1}, 'Y_cp',                               % 9     NOTE: NOT USED:  yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)   1/55*112/1 = 1
    {0.0002,0.0002}, 'k_twty',                      % 10    NOTE: NOT USED: specific Chl a to P transformation rate (1/day) at 20 deg C
    {0,0}, 'dop_twty',                                % 11    NOTE: NOT USED: specific DOP to P transformation rate (day-1) at 20 deg C
    {0.01,0.01}, 'oc_DOC',                            % 12    Optical cross-section of DOC (m2/mg DOC)
    {0.1,0.1}, 'qy_DOC',                               % 13    Quantum yield (mg DOC degraded/mol quanta)
    {0.1,0.1}, 'k_BOD',                               % 14    NOTE: NOT USED: Organic decomposition rate (1/d)
    {5,5}, 'w_CH4',                                 % 15    Methane gas rising velocity (m/d)
    {1.047,1.047}, 'theta_bod',                     % 16    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T ? 10 °C
    {1.13,1.13}, 'theta_bod_ice',                    % 17    NOTE: NOT USED: Temperature adjustment coefficient for BOD, T < 10 °C
    {1,1}, 'theta_sod',                             % 18    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T ? 10 °C
    {1,1}, 'theta_sod_ice',                          % 19    NOTE: NOT USED: Temperature adjustment coefficient for SOD, T < 10 °C
    {4,4}, 'BOD_temp_switch',                        % 20    NOTE: NOT USED: Threshold for bod or bod_ice °C
    {7.5,7.5}, 'pH',                                  % 21    Lake water pH
    {2,2}, 'Q10_wc',                                % 22    Q10 for reactions of respiration
    {1,1}, 'wc_factor',                             % 23    Scaling factor for rates in WC
    {4.84970,4.84970}, 'T_ref_wc'                  % 24    Reference Temperature for rates
    };       

% TP test diatom, Green, cyano

% TP ('Diatom grazed', 'Diatom UN-grazed','Green/Other grazed','Green/Other UN-grazed', 'Cyano N-limited', 'Cyano N-fixer')
% algal_params = {
%     {0.000025, 0.000025, 0.0005, 0.0005, 0.00075, 0.00075}, 'PAR_sat',        % 1    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
%     {0.03,   0.03,     0.005, 0.005, 0.01, 0.01}, 'beta_chl',         % 2    Optical cross_section of chlorophyll (m2 mg-1)
%     {0.3,    0.3,      0.05, 0.05, 0.005, 0.005}, 'w_chl',              % 3    Settling velocity for Chl a (m day-1) % TP 0.05
%     {0.18,   0.13,     0.05, 0.025, 0.025, 0.0125}, 'm_twty',       % 4    Loss rate (1/day) at 20 deg C  % TP 0.23063 {0.15, 0.1, 0.1, 0.05, 0.05, 0.025}
%     {1.1,    1.1,      2 ,2, 1.8, 1.8}, 'g_twty',                   % 5    Specific growth rate (1/day) at 20 deg C
%     {5,      5,         5, 5, 5, 5}, 'P_half',                 % 6    Half saturation growth P level (mg/m3) - per species
%     {80,     80,        80, 80, 80, 0.1}, 'N_half',                 % 7    Half saturation growth N level (mg/m3)TP 1.8 check units - per species
%     {550,    550,       NaN, NaN, NaN, NaN}, 'Si_half',              % 8    Half saturation growth Si level (mg/m3) - per species
%     {1,      1,         1,   1, 1, 0}, 'N-Limited'                 % 9    Define if algae is N-Limited (1) or not (0) 
%     {0.5,    0.5,       0.4, 0.4, 0.1, 0.1}, 'I_scChl'                 % 10   Scaling factor for inflow concentration of Chl (-)
%     {1,      1,         0,   0, 0, 0}, 'Si-Limited'                % 11   Define if algae is Si-Limited (1) or not (0) 
% };

% TP ('Diatom grazed', 'Diatom UN-grazed','Green/Other grazed','Green/Other UN-grazed', 'Cyano N-limited', 'Cyano N-fixer')
algal_params = {
    {0.000005, 0.000005, 0.0005, 0.0005, 0.00075, 0.00075}, 'PAR_sat',        % 1    PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
    {0.03,   0.03,     0.005, 0.005, 0.01, 0.01}, 'beta_chl',                   % 2    Optical cross_section of chlorophyll (m2 mg-1)
    {0.25,    0.25,      0.05, 0.05, 0.0005, 0.0005}, 'w_chl',                      % 3    Settling velocity for Chl a (m day-1) % TP 0.05
    {0.1,   0.1,     0.025, 0.025, 0.025, 0.025}, 'm_twty',                   % 4    Loss rate (1/day) at 20 deg C  % TP 0.23063 {0.15, 0.1, 0.1, 0.05, 0.05, 0.025}
    {0.85,    0.85,      1.0 ,1.2, 1.3, 1.3}, 'g_twty',                           % 5    Specific growth rate (1/day) at 20 deg C
    {5,      5,         5, 5, 5, 5}, 'P_half',                                  % 6    Half saturation growth P level (mg/m3) - per species
    {80,     80,        80, 80, 10, 10}, 'N_half',                             % 7    Half saturation growth N level (mg/m3)TP 1.8 check units - per species
    {550,    550,       NaN, NaN, NaN, NaN}, 'Si_half',              % 8    Half saturation growth Si level (mg/m3) - per species
    {1,      1,         1,   1, 1, 0}, 'N-Limited'                 % 9    Define if algae is N-Limited (1) or not (0) 
    {0.167,    0.167,       0.167, 0.167, 0.167, 0.167}, 'I_scChl'                 % 10   Scaling factor for inflow concentration of Chl (-)
    {1,      1,         0,   0, 0, 0}, 'Si-Limited'                % 11   Define if algae is Si-Limited (1) or not (0) 
};

sediment_params = {
    1.0549e-01,   'k_Chl',                 % 1       % 1
    1.2624e-02,  'k_POP',                 % 2       % 1
    5.2341e-02, 'k_POC',                  % 3       % 0.01
    1.2941e-02,  'k_DOP',                 % 4       % 1
    8.7662e-02, 'k_DOC',                  % 5       % 1
    0.008, 'Km_O2',                 % 6       % Canavan, R. W (2006) rho=2.5
    0.01,  'Km_NO3',                % 7       % Canavan, R. W (2006) rho=2.5
    100/2.5,  'Km_Fe(OH)3',         % 8       % Cappellen, 1996  rho=2.5
    100/2.5,  'Km_FeOOH',           % 9       % Cappellen, 1996  rho=2.5
    1.5,  'Km_SO4',                 % 10       % Cappellen, 1996 rho=2.5
    0.001,'Km_oxao',                % 11       NOTE note used % the same as Km rho=2.5
    0.1,  'Km_amao',                % 12       NOTE note used, % the same as Km rho=2.5
    0.008, 'Kin_O2',                % 13       % the same as Km rho=2.5
    0.01,  'Kin_NO3',               % 14       % the same as Km rho=2.5
    0.08,   'Kin_FeOH3',         % 15       % the same as Km rho=2.5
    0.08,   'Kin_FeOOH',         % 16       % the same as Km rho=2.5
    20,    'k_amox',                % 17       % Canavan, R. W (2006)
    50000,  'k_Feox',                % 18       % Canavan, R. W (2006)
    0.1,   'k_Sdis',                % 19       %
    2500,  'k_Spre',                % 20       %
    3.3,   'k_FeS2pre',             % 21       % Canavan (2006)
    0.1,   'k_alum',                % 22   % NOTE not used in the model
    6.3601e+00,     'k_pdesorb_a',         % 23
    1.1171e+01,     'k_pdesorb_b',         % 24
    20000,  'k_fesox',              % 25        % R23 %Canava
    1000,   'k_fes2ox',             % 26        % R23 % Katsev (2013)
    8,      'k_tS_Fe',              % 27      % Cappellen (1996) in Canavan, R. W (2006) the reaction is different
    9600,  'Ks_FeS',                % 28      % Canavan, R. W (2006)
    0.001, 'k_Fe_dis',              % 29      % Canavan, R. W
    0.04,'k_Fe_pre',             % 30        % Katsev, R. W (2013)
    0.00037,  'k_apa_pre',          % 31
    0.37,     'k_apa_dis',          % 32
    5.95799315598138e-05,  'K_apa',         % 33      % linl.dat PHREEQC
    0.04,  'k_CaCO3_pre',        % 34      % Katsev (2013)
    0.05,  'k_CaCO3_dis',           % 35      % Katsev (2013)
    5.0e-09,  'K_CaCO3',               % 36      %
    180,  'k_FeCO3_pre',        % 37      % Cappellen (1996)
    0.25,     'k_FeCO3_dis',        % 38      % Cappellen (1996)
    3.98107170553497e-09,  'K_FeCO3',          % 39      % Cappellen (1996)
    0.00037,  'k_viv_pre',          % 40
    0.37,  'k_viv_dis',             % 41
    1.88929597786807e-05, 'K_viv',          % 42     % llnl.dat PHREEQC
    1.0e-06,  'k_oms',                 % 43
    10000,   'k_tsox',                % 44     % Canavan, R. W (2006)
    0.12, 'k_FeSpre',            % 45     % from "Non-steady state diagenesis of organic and inorganic sulfur in lake sediments Raoul-Marie Couture, Rachele Fischer b, Philippe Van Cappellen b, Charles Gobeil c
    10000000,   'k_ch4_o2',              % 46     % Canavan, R. W (2006)
    0.1,  'k_ch4_so4',             % 47     % Canavan, R. W (2006)
    0.0015,  'Kh_CH4',              % 48     % Henry cobstant M/atm
    1000,   'k_ch4_dis',             % 49
    32.5,     'w_CH4g',                % 50     % Rising velocity of methane
    0.034,  'Kh_CO2',               % 51     % Henry cobstant M/atm
    32.5,  'accel',                 % 52
    0.15*0.5,   'Kd_fe2',           % 53
    1.35,   'k_pdesorb_c',          % 54
    0.95,   'fi_in',                % 55
    0.85,   'fi_f',                 % 56
    0.5,    'X_b',                  % 57
    1,      'tortuosity',           % 58
    0.1,    'w',                    % 59
    301,    'n',                    % 60
    30,     'depth',                % 61
    7.2,   'alfa0',                % 62
    106,    'Cx1',                  % 63           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    16,     'Ny1',                  % 64           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz1',                  % 65           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    200,    'Cx2',                  % 66           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    20,     'Ny2',                  % 67           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Pz2',                  % 68           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    1,      'Cx3',                  % 69           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0.1,    'Ny3',                  % 70           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    0,      'Pz3',                  % 71           % OM composition, it also defines rates of reaction (lower number - slower the reaction)
    [3,30,30,30,25,10],     'effective_depth',     % 72           % depth below which the lake is affected by sediments, [m], if -1 (experimental) , then sediments below pycnocline
    30,     'effective_depth',     % 72           % depth below which the lake is affected by sediments, [m], if -1 (experimental) , then sediments below pycnocline
    192,     'n_ts',                 % 73           % (96 is the minimum, 48 for calibration) number of time steps during 1 day (fixed time step of MyLake) for chemical and sediment module (the modules should be in sync)
    0,      'pH algorithm',         % 74           % 0. Disabled pH=7 % 1. Phreeqc  % 2. Electro-neutrality Equation
    0.1,    'SO4 flux',             % 75           % default 0; flux of sulphate from bottom of the sediment. Custom boundary condition for Lake Vansjo only
    };

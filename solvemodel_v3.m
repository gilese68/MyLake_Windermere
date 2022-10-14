function [MyLake_results, sediment_results] = solvemodel_v3(M_start, M_stop, Initfile, Inputfile, Parafile, tanks_params, algal_params, varargin)

% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
%
% VERSION 1.2.1 (two phytoplankton groups are included; variable Cz denotes
% this second group now. Frazil ice included + some small bug-fixes and code 
% rearrangements. Using convection_v2.m code)
%
% Main module
% Code checked by TSA, xx.xx.200x
% Last modified by TSA, 21.08.2007

% Modified to include Fokema-module by Kai Rasmus. 16.5.2007
% Modified to include the latest Fokema module 30.12.2010 by PK
% New matrices: DOCzt1,DOCzt2,DOCzt3,Daily_BB1t,Daily_BB2t,Daily_BB3t,Daily_PBt

% New DIC variable 29.12.2010 (incl. inflow, convection, diffusion) by PK
% New O2 variable 10.2.2011 by PK

% Major Revisions JBA Consulting and Lancaster University 2020 to include
% multiple "tanks" and any number of algae

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs (to function)
%       M_start     : Model start date [year, month, day]
%       M_stop      : Model stop date [year, month, day]
%               + Input filenames and sheetnames

% Inputs (received from input module):
%		tt		: Solution time domain (day)
%       In.Z    : Depths read from initial profiles file (m)
%       In.Az   : Areas read from initial profiles file (m2)
%       In.Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In.Cz   : Initial chlorophyll (group 2) profile read from initial profiles file (-)
%       In.POCz   : Initial sedimenting tracer (or suspended inorganic matter) profile read from initial profiles file (kg m-3)
%       In.TPz  : Initial total P profile read from initial profiles file (incl. DOP & Chla & Cz) (mg m-3)
%       In.DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In.Chlz : Initial chlorophyll (group 1) profile read from initial profiles file (mg m-3)
%       In.DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In.DICz  : Initial DIC profile read from initial profiles file (mg m-3) (PK)
%       In.O2z  : Initial oxygen profile read from initial profiles file (mg m-3) (PK)
%       In.TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In.Chlz_sed : Initial chlorophyll profile  (groups 1+2) in the sediment compartments read from initial profiles file (mg m-3)
%       In.FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 23 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (23 * 2)
%       Bio_par_names   : Names for Bio_par

% Outputs (other than Inputs from input module):
%		Qst : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
%		Kzt	: Predicted vertical diffusion coefficient (tt * zz) (m2 d-1)
%		Tzt	: Predicted temperature profile (tt * zz) (deg C)
%		Czt	: Predicted chlorophyll (group 2) profile (tt * zz) (-)
%		POCzt	: Predicted passive sedimenting tracer (or suspended inorganic matter) profile (tt * zz) (kg m-3)=(g L-1)
%		Pzt	: Predicted dissolved inorganic phosphorus profile (tt * zz) (mg m-3)
%		Chlzt	    : Predicted chlorophyll (group 1) profileo (tt * zz) (mg m-3)
%		PPzt	    : Predicted particulate inorganic phosphorus profile (tt * zz) (mg m-3)
%		DOPzt	    : Predicted dissolved organic phosphorus profile (tt * zz) (mg m-3)
%		DOCzt	    : Predicted dissolved organic carbon (DOC) profile (tt * zz) (mg m-3)
%		DICzt	    : Predicted dissolved inorganic carbon (DIC) profile (tt * zz) (mg m-3) (PK)
%		CO2zt	    : Predicted dissolved carbon dioxide profile (tt * zz) (mg m-3) (PK)
%		O2zt	    : Predicted dissolved oxygen profile (tt * zz) (mg m-3) (PK)
%       O2_sat_rel  : Predicted relative oxygen saturation (PK)
%       O2_sat_abs  : Predicted absolute oxygen saturation (PK)
%		Qz_sed      : Predicted  sediment-water heat flux (tt * zz) (W m-2, normalised to lake surface area)
%       lambdazt    : Predicted average total light attenuation coefficient down to depth z (tt * zz) (m-1)
%       P3zt_sed    : Predicted P conc. in sediment for P (mg m-3), PP(mg kg-1 dry w.) and Chl (mg kg-1 dry w.) (tt * zz * 3)
%       P3zt_sed_sc : Predicted P source from sediment for P, PP and Chl (mg m-3 day-1) (tt * zz * 3)
%       His         : Ice information matrix ([Hi Hs Hsi Tice Tair rho_snow IceIndicator] * tt)
%       DoF, DoM    : Days of freezing and melting (model timestep number)
%       MixStat     : Temporary variables used in model testing, see code (N * tt)

% Fokema outputs
%       CDOMzt      : Coloured dissolved organic matter absorption m-1
%                   : (tt * zz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


warning off MATLAB:fzero:UndeterminedSyntax %suppressing a warning message

% These variables are still global and not transferred by functions
global ies80


%disp(['Running MyLake-Sediment from ' datestr(datenum(M_start)) ' to ' datestr(datenum(M_stop)) ' ...']);

% DO NOT CHANGE!
dt=1.0; %model time step = 1 day (DO NOT CHANGE!)
% DO NOT CHANGE!

%Read input data
[In, tt, Ice0, Wt, Inflw, Phys_par, Bio_par, Tanks_par] = modelinputs_v3(M_start, M_stop, Initfile, Inputfile, Parafile, dt, tanks_params);


load albedot1.mat; %load albedot1 table, in order to save execution time

mylake_params.floculation_switch = switches.floculation_switch; 
mylake_params.rate_estimator_switch = switches.rate_estimator_switch;
mylake_params.dt = dt;
% Algal parameters common to all tanks
mylake_params.PAR_sat = cell2mat(algal_params{1});
mylake_params.beta_chl = cell2mat(algal_params{2});
mylake_params.w_chl = cell2mat(algal_params{3});
mylake_params.m_twty = cell2mat(algal_params{4});
mylake_params.g_twty = cell2mat(algal_params{5});
mylake_params.P_half = cell2mat(algal_params{6});
mylake_params.N_half = cell2mat(algal_params{7});
mylake_params.Si_half = cell2mat(algal_params{8});
mylake_params.N_limited = cell2mat(algal_params{9});
mylake_params.I_scChl = cell2mat(algal_params{10});
mylake_params.Si_limited = cell2mat(algal_params{11}); % TP added Si limitation indicator

theta_m = exp(0.1*log(2));    % loss and growth rate parameter base, ~1.072
e_par = 240800;               % Average energy of PAR photons (J mol-1)
Tf = 0;                       % water freezing point temperature (deg C)
F_OM = 1e+6*0.012;            % mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)
Fstable = 655;                % Inactive P conc. in inorg. particles (mg/kg dw);

% Set up multiple tanks: depths, pars etc.
Nt = Tanks_par.Nt;
zz = [];
for t = 1:Nt
    tanks(t) = tank(Ice0(min(length(In), t), :), t, In(min(length(In), t)).Tz, switches.rate_estimator_switch,...
        theta_m, In(min(length(In), t)), Phys_par, Bio_par, Tanks_par.tank_areas, mylake_params);
    Nz_tanks(t) = tanks(t).Nz_tank;
    
    if length(tanks(t).zz_tank) > length(zz)
        zz = tanks(t).zz_tank;
    end
    tanks(t).dz = Phys_par(1,1); % TP added
end

% Initialise horizontal mixing
ml_horizontal_mix = horizontal_mix(Tanks_par.tank_advection_matrix, Tanks_par.tank_diffusion_matrix, ...
    Tanks_par.decline_parameter_advection, Tanks_par.decline_parameter_diffusion, Tanks_par.river_inflow_tank_id, ...
    Tanks_par.river_inflow_switch, zz, tanks, Nt, Tanks_par.river_in_coef, Tanks_par.river_out_coef);

% Allocate and initialise output data matrices
ml_output_matrices = output_matrices(tt, Nz_tanks, Nt, length(mylake_params.beta_chl));

%% >>>>>> Start of the time loop >>>>>>
for i = 1:length(tt);
     disp(num2str(i))
     global timesteptp
     timesteptp = i;
     
     
    for t = 1:Nt %loop on tanks 
        [tanks(t)] = tanks(t).run_time_step(i, tt, Wt, albedot1, ies80, dt, e_par);
    end %for t = 1:Nt
    
    %% Add river inflow    
    [tanks] = ml_horizontal_mix.tanks_river_inflow(Inflw, Tf, i, ies80, tanks, Nt, length(mylake_params.beta_chl));
    
    for t = 1:Nt %loop on tanks
        [tanks(t)] = tanks(t).run_time_step_2(Tf, M_start, i, ies80, dt, Wt);
    end
    
    %% Output MyLake matrices
    ml_output_matrices = ml_output_matrices.update_output_matrices(i, switches.photobleaching, tanks, Nt, dt, length(mylake_params.beta_chl));

end %for i = 1:length(tt)
% <<<<<< End of the time loop <<<<<<

stamp = datetime('now');

%Saving sediment values
if switches.sediment_module           % MATSEDLAB sediment module
    for t=1:Nt
        sediment_results.tanks(t).concentrations = tanks(t).sediment_concentrations_zt;
        sediment_results.tanks(t).z = tanks(t).sediment_params.x'; 
        sediment_results.tanks(t).Bioirrigation_fx_zt = tanks(t).sediment_bioirrigation_fluxes_zt;
        sediment_results.tanks(t).params = tanks(t).sediment_params;
        sediment_results.tanks(t).sediment_transport_fluxes = tanks(t).sediment_transport_fluxes_zt;
        sediment_results.tanks(t).rates = tanks(t).sediment_additional_results_zt.rates;
        sediment_results.tanks(t).mylake_params = tanks(t).mylake_params;
    end
    
    sediment_results.days = datenum(M_start):datenum(M_stop);
    sediment_results.m_start = M_start;
    sediment_results.m_stop = M_stop;
    sediment_results.date_of_run = stamp;
else
    sediment_results = {};
end

MyLake_results = ml_output_matrices.save_results(Nt);
MyLake_results.z = zz;
MyLake_results.days = datenum(M_start):datenum(M_stop);
MyLake_results.params = mylake_params;
MyLake_results.params.Phys_par = Phys_par;
MyLake_results.params.Bio_par = Bio_par;
MyLake_results.Inflw = Inflw;
MyLake_results.Wt = Wt;
MyLake_results.m_start = M_start;
MyLake_results.m_stop = M_stop;
MyLake_results.date_of_run = stamp;


end
%end of function


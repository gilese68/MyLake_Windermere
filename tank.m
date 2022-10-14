classdef tank
    %TANK regroups tank specific properties and methods
    %   One tank object should be initialised for each tank before the
    %   start of the time step. Calculation of variables during the time
    %   step is done in two phases, between which the horizontal mix (river
    %   inflow, advection and eddy diffusion between tanks) is added.
    
    properties
        tank_id
        Nz_tank             % Number of vertical layers in the tank
        Az_tank             % Tank areas at each layer
        Vz_tank             % Tank volume at each layer
        zz_tank             % Tank layer depth
        zm_tank
        profiles            % Structure containing species profiles
        Tzy_sed, Qz_sed
        S_resusp
        ice_snow
        sediment_concentrations
        sediment_params
        sediment_matrix_templates
        results
        mylake_params
        sediment_concentrations_zt
        sediment_transport_fluxes_zt
        sediment_bioirrigation_fluxes_zt
        sediment_additional_results_zt
        MyLake_results
        rho, tau, Qsw, Qlw, Qsl
        Inflow, lvlD
        Tprof_prev
        lambdaz_wtot
        lambdaz_wtot_avg, surfflux
        DayFrac
        H_sw_z
        Kz
        tank_vertical_diffusion
        tank_sediments
        I_sc
        N_algae         % Number of algal species
        dz % TP added
    end
    
    methods
        function obj = tank(Ice0, t, In_Tz, rate_estimator_switch, ...
                            theta_m, In, Phys_par, Bio_par, tank_areas, ...
                            mylake_params)
            %TANK Construct an instance of this class
            %   Set tank parameters: volume, area, depth, initial profiles,
            %   ice snow object, sediment initialisation
            
            global sed_par_file
            
            
            col = min(size(Phys_par, 2), t);
            
            % Unpack the more fixed parameter values from input array "Phys_par"
            dz = Phys_par(1, col); %grid step size (m)

            zm = In.Z(end); %max depth
            zz = [0:dz:zm-dz]'; %solution depth domain

            obj.mylake_params = mylake_params;
            obj.mylake_params.dz = Phys_par(1, col); % layer depth (m)
            obj.mylake_params.Kz_K1 = Phys_par(2, col); % open water diffusion parameter (-)
            obj.mylake_params.Kz_K1_ice = Phys_par(3, col); % under ice diffusion parameter (-)
            obj.mylake_params.Kz_N0 = Phys_par(4, col); % min. stability frequency (s-2)
            
            obj.tank_vertical_diffusion = vertical_diffusion(obj.mylake_params.Kz_K1, obj.mylake_params.Kz_K1_ice, obj.mylake_params.Kz_N0);

            obj.mylake_params.C_shelter = Phys_par(5, col); % wind shelter parameter (-)
            obj.mylake_params.C_solar = Phys_par(36, col); % TP add solar scaling parameter (-)
            obj.mylake_params.C_air_temp = Phys_par(37, col); % TP add solar scaling parameter (-)
            obj.mylake_params.lat = Phys_par(6, col); %latitude (decimal degrees)
            obj.mylake_params.lon = Phys_par(7, col); %longitude (decimal degrees)
            obj.mylake_params.alb_melt_ice = Phys_par(8, col);   %albedo of melting ice (-)
            obj.mylake_params.alb_melt_snow = Phys_par(9, col); %albedo of melting snow (-)
%             obj.mylake_params.PAR_sat = Phys_par(10, col);         %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1)
            obj.mylake_params.f_par = Phys_par(10, col);           %Fraction of PAR in incoming solar radiation (-)
%             obj.mylake_params.beta_chl = Phys_par(12, col);        %Optical cross_section of chlorophyll (m2 mg-1)
            obj.mylake_params.lambda_i = Phys_par(11, col);       %PAR light attenuation coefficient for ice (m-1)
            obj.mylake_params.lambda_s = Phys_par(12, col);       %PAR light attenuation coefficient for snow (m-1)
            obj.mylake_params.F_sed_sld = Phys_par(13, col);      %volume fraction of solids in sediment (= 1-porosity)
            obj.I_sc.V = Phys_par(14, col); %scaling factor for inflow volume (-)
            obj.I_sc.T = Phys_par(15, col); %scaling coefficient for inflow temperature (-)
%             obj.I_sc.C = Phys_par(16, col); %scaling factor for inflow concentration of C (-)
            obj.I_sc.POC = Phys_par(16, col); %scaling factor for inflow concentration of POC (-)
            obj.I_sc.TP = Phys_par(17, col); %scaling factor for inflow concentration of total P (-)
            obj.I_sc.DOP = Phys_par(18, col); %scaling factor for inflow concentration of diss. organic P (-)
%             obj.I_sc.Chl = Phys_par(20, col); %scaling factor for inflow concentration of Chl a (-)
            obj.I_sc.DOC = Phys_par(19, col); %scaling factor for inflow concentration of DOC  (-)
            obj.I_sc.POP = Phys_par(20, col); %scaling factor for inflow concentration of POP  (-)
            obj.I_sc.O = Phys_par(21, col); %scaling factor for inflow concentration of O2  (-)
            obj.I_sc.DIC = Phys_par(22, col);   %Scaling factor for inflow concentration of DIC  (-)
            obj.I_sc.NO3 = Phys_par(23, col);   %scaling factor for inflow concentration   (-)
            obj.I_sc.NH4 = Phys_par(24, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.SO4 = Phys_par(25, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.Fe2 = Phys_par(26, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.Ca2 = Phys_par(27, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.pH = Phys_par(28, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.CH4aq = Phys_par(29, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.Fe3 = Phys_par(30, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.Al3 = Phys_par(31, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.FeS = Phys_par(32, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.CaCO3 = Phys_par(33, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.CH4g = Phys_par(34, col);    %scaling factor for inflow concentration   (-)
            obj.I_sc.Chl = mylake_params.I_scChl; %scaling factor for inflow concentration of Chl (-)
            obj.I_sc.Si = Phys_par(35, col); %TP: %scaling factor for inflow concentration of Si (-)
            % Unpack the more site specific parameter values from input array "Bio_par"

            obj.mylake_params.swa_b0 = Bio_par(1, col); % non-PAR light attenuation coeff. (m-1)
            obj.mylake_params.swa_b1 = Bio_par(2, col); %  PAR light attenuation coeff. (m-1)
            obj.mylake_params.S_res_epi = Bio_par(3, col);      %Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
            obj.mylake_params.S_res_hypo = Bio_par(4, col);     %Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
            obj.mylake_params.H_sed = Bio_par(5, col);          %height of active sediment layer (m, wet mass)
            obj.mylake_params.Psat_L = Bio_par(6, col);           %Half saturation parameter for Langmuir isotherm
            obj.mylake_params.Fmax_L = Bio_par(7, col);    %Scaling parameter for Langmuir isotherm !!!!!!!!!!!!

            obj.mylake_params.w_s = Bio_par(8, col);              %settling velocity for S (m day-1)
            obj.mylake_params.Y_cp = Bio_par(9, col);            %yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)
            obj.mylake_params.k_twty = Bio_par(10, col);          %specific Chl a to P transformation rate (1/day) at 20 deg C
            obj.mylake_params.dop_twty = Bio_par(11, col);        %specific DOP to P transformation rate (day-1) at 20 deg C

            obj.mylake_params.oc_DOC = Bio_par(12, col);           %Optical cross-section of DOC (m2/mg DOC)
            obj.mylake_params.qy_DOC = Bio_par(13, col);           %Quantum yield (mg DOC degraded/mol quanta)
            %===========

            % Parameters for oxygen

            obj.mylake_params.k_BOD = Bio_par(14, col);            %Organic decomposition rate (1/d)
            obj.mylake_params.w_CH4 = Bio_par(15, col);            % Methane gas rising velocity (m d-1)
            obj.mylake_params.theta_bod = Bio_par(16, col);        %Temperature adjustment coefficient for BOD, T ? 10 °C
            obj.mylake_params.theta_bod_ice = Bio_par(17, col);    %Temperature adjustment coefficient for BOD, T < 10 °C
            obj.mylake_params.theta_sod = Bio_par(18, col);        %Temperature adjustment coefficient for SOD, T ? 10 °C
            obj.mylake_params.theta_sod_ice = Bio_par(19, col);    %Temperature adjustment coefficient for SOD, T < 10 °C
            obj.mylake_params.BOD_temp_switch = Bio_par(20, col);             %Threshold for bod or bod_ice °C
            obj.mylake_params.pH = Bio_par(21, col);               %Lake water pH
            % WC chemistry:
            obj.mylake_params.Q10 = Bio_par(22, col);
            obj.mylake_params.wc_factor = Bio_par(23, col);
            obj.mylake_params.T_ref = Bio_par(24, col);
            obj.mylake_params.theta_m = theta_m;

            % Sediments class initialisation
            obj.tank_sediments = sediments(obj.mylake_params.Fmax_L, obj.mylake_params.H_sed, obj.mylake_params.F_sed_sld);

            obj.N_algae = length(mylake_params.beta_chl);
            
            %% Initial profiles
            if ~isnan(tank_areas)
                Az = interp1(In.Z, tank_areas(:, t), zz);
            else
                Az = interp1(In.Z, In.Az, zz);
            end
            Vz = dz * (Az + [Az(2:end, :); zeros(1, 1)]) / 2;
            
            Nz = find(Vz>0, 1, 'last' );
            Az = Az(1:Nz);
            Vz = Vz(1:Nz);
            zz = zz(1:Nz);

            %% Passing MyLake parameters in Chemical module
            obj.mylake_params.dz = dz; 
            obj.mylake_params.zm = zm; 
            obj.mylake_params.zz = zz; 

            obj.mylake_params.I_scV = obj.I_sc.V; 
            obj.mylake_params.I_scT = obj.I_sc.T; 
%             obj.mylake_params.I_scC = obj.I_sc.C; 
            obj.mylake_params.I_scPOC = obj.I_sc.POC; 
            obj.mylake_params.I_scTP = obj.I_sc.TP; 
            obj.mylake_params.I_scDOP = obj.I_sc.DOP; 
%             obj.mylake_params.I_scChl = obj.I_sc.Chl; 
            obj.mylake_params.I_scDOC = obj.I_sc.DOC; 
            obj.mylake_params.I_scPOP = obj.I_sc.POP; 
            obj.mylake_params.I_scO = obj.I_sc.O; 
            obj.mylake_params.I_scDIC = obj.I_sc.DIC; 
            obj.mylake_params.I_scNO3 = obj.I_sc.NO3; 
            obj.mylake_params.I_scNH4 = obj.I_sc.NH4; 
            obj.mylake_params.I_scSO4 = obj.I_sc.SO4; 
            obj.mylake_params.I_scFe2 = obj.I_sc.Fe2; 
            obj.mylake_params.I_scCa2 = obj.I_sc.Ca2; 
            obj.mylake_params.I_scpH = obj.I_sc.pH; 
            obj.mylake_params.I_scCH4aq = obj.I_sc.CH4aq; 
            obj.mylake_params.I_scFe3 = obj.I_sc.Fe3; 
            obj.mylake_params.I_scAl3 = obj.I_sc.Al3; 
            obj.mylake_params.I_scFeS = obj.I_sc.FeS; 
            obj.mylake_params.I_scCaCO3 = obj.I_sc.CaCO3; 
            obj.mylake_params.I_scCH4g = obj.I_sc.CH4g; 

            obj.tank_id = t;
            obj.Nz_tank = Nz;
            obj.Az_tank = Az;
            obj.Vz_tank = Vz;
            obj.zz_tank = zz;
            obj.zm_tank = zz(end) + dz;
            
            obj.mylake_params.Az = Az;
            obj.mylake_params.Vz = Vz;
            obj.mylake_params.zz = zz;
            obj.mylake_params.zm = zz(end) + dz;
            
            obj.profiles = initialise_profiles(In, obj.mylake_params.zz, obj.mylake_params.dz, obj.N_algae);
            
            for j = 1:Nz
                obj.Tzy_sed(:,j) = interp1([0.2 10], [obj.profiles.Tz(j) 4], [0.2:0.2:2 2.5:0.5:10])';
            end

            obj.S_resusp = obj.mylake_params.S_res_hypo * ones(Nz,1); %hypolimnion resuspension assumed on the first time step

            % Initialise ice snow cover module
            obj.ice_snow = ice_snow_cover(Ice0, t, obj.mylake_params.lambda_i, obj.mylake_params.lambda_s);

            % ============ sediment module ============
            spec = '%f';
            for ii=1:6
                spec = strcat(spec, '%f');
            end
            
            f=fopen(sed_par_file);
            data = textscan(f,spec, 75,'Delimiter', '\t');
            fclose(f);
            
            sed_params = [data{1}, data{t+1}];
            temp_sed = tempname;
            dlmwrite(temp_sed, sed_params,'delimiter','\t');
            
            og_sed_par_file = sed_par_file;
            sed_par_file = temp_sed;
            
            % Allocation and initial sediment profiles concentrations and reading initial concentrations for sediment from file
            [obj.sediment_concentrations, obj.sediment_params, obj.sediment_matrix_templates] = sediment_init(obj.mylake_params.pH, obj.zm_tank, In_Tz(end), rate_estimator_switch);
            
            sed_par_file = og_sed_par_file;
            
        end
        
        function [obj] = run_time_step(obj, i, tt, Wt, albedot1, ies80, dt, e_par)
            %RUN_TIME_STEP First phase of the time step calculation.
            
            Nz = obj.Nz_tank;
            Az = obj.Az_tank;
            Vz = obj.Vz_tank;
            zz = obj.zz_tank;
            
            lat = obj.mylake_params.lat;
            lon = obj.mylake_params.lon;
            alb_melt_ice = obj.mylake_params.alb_melt_ice;
            alb_melt_snow = obj.mylake_params.alb_melt_snow;
            swa_b0 = obj.mylake_params.swa_b0;
            swa_b1 = obj.mylake_params.swa_b1;
            beta_chl = obj.mylake_params.beta_chl;
            dz = obj.mylake_params.dz;
            f_par = obj.mylake_params.f_par;
            PAR_sat = obj.mylake_params.PAR_sat;
            w_s = obj.mylake_params.w_s;
            w_CH4 = obj.mylake_params.w_CH4;
            w_chl = obj.mylake_params.w_chl;
            oc_DOC = obj.mylake_params.oc_DOC;
            qy_DOC = obj.mylake_params.qy_DOC;
            pH = obj.mylake_params.pH;
            C_shelter = obj.mylake_params.C_shelter;
            C_solar = obj.mylake_params.C_solar; %TP added to modify solar heating and PAR
            C_air_temp = obj.mylake_params.C_air_temp; %TP added to modify air temp
            
            
            
            
            %% Surface heat fluxes (W m-2), wind stress (N m-2) & daylight fraction (-), based on Air-Sea Toolbox
                       
            GRad = Wt(i,1) * C_solar;% TP Wt(i,1) is global radiation - apply factor to this for obj.tank_id
            ATemp = Wt(i,3) * C_air_temp;% TP Wt(3) is air temperature - apply factor to this for obj.tank_id
            [obj.Qsw, obj.Qlw, obj.Qsl, obj.tau, obj.DayFrac, DayFracHeating] = heatflux_v12(tt(i),  GRad, Wt(i,2), ATemp, Wt(i,4), Wt(i,5), Wt(i,6), obj.profiles.Tz(1), ...
                lat, lon, obj.ice_snow.WEQs, obj.ice_snow.Hi, alb_melt_ice, alb_melt_snow, albedot1);     %Qlw and Qsl are functions of Tz(1)

            %% Calculate total mean PAR and non-PAR light extinction coefficient in water (background + due to Chl a)
            [obj.lambdaz_wtot_avg, lambdaz_NP_wtot_avg] = calculate_average_light_extinction_coefficients(Nz, swa_b0, swa_b1, beta_chl, ...
                obj.profiles.Chlz, switches.selfshading_switch);

            obj.ice_snow = obj.ice_snow.set_attenuation_coefficient();

            obj.rho = polyval(ies80, max(0, obj.profiles.Tz)) + min(obj.profiles.Tz,0);  % Density (kg/m3)

            %% Sediment vertical heat flux, Q_sed
            % (averaged over the whole top area of the layer, although actually coming only from the "sides")
            [obj.Qz_sed, obj.Tzy_sed] = obj.tank_sediments.calculate_sediment_heat_flux(switches.sediment_heatflux_switch, obj.Tzy_sed, obj.profiles.Tz, dt, Nz, Az);

            Cw = 4.18e+6;	% Volumetric heat capacity of water (J K-1 m-3)

            %% Heat sources/sinks:
            %Total attenuation coefficient profile, two-band extinction, PAR & non-PAR
            Par_Attenuation = exp([0; -obj.lambdaz_wtot_avg] .* [zz; zz(end)+dz]);
            NonPar_Attenuation = exp([0; -lambdaz_NP_wtot_avg] .* [zz; zz(end)+dz]);

            Attenuation_z = (-f_par * diff([1; ([Az(2:end); 0] ./ Az) .* Par_Attenuation(2:end)]) + ...
                (-(1-f_par)) * diff([1; ([Az(2:end); 0] ./ Az) .* NonPar_Attenuation(2:end)])); %NEW (corrected 210807)

            [dT, obj.ice_snow, Qz] = calculate_heat_source(obj.ice_snow, obj.Qsw, obj.Qlw, obj.Qsl, Attenuation_z, DayFracHeating, Az, Vz, dt, obj.Qz_sed, Cw);

            obj.Tprof_prev = obj.profiles.Tz; %temperature profile at previous time step (for convection_v2.m)
            obj.profiles.Tz = obj.profiles.Tz + dT;        %Temperature change after daytime surface heatfluxes (or whole day in ice covered period)

            %% Convective mixing adjustment (mix successive layers until stable density profile)
            % and
            % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
            [obj.profiles] = convection_v2(obj.profiles, obj.Tprof_prev, Vz, Cw, f_par, obj.lambdaz_wtot_avg, zz, swa_b0, switches.tracer_switch, 1);
            obj.Tprof_prev = obj.profiles.Tz; %NEW!!! Update Tprof_prev
            
            if(obj.ice_snow.IceIndicator == 0)
                [obj.ice_snow, obj.profiles, obj.Qlw, obj.Qsl, obj.DayFrac, obj.Qsw, obj.tau] = calculate_nightime_heating(obj.ice_snow, tt, Wt, obj.profiles, lat, lon, ...
                alb_melt_ice, alb_melt_snow, albedot1, Az, Vz, Cw, obj.Tprof_prev, i, f_par, obj.Qz_sed, obj.Qlw, obj.Qsl, dt, ...
                obj.lambdaz_wtot_avg, zz, swa_b0, switches.tracer_switch, Qz, obj.mylake_params.C_solar);  
            end
            
            %% Vertical turbulent diffusion
            g = 9.81; % Gravity acceleration (m s-2)
            [Fi, obj.Kz] = obj.tank_vertical_diffusion.calculate_turbulent_vertical_diffusion(obj.profiles.Tz, zz, obj.ice_snow, Vz, Az, dz, dt, ies80, g);

            obj.profiles.Tz = Fi \ (obj.profiles.Tz);        %Solving new temperature profile (diffusion, sources/sinks already added to Tz above)

            %% Convective mixing adjustment (mix successive layers until stable density profile)
            % (don't allow temperature jumps over temperature of maximum density, no summer/autumn turnover here!)
            [obj.profiles] = convection_v2(obj.profiles, obj.Tprof_prev, Vz, Cw, f_par, obj.lambdaz_wtot_avg, zz, swa_b0, switches.tracer_switch, 0);

            % Atmospheric Deposition
            if switches.deposition_switch == 1
                obj.profiles.Pz = obj.profiles.Pz + (Deposition(i,5) ./ sum(Vz(1))) ; % qty added in mg to the top layer z
            end

            % Calculate again the total mean PAR light extinction coefficient in water (background + due to Chl a)
            obj.lambdaz_wtot_avg = zeros(Nz,1);
            lamm_t = zeros(Nz,obj.N_algae);

            if (switches.selfshading_switch == 1)
                 
                obj.lambdaz_wtot =  swa_b1 * ones(Nz,1);
                for ii = 1:obj.N_algae% TP loop around algal species
                  obj.lambdaz_wtot = obj.lambdaz_wtot + beta_chl(ii) * obj.profiles.Chlz{ii} ; 
                  lamm_t(:,ii) = beta_chl(ii) * obj.profiles.Chlz{ii};% TP temp variable for multiple algal species
                end
                
                for j=1:Nz 
                    obj.lambdaz_wtot_avg(j) = mean(swa_b1 * ones(j,1) +  sum(lamm_t(1:j,:)',1)');
                end
  
            else %constant with depth
                obj.lambdaz_wtot = swa_b1 * ones(Nz,1);
                obj.lambdaz_wtot_avg = swa_b1 * ones(Nz,1);
            end %if selfshading...

            [obj.H_sw_z, obj.profiles.PAR_z] = calculate_photosynthetic_active_radiation(Nz, obj.ice_snow, obj.DayFrac, e_par, f_par, obj.Qsw, ...
            obj.lambdaz_wtot_avg, zz, PAR_sat, C_solar);

            % NOTE: All reactions are moved in "rates" method below on 26.03.2017
            % Corrected on 26.06.2017: diffusion and advection for all species in WC
            [obj.profiles, obj.surfflux, fokema_variables] = calculate_advection_and_diffusion_in_profiles(obj.profiles, obj.Kz, Vz, Az, dz, dt, w_s, w_CH4, w_chl, tt, i, ...
                switches.photobleaching, obj.Qsw, zz, Fi, oc_DOC, qy_DOC, f_par, e_par, Attenuation_z, pH, obj.ice_snow, C_shelter, Wt);

            %Sediment-water exchange (DOP source neglected)
            %-porewater to water

            if switches.resuspension_enabled == 0
                sediments.ksw = 0;
            end
            
            if (switches.photobleaching == 1)
                obj.results.DOCz1 = fokema_variables.DOCz1_new; %Fokema-model DOC subpool 1
                obj.results.DOCz2 = fokema_variables.DOCz2_new; %Fokema-model DOC subpool 2
                obj.results.DOCz3 = fokema_variables.DOCz3_new; %Fokema-model DOC subpool 3
                obj.results.DOC1frac = fokema_variables.DOC1frac; %Fokema-model subpool 1 fraction
                obj.results.DOC2frac = fokema_variables.DOC2frac; %Fokema-model subpool 2 fraction
                obj.results.DOC3frac = fokema_variables.DOC3frac; %Fokema-model subpool 3 fraction
                obj.results.Daily_BB1 = fokema_variables.Daily_BB1; %Fokema-model subpool 1 daily bacterial decomposition
                obj.results.Daily_BB2 = fokema_variables.Daily_BB2; %Fokema-model subpool 2 daily bacterial decomposition
                obj.results.Daily_BB3 = fokema_variables.Daily_BB3; %Fokema-model subpool 3 daily bacterial decomposition
                obj.results.Daily_PB = fokema_variables.Daily_PB; %Fokema-model daily photobleaching
            end
        end
            
        function [obj] = run_time_step_2(obj, Tf, M_start, i, ies80, dt, Wt)
            % RUN_TIME_STEP_2 Second phase of the time step calculation.
            global timesteptp
            Cw = 4.18e+6;	% Volumetric heat capacity of water (J K-1 m-3)
            g = 9.81; % Gravity acceleration (m s-2)
            t = obj.tank_id;
            Az = obj.Az_tank;
            Vz = obj.Vz_tank;
            zz = obj.zz_tank;
            S_res_epi = obj.mylake_params.S_res_epi;
            S_res_hypo = obj.mylake_params.S_res_hypo;
            f_par = obj.mylake_params.f_par;
            swa_b0 = obj.mylake_params.swa_b0;
            C_shelter = obj.mylake_params.C_shelter;
            dz = obj.mylake_params.dz;
            pH = obj.mylake_params.pH;
            
            % Convective mixing adjustment (mix successive layers until stable density profile,  no summer/autumn turnover here!)
            [obj.profiles] = convection_v2(obj.profiles, obj.Tprof_prev, Vz, Cw, f_par, obj.lambdaz_wtot_avg, zz, swa_b0, switches.tracer_switch, 0);

            %% Ice cover module
            if (obj.ice_snow.IceIndicator == 0)

                [obj.profiles] = mixing_from_wind_energy(C_shelter, Az, Vz, obj.tau, obj.rho, dt, dz, zz, g, obj.profiles, ies80);

            else % ice cover module
                obj.ice_snow.XE_surf = (obj.profiles.Tz(1) - Tf) * Cw * dz; %Daily heat accumulation into the first water layer (J m-2)
                obj.profiles.Tz(1) = Tf; %Ensure that temperature of the first water layer is kept at freezing point

                if (Wt(i,3) < Tf) %if air temperature is below freezing
                    [obj.ice_snow, Hi_new, dWEQs, dWEQnews] = obj.ice_snow.ice_and_snow_formation(i, Tf, Wt);

                else %if air temperature is NOT below freezing
                    [obj.ice_snow, Hi_new, dWEQs, dWEQnews] = obj.ice_snow.ice_and_snow_melting(Tf, obj.Qsw, obj.Qlw, obj.Qsl);

                end %if air temperature is or isn't below freezing

                %% Update ice and snow thicknesses
                [obj.ice_snow] = obj.ice_snow.update_ice_and_snow_thicknesses(Hi_new, dWEQs);

                [obj.ice_snow] = obj.ice_snow.update_snow_density(dWEQnews, switches.snow_compaction_switch, Wt, Tf, i);

                [obj.ice_snow] = obj.ice_snow.reset_ice_snow(t, i, M_start);
                
            end %of ice cover module

            %== P-partitioning in water==
            % TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
            % [Pz, trash]=Ppart(POCz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
            % PPz=TIPz-Pz;

            %DIC-partitioning in water
            [obj.profiles.CO2z, obj.profiles.HCO3z, obj.profiles.CO3z, ~, ~, ~] = carbonequilibrium(obj.profiles.DICz, obj.profiles.Tz, pH);
            % Relative dissolved oxygen concentration
            [O2_sat_rel, O2_sat_abs] = relative_oxygen(obj.profiles.O2z, obj.profiles.Tz, Wt(i,5),dz);

            %% Initial freezing
            [obj.ice_snow, obj.profiles.Tz] = obj.ice_snow.initial_freezing(obj.profiles.Tz, Tf, Vz, Cw, Az, i, M_start);

            %% Calculate pycnocline depth
            pycno_thres = 0.05;  %treshold density gradient value (kg m-3 m-1)
            obj.rho = polyval(ies80,max(0, obj.profiles.Tz)) + min(obj.profiles.Tz,0);
            dRdz = [NaN; abs(diff(obj.rho))];
            di = find((dRdz < (pycno_thres * dz)) | isnan(dRdz));
            %dRdz(di)=NaN;
            %TCz = nansum(zz .* dRdz) ./ nansum(dRdz);
            dRdz(di) = 0; %modified for MATLAB version 7
            
            
            % TP remove old pycnocline code
            %TCz = sum(zz .* dRdz) ./ sum(dRdz);
            % Estimate thermocline
            TCz = FindThermoDepth(obj.rho,zz,0.1);
            % Estimate top of metalimnion
            TCz = FindMetaTop(abs(diff(obj.rho)),TCz,zz,0.1);
            
            if TCz < 0.51
                TCz = NaN;% To be consistent with MyLake change to NaN if TCz not identified
            end
%             
            
%             plot(diff(obj.rho))
%             plot(obj.profiles.Tz)
%             plot(obj.rho)
%             if timesteptp > 150
%                 'trev'
%             end
%             
            
            
            %vector with S_res_epi above, and S_res_hypo below the pycnocline
            inx = find(zz <= TCz);
            obj.S_resusp(inx) = S_res_epi;
            inx=find(zz > TCz);
            obj.S_resusp(inx) = S_res_hypo;

            if (obj.ice_snow.IceIndicator == 1)
                obj.S_resusp(:) = S_res_hypo;  %only hypolimnetic type of resuspension allowed under ice
            end

            if (isnan(TCz) & (obj.ice_snow.IceIndicator==0))
                obj.S_resusp(:) = S_res_epi;   %if no pycnocline and open water, resuspension allowed from top to bottom
            end


            %% WC chemistry:
            if any(isnan(obj.profiles.O2z)) | any(isnan(obj.profiles.Chlz{1})) | any(isnan(obj.profiles.DOCz)) | any(isnan(obj.profiles.NO3z)) | any(isnan(obj.profiles.Fe3z)) | any(isnan(obj.profiles.SO4z)) | any(isnan(obj.profiles.NH4z)) | any(isnan(obj.profiles.Fe2z)) | any(isnan(obj.profiles.H2Sz)) | any(isnan(obj.profiles.HSz)) | any(isnan(obj.profiles.Pz)) | any(isnan(obj.profiles.Al3z)) | any(isnan(obj.profiles.PPz)) | any(isnan(obj.profiles.Ca2z)) | any(isnan(obj.profiles.CO2z))
                error('NaN')
            end

            if any(isnan(obj.profiles.O2z)) | any(isnan(obj.profiles.Pz)) | any(isnan(obj.profiles.Fe2z)) | any(isnan(obj.profiles.NO3z)) | any(isnan(obj.profiles.NH4z))
                error('NaN')
            end

            if switches.wc_chemistry_module
                [obj.profiles, C_new, wc_rates_av, mylake_temp_results] = run_wc_chemistry_module(obj.profiles, obj.DayFrac, obj.lambdaz_wtot, obj.H_sw_z, obj.mylake_params, obj.sediment_params, dt, switches.wc_int_method, obj.N_algae);
            end

            if any(isnan(C_new))
                error('NaN')
            end

            if any(isnan(obj.profiles.O2z)) | any(isnan(obj.profiles.Pz)) | any(isnan(obj.profiles.Fe2z)) | any(isnan(obj.profiles.NO3z)) | any(isnan(obj.profiles.NH4z))
                error('NaN')
            end


            %% sediment module
            if switches.sediment_module
                    % Making cells of params for using during coupling
                [sediment_bioirrigation_fluxes, sediment_transport_fluxes, obj.sediment_concentrations, sediment_additional_results, obj.profiles] = ...
                    run_sediment_module(obj.profiles, mylake_temp_results, obj.sediment_params, obj.mylake_params, obj.sediment_concentrations, obj.sediment_matrix_templates, obj.N_algae);

                fields = fieldnames(obj.sediment_concentrations);
                for fd_idx = 1:numel(fields)
                    obj.sediment_concentrations_zt.(fields{fd_idx})(:,i) = obj.sediment_concentrations.(fields{fd_idx});
                end

                fields = fieldnames(sediment_transport_fluxes);
                for fd_idx = 1:numel(fields)
                    obj.sediment_transport_fluxes_zt.(fields{fd_idx})(:,i) = sediment_transport_fluxes.(fields{fd_idx});
                end

                fields = fieldnames(sediment_bioirrigation_fluxes);
                for fd_idx = 1:numel(fields)
                    obj.sediment_bioirrigation_fluxes_zt.(fields{fd_idx})(:,i) = sediment_bioirrigation_fluxes.(fields{fd_idx});
                end

                if obj.mylake_params.rate_estimator_switch

                    fields = fieldnames(sediment_additional_results.rates);
                    for fd_idx = 1:numel(fields)
                        obj.sediment_additional_results_zt.rates.(fields{fd_idx})(:,i) = sediment_additional_results.rates.(fields{fd_idx});
                    end

                    fields = fieldnames(wc_rates_av);
                    for fd_idx = 1:numel(fields)
                        obj.MyLake_results.rates.(fields{fd_idx})(:,i) = wc_rates_av.(fields{fd_idx});
                    end
                else
                    obj.sediment_additional_results_zt.rates = false;
                end

            end

            % Save variables for output matrices by tanks
            obj.results.Qs = [obj.Qsw obj.Qlw obj.Qsl]';
            obj.results.rho = obj.rho;
            obj.results.Kz = [0; obj.Kz];
            obj.results.Qz_sed = obj.Qz_sed;
            obj.results.lvlDz = obj.lvlD;
            obj.results.O2_sat_rel = O2_sat_rel;
            obj.results.O2_sat_abs = O2_sat_abs;
            obj.results.lambdaz_wtot_avg = obj.lambdaz_wtot_avg;
            obj.results.surfflux = obj.surfflux;
            obj.results.MixStat(3) = obj.lambdaz_wtot(2);%Inflow.DOC;
            obj.results.H_sw_z = diff([-obj.H_sw_z; zeros(1, size(obj.H_sw_z, 2))]);
            obj.results.His(1) = obj.ice_snow.Hi;
            obj.results.His(2) = (obj.ice_snow.rho_fw/obj.ice_snow.rho_snow)*obj.ice_snow.WEQs;
            obj.results.His(3) = obj.ice_snow.Hsi;
            obj.results.His(4) = obj.ice_snow.Tice;
            obj.results.His(5) = Wt(i,3);
            obj.results.His(6) = obj.ice_snow.rho_snow;
            obj.results.His(7) = obj.ice_snow.IceIndicator;
            obj.results.His(8) = obj.ice_snow.HFrazil; %NEW!!!
            obj.results.MixStat(1) = obj.Inflow.POC;
            obj.results.MixStat(2) = obj.Inflow.TP;
            Growth_bioz = 0; % NOTE: This part moved to reaction module
            Loss_bioz = 0; % NOTE: This part moved to reaction module
            obj.results.MixStat(4) = mean(Growth_bioz); %Only for chlorophyll group 1 (a)
            obj.results.MixStat(5) = mean(Loss_bioz);  %Only for chlorophyll group 1 (a)
            obj.results.MixStat(6) = obj.Inflow.V;
            if (obj.ice_snow.IceIndicator == 1)
                obj.results.MixStat(7:11) = NaN;
            else
                dum=interp1(zz,obj.profiles.Pz,[0:0.1:4]);
                obj.results.MixStat(7) = mean(dum); %diss-P conc. 0-4m in ice-free period

                dum=interp1(zz,obj.profiles.Chlz{1},[0:0.1:4]);  % TP - this is only species 1 but is a stat perhaps not used in our project??
                obj.results.MixStat(8) = mean(dum); %Chla conc. 0-4m in ice-free period

                dum=interp1(zz,obj.profiles.PPz,[0:0.1:4]);
                obj.results.MixStat(9) = mean(dum); %particulate inorg. P conc. 0-4m in ice-free period

                dum=interp1(zz,obj.profiles.DOPz,[0:0.1:4]);
                obj.results.MixStat(10) = mean(dum); %dissolved organic P conc. 0-4m in ice-free period

                dum=interp1(zz,obj.profiles.POCz,[0:0.1:4]);
                obj.results.MixStat(11) = mean(dum); %particulate matter conc. 0-4m in ice-free period
            end

            obj.results.MixStat(12) = TCz; %pycnocline depth

            obj.results.MixStat(13) = 1e-6*obj.Inflow.V*obj.Inflow.TP; %total P inflow (kg day-1)
            %         if (Inflow.V>Vz(1))
            %             disp('Large inflow!!')
            %         end
            obj.results.MixStat(14) = 1e-6*obj.Inflow.V*(obj.profiles.Pz(1) + obj.profiles.PPz(1) + obj.profiles.DOPz(1) + obj.profiles.Chlz{1}(1) + obj.profiles.Chlz{2}(1)); %total P outflow (kg day-1)

            if (obj.ice_snow.IceIndicator == 1)
                obj.results.MixStat(21) = NaN;
            else
                dum = interp1(zz,obj.profiles.Chlz{2}, [0:0.1:4]);
                obj.results.MixStat(21) = mean(dum); %Chl group 2 conc. 0-4m in ice-free period
            end
            % TP should be multi-species but appears to be unused
            Growth_bioz_2 = 0;% NOTE: This part moved to reaction module
            Loss_bioz_2 = 0;  % NOTE: This part moved to reaction module
            obj.results.MixStat(22) = mean(Growth_bioz_2); %For chlorophyll group 2
            obj.results.MixStat(23) = mean(Loss_bioz_2);  %For chlorophyll group 2
        end
    end
end


function [lambdaz_wtot_avg, lambdaz_NP_wtot_avg] = calculate_average_light_extinction_coefficients(Nz, swa_b0, swa_b1, beta_chl, ...
    Chlz, selfshading_switch)

% Calculate total mean PAR and non-PAR light extinction coefficient in water (background + due to Chl a)
lambdaz_wtot_avg = zeros(Nz,1);
lambdaz_NP_wtot_avg = zeros(Nz,1);

if (selfshading_switch == 1)
    % TP added for generalising to multiple algal species
    lambdaz_wtot =  swa_b1 * ones(Nz,1);
    lambdaz_wtot_NP =  swa_b0 * ones(Nz,1);
    for ii = 1:length(beta_chl)% TP loop around algal species  
        lamm_t(:,ii) = beta_chl(ii) * Chlz{ii};% TP temp variable for multiple algal species
        lamm_t_NP(:,ii) = beta_chl(ii) * Chlz{ii};% TP temp variable for multiple algal species
    end
    
    for j=1:Nz % TP modified for 1: N algae
        lambdaz_wtot_avg(j) = mean(swa_b1 * ones(j,1) +  sum(lamm_t(1:j,:)',1)');
        lambdaz_NP_wtot_avg(j) = mean(swa_b0 * ones(j,1) +  sum(lamm_t_NP(1:j,:)',1)');
    end
    
else %constant with depth
    lambdaz_wtot_avg = swa_b1 * ones(Nz,1);
    lambdaz_NP_wtot_avg = swa_b0 * ones(Nz,1);
end %if selfshading...
end    
    
function [dT, ice_snow, Qz] = calculate_heat_source(ice_snow, Qsw, Qlw, Qsl, Attenuation_z, DayFracHeating, Az, Vz, dt, Qz_sed, Cw)

    if(ice_snow.IceIndicator == 0)
        % 1) Vertical heating profile for open water periods (during daytime heating)
        Qz = (Qsw + ice_snow.XE_melt) * Attenuation_z; %(W m-2)
        Qz(1) = Qz(1) + DayFracHeating * (Qlw + Qsl); %surface layer heating
        ice_snow.XE_melt=0; %Reset
        dT = Az .* ((60*60*24*dt) * Qz + DayFracHeating * Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (daytime heating, ice melt, sediment);

        % === Frazil ice melting, NEW!!! === %
        [ice_snow, dT] = ice_snow.frazil_ice_melting(dT, Az, Cw, Vz);
        
    else
        % Vertical heating profile for ice-covered periods (both day- and nighttime)
        Qz = Qsw * ice_snow.IceSnowAttenuationCoefficient * Attenuation_z; %(W/m2)
        dT = Az .* ((60*60*24*dt) * Qz + Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (solar rad., sediment);
    end
end
    
function [ice_snow, profiles, Qlw, Qsl, DayFrac, Qsw, tau] = calculate_nightime_heating(ice_snow, tt, Wt, profiles, lat, lon, ...
    alb_melt_ice, alb_melt_snow, albedot1, Az, Vz, Cw, Tprof_prev, i, f_par, Qz_sed, Qlw, Qsl, dt, ...
    lambdaz_wtot_avg, zz, swa_b0, tracer_switch, Qz, C_solar)

    % 2) Vertical heating profile for open water periods (during nighttime heating)
    GRad = Wt(i,1) * C_solar;% TP Wt(i,1) is global radiation - apply factor to this for obj.tank_id
    [Qsw, Qlw_2, Qsl_2, tau, DayFrac, DayFracHeating] = heatflux_v12(tt(i), GRad, Wt(i,2), Wt(i,3), Wt(i,4), Wt(i,5), Wt(i,6), profiles.Tz(1), ...
        lat, lon, ice_snow.WEQs, ice_snow.Hi, alb_melt_ice, alb_melt_snow, albedot1); %Qlw and Qsl are functions of Tz(1)
    Qz(1) = (1 - DayFracHeating) * (Qlw_2 + Qsl_2); %surface layer heating
    Qz(2:end) = 0; %No other heating below surface layer
    dT = Az .* ((60*60*24*dt) * Qz + (1-DayFracHeating)*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (longwave & turbulent fluxes);

    % === NEW!!! frazil ice melting === %
    [ice_snow, dT] = ice_snow.frazil_ice_melting(dT, Az, Cw, Vz);

    profiles.Tz = profiles.Tz + dT;         %Temperature change after nighttime surface heatfluxes
    
    % Convective mixing adjustment (mix successive layers until stable density profile)
    % and
    % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
    [profiles] = convection_v2(profiles, Tprof_prev, Vz, Cw, f_par, lambdaz_wtot_avg, zz, swa_b0, tracer_switch, 1);
    
    Qlw = DayFracHeating * Qlw + (1 - DayFracHeating) * Qlw_2; %total amounts, only for output purposes
    Qsl = DayFracHeating * Qsl + (1 - DayFracHeating) * Qsl_2; %total amounts, only for output purposes
       
end

function [H_sw_z, PAR_z] = calculate_photosynthetic_active_radiation(Nz, ice_snow, DayFrac, e_par, f_par, Qsw, ...
    lambdaz_wtot_avg, zz, PAR_sat, C_solar)
    % e_par is average energy of PAR photons (J mol-1)
    % lambdaz_wtot_avg is average PAR light extinction coefficient
    %Photosynthetically Active Radiation 
    H_sw_z = NaN * zeros(Nz, length(PAR_sat));

    if ((ice_snow.IceIndicator == 0) && (DayFrac > 0))
        PAR_z = ((3/2) / (e_par * DayFrac)) * f_par * Qsw  * exp(-lambdaz_wtot_avg .* zz);
        %Irradiance at noon (mol m-2 s-1) at levels zz
    elseif ((ice_snow.IceIndicator == 1) && (DayFrac > 0))    %extra light attenuation due to ice and snow
        PAR_z = ((3/2) / (e_par * DayFrac)) * ice_snow.IceSnowAttenuationCoefficient * f_par *...
            Qsw  * exp(-lambdaz_wtot_avg .* zz);
    else
        PAR_z = zeros(Nz,1); %DayFrac==0, polar night
    end
    % =====
    PAR_z = PAR_z * C_solar + 1e-10; % TP modify for reduction in PAR associated with solar panel
    U_sw_z = repmat(PAR_z, 1, length(PAR_sat)) ./ PAR_sat; %scaled irradiance at levels zz
    inx_u = find(U_sw_z <= 1); %undersaturated
    inx_s = find(U_sw_z > 1);  %saturated

    H_sw_z(inx_u) = (2/3) * U_sw_z(inx_u);  %undersaturated

    dum_a = sqrt(U_sw_z);
    dum_b = sqrt(U_sw_z - 1);
    H_sw_z(inx_s) = (2/3) * U_sw_z(inx_s) + log((dum_a(inx_s) + dum_b(inx_s)) ./ (dum_a(inx_s) ...  %saturated
        - dum_b(inx_s))) - (2/3) * (U_sw_z(inx_s) + 2) .* (dum_b(inx_s) ./ dum_a(inx_s));
   
end
    
function [profiles, surfflux, fokema_variables] = calculate_advection_and_diffusion_in_profiles(profiles, Kz, Vz, Az, dz, dt, w_s, w_CH4, w_chl, tt, i, ...
    photobleaching, Qsw, zz, Fi, oc_DOC, qy_DOC, f_par, e_par, Attenuation_z, pH, ice_snow, C_shelter, Wt)

    Fi_ad_w_s = tridiag_HAD_v11([NaN; Kz], w_s, Vz, Az, dz, dt); %Tridiagonal matrix for advection and diffusion
    Fi_ad_w_CH4 = tridiag_HAD_v11([NaN; Kz], w_CH4, Vz, Az, dz, dt); %Tridiagonal matrix for advection and diffusion
   
    % TP generalise for multiple species
    for ii = 1:size(profiles.Chlz,2)
        Fi_ad_w_chl = tridiag_HAD_v11([NaN; Kz], w_chl(1), Vz, Az, dz, dt); %Tridiagonal matrix for advection and diffusion
        Fi_ad_w_chl_2 = tridiag_HAD_v11([NaN; Kz],w_chl(2), Vz, Az, dz, dt); %Tridiagonal matrix for advection and diffusion
        profiles.Chlz{ii} = Fi_ad_w_chl_2 \ (profiles.Chlz{ii});  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
        profiles.Chlz{ii} = Fi_ad_w_chl \ (profiles.Chlz{ii});  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level) 
    end
    
    profiles.Pz = Fi \ (profiles.Pz); %Solving new dissolved inorganic P profile (diffusion)
    profiles.DOPz = Fi \ (profiles.DOPz);
    profiles.NO3z = Fi \ profiles.NO3z;
    profiles.NH4z = Fi \ profiles.NH4z;
    profiles.SO4z = Fi \ profiles.SO4z;
    profiles.HSz = Fi \ profiles.HSz;
    profiles.H2Sz = Fi \ profiles.H2Sz;
    profiles.Fe2z = Fi \ profiles.Fe2z;
    profiles.Ca2z = Fi \ profiles.Ca2z;
    
    profiles.POCz = Fi_ad_w_s \ (profiles.POCz);           %Solving new suspended solids profile (advection + diffusion)
    profiles.PPz = Fi_ad_w_s \ (profiles.PPz);     %Solving new suspended particulate inorganic P profile (advection + diffusion)
    profiles.Fe3z = Fi_ad_w_s \ (profiles.Fe3z);
    profiles.Al3z = Fi_ad_w_s \ (profiles.Al3z);
    profiles.FeS = Fi_ad_w_s \ (profiles.FeSz);
    profiles.CaCO3z = Fi_ad_w_s \ (profiles.CaCO3z);
    profiles.POPz = Fi_ad_w_s \ (profiles.POPz);

    profiles.CH4gz = Fi_ad_w_CH4 \ (profiles.CH4gz);
    
    %Dissolved organic carbon
    % - current version
    Kd_old = 0;
    Theeta = 0;
    Currdate = datevec(tt(i)); %Date
    Date = Currdate(1,2); %Month number
    fokema_variables = struct;
    if (photobleaching == 1) %Fokema
        %[DOCz,Kd_new] = fokema(DOCz,Kd_old,Qsw,Tz,Theeta,Date,zz);
        profiles.DOCz1 = 0.0775 .* profiles.DOCz; %Subpools
        profiles.DOCz2 = 0.1486 .* profiles.DOCz;
        profiles.DOCz3 = 0.7739 .* profiles.DOCz;
        [fokema_variables.DOCz1_new, fokema_variables.DOCz2_new, fokema_variables.DOCz3_new, fokema_variables.DOC1frac, fokema_variables.DOC2frac, ...
            fokema_variables.DOC3frac, Kd_new, fokema_variables.Daily_BB1, fokema_variables.Daily_BB2, ...
            fokema_variables.Daily_BB3, fokema_variables.Daily_PB] = ...
            fokema_new(profiles.DOCz1, profiles.DOCz2, profiles.DOCz3, Kd_old, Qsw, profiles.Tz, Theeta, Date, zz);
        profiles.DOCz = DOCz1_new + DOCz2_new + DOCz3_new; %Total DOC
        profiles.DOCz = Fi \ profiles.DOCz; %Solving new dissolved organic C profile (diffusion)
        %DOCz1_new = Fi \ DOCz1_new; %Solving new dissolved organic C profile (diffusion)
        %DOCz2_new = Fi \ DOCz2_new; %Solving new dissolved organic C profile (diffusion)
        %DOCz3_new = Fi \ DOCz3_new; %Solving new dissolved organic C profile (diffusion)

    else %TSA model
        dDOC = -oc_DOC * qy_DOC * f_par * (1/e_par) * (60*60*24*dt) * Qsw * Attenuation_z; %photochemical degradation
        %[m2/mg_doc]*[mg_doc/mol_qnt]*[-]*[mol_qnt/J]*[s/day]*[J/s/m2]*[-] = [1/day]
        profiles.DOCz = Fi \ (profiles.DOCz + dDOC .* profiles.DOCz); %Solving new dissolved organic C profile (diffusion)
    end

    %Oxygen surface flux
    if(ice_snow.IceIndicator == 0)
        [profiles.O2z(1), O2flux , O2_eq, K0_O2] = oxygenflux(profiles.O2z(1), C_shelter^(1/3)*Wt(i,6), Wt(i,5), profiles.Tz(1), dz);
    else
        O2flux = 0;
    end

    profiles.O2z = Fi \ profiles.O2z; %Solving new dissolved oxygen profile (diffusion)
    profiles.DOCz = Fi \ profiles.DOCz;

    %Dissolved inorganic carbon
    %DIC partitioning in water
    [profiles.CO2z, profiles.HCO3z, profiles.CO3z, CO2frac, HCO3frac, CO3frac] = carbonequilibrium(profiles.DICz, profiles.Tz, pH);

    % CO2 production by degraded DOC
    % CO2z = max(0,CO2z + 1.375.*(-O2_diff));
    profiles.DICz = profiles.CO2z ./ CO2frac;
    %TC = Tz(1); %For monitoring only

    %Carbon dioxide surface flux
    if(ice_snow.IceIndicator==0)
        [profiles.CO2z(1), surfflux, CO2_eq, K0, CO2_ppm] = carbondioxideflux(profiles.CO2z(1), C_shelter^(1/3)*Wt(i,6), Wt(i,5), profiles.Tz(1), dz, tt(i));
        profiles.DICz(1) = profiles.CO2z(1) / CO2frac(1);
    else
        surfflux=0;
    end

    profiles.DICz = Fi \ profiles.DICz; %Solving new DIC profile (diffusion)

end

%=== ADVECTIVE-DIFFUSIVE EQUATION ===
function [Fi] = tridiag_HAD_v11(Kz,U,Vz,Az,dz,dt)

if (U<0)
    error('only positive (downward) velocities allowed')
end

if (U==0)
    U=eps; %set Vz next to nothing (=2.2204e-016) in order to avoid division by zero
end

Nz=length(Vz); %number of grid points/layers

theta=U*(dt/dz);

az = theta.*(1 + (1./(exp( (U*Vz)./(Kz.*Az) ) - 1)));                   %coefficient for i-1
cz = theta./(exp( (U*Vz)./([Kz(2:end); NaN].*[Az(2:end); NaN]) ) - 1);  %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1) = 1 + theta + cz(1);

%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function
end

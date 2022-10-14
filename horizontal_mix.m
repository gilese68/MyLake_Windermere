classdef horizontal_mix
    %HORIZONTAL_MIX calculates the river inflow into the tanks and the
    %advection and eddy diffusion between tanks
    %   Advection and eddy diffusion between tanks is defined by two input
    %   matrices. Exponential decline is applied vertically, defined by the
    %   decline parameter. River inflow can be calculated with two methods: 
    %   depending on layer density or using the decline profile
    
    properties
        advection_matrix        	% Advection coefficients between tanks
        diffusion_matrix            % Eddy diffusion coefficients between tanks
        decline_parameter_advection % Parameter for exponential decline function
        decline_profile_advection   % Decline profile calculated with decline_parameter_advection
        decline_parameter_diffusion % Parameter for exponential decline function
        decline_profile_diffusion   % Decline profile calculated with decline_parameter_advection
        tanks_decline_profile_advection       % Decline profiles for each couple of tanks, depends on tanks depth
        tanks_decline_profile_diffusion       % Decline profiles for each couple of tanks, depends on tanks depth
        river_inflow_tank_id        % ID of the tank with river inflow
        river_inflow_switch         % Switch to choose between the two method of river inflow (by density or decline profile)
        check_volume
        Nz_tanks
        Vz_tanks
        river_in_coef               % River inflow coefficient per tank
        river_out_coef              % River outflow coefficient per tank
    end
    
    methods
        function obj = horizontal_mix(advection_matrix, diffusion_matrix, decline_parameter_advection, decline_parameter_diffusion, ...
                river_inflow_tank_id, river_inflow_switch, zz, tanks, Nt, river_in_coef, river_out_coef)
            %HORIZONTAL_MIX Construct an instance of this class
            %   Initilise properties using user inputs
            obj.advection_matrix = advection_matrix;
            obj.diffusion_matrix = diffusion_matrix;
            obj.decline_parameter_advection = decline_parameter_advection;
            obj.decline_parameter_diffusion = decline_parameter_diffusion;
            obj.river_inflow_tank_id = river_inflow_tank_id;
            obj.river_inflow_switch = river_inflow_switch;
            obj.decline_profile_advection = exp(-zz * decline_parameter_advection);
            obj.decline_profile_diffusion = exp(-zz * decline_parameter_diffusion);
            obj.river_in_coef = river_in_coef;
            obj.river_out_coef = river_out_coef;
            obj.Nz_tanks = [tanks.Nz_tank];
            obj.Vz_tanks = zeros(max(obj.Nz_tanks), Nt);
            for t = 1:Nt
               obj.Vz_tanks(1:obj.Nz_tanks(t), t) = tanks(t).Vz_tank;
            end
            obj.tanks_decline_profile_advection = decline_profiles_advection(obj.decline_profile_advection, obj.Nz_tanks);
            obj.tanks_decline_profile_diffusion = decline_profiles_diffusion(obj.decline_profile_diffusion, obj.Nz_tanks);
            check_advection_coef = sum(obj.river_in_coef) - sum(obj.river_out_coef);
            for i=1:Nt
                for j=1:Nt
                    if i>j
                        check_advection_coef = check_advection_coef - obj.advection_matrix(i, j, 1);
                    else
                        check_advection_coef = check_advection_coef + obj.advection_matrix(j, i, 1);
                    end
                end
            end

            if check_advection_coef~=0
                disp('Warning: advection coefficients in the advection matrix and the river in and out coefficients do not ensure volume conservation - see code in horizontal_mix.m');
            end
            
        end
        
        function [tanks] = tanks_river_inflow(obj, Inflw, Tf, i, ...
                    ies80, tanks, Nt, N_algae)
            %TANK_RIVER_INFLOW calculates river inflow, advection and
            %diffusion between tanks
            %   Store current profiles, loop on tanks to update profiles
            %   with river inflow and advection between tanks. 
            %%%%%%%%%%% Limitation %%%%%%
            %   if river inflow method is by layer density, no flow can
            %   come from other tanks to the tank with river inflow.
                        
            obj.check_volume = zeros(max(obj.Nz_tanks), 1);

            
            
            for t = 1:Nt
                % Store profiles
                prev_profiles(t) = tanks(t).profiles;
                tanks(t).Inflow = get_inflow(tanks(t).I_sc, Inflw, i, Tf, obj.river_inflow_switch, N_algae);
            end
         
            for t = 1:Nt
                if obj.river_inflow_tank_id == t && obj.river_inflow_switch == 1
                    [tanks(t)] = river_inflow(Inflw, tanks(t).I_sc, Tf, i, ...
                    ies80, tanks(t).zz_tank, tanks(t).dz, tanks(t), prev_profiles(t));
                    for tt = 1:Nt
                        prev_profiles(tt) = tanks(tt).profiles;
                    end
                     %obj.check_volume(1) = obj.check_volume(1) + tanks(t).Inflow.V;
                else
                    tanks(t).lvlD = NaN;                   
                end
                [tanks(t).profiles, obj] = obj.advection_diffusion(tanks(t).profiles, prev_profiles, ...
                    t, tanks(t).Vz_tank, tanks(t).Nz_tank, tanks(t).Inflow, Nt, i);
            end
            if i==1
                disp(['Conservation of volume check per time step: ' num2str(sum(obj.check_volume)) 'm^3, namely ' num2str(sum(obj.check_volume)/sum(obj.Vz_tanks, 'all') * 100) '% of total volume.']);
                if sum(obj.check_volume) > 1 % & obj.river_inflow_switch == 2 
                   'You have a flow volume error - check advection and diffusion matricies' 
                   crash
                end
            end
       
        end
        
        function [profiles, obj] = advection_diffusion(obj, profiles, tanks_prev_profiles, t, Vz, Nz, Inflow, Nt, i)
           %ADVECTION_DIFFUSION calculates advection and diffusion between tanks, and river inflow
           %if the decline profile method has been chosen
            tanks_prev_profiles = resize(tanks_prev_profiles, obj.Nz_tanks, Nz);
            tank_decline_profile_advection = squeeze(obj.tanks_decline_profile_advection(t, :, 1:Nz));
            tank_decline_profile_diffusion = squeeze(obj.tanks_decline_profile_diffusion(t, :, 1:Nz));
            
            fn = fieldnames(profiles);
            for k = 1:numel(fn)
                if strcmp(fn{k}, 'Chlz') == 0
                    if isfield(Inflow, fn{k}(1:end-1))
                        advection_inflow = Inflow.(fn{k}(1:end-1));
                    else
                        advection_inflow = NaN;
                    end

                    if size(obj.advection_matrix, 3) == 1 %This is the case as can also be specified as a timeseries
                        [advection_in, advection_out, V_in_ad] = obj.advective_terms(advection_inflow, tanks_prev_profiles, Nz, t, fn{k}, tank_decline_profile_advection, Nt, 1, 0,obj.river_inflow_switch);
                    else
                        [advection_in, advection_out, V_in_ad] = obj.advective_terms(advection_inflow, tanks_prev_profiles, Nz, t, fn{k}, tank_decline_profile_advection, Nt, i, 0,obj.river_inflow_switch);
                    end

                    if size(obj.diffusion_matrix, 3) == 1
                        [diffusion_in, diffusion_out, V_in_diff] = obj.diffusive_terms(tanks_prev_profiles, Nz, t, fn{k}, tank_decline_profile_diffusion, Nt, 1, 0);
                    else
                        [diffusion_in, diffusion_out, V_in_diff] = obj.diffusive_terms(tanks_prev_profiles, Nz, t, fn{k}, tank_decline_profile_diffusion, Nt, i, 0);
                    end

                    
                    %TP Added this switch otherwise the volumes assocuated with tank concentrations are wrong
                    if obj.river_inflow_switch == 1
                        % Changes needed as not originally set up for the
                        % original V2 Mylake function which incorporates
                        % the river flow into the a given layer determined
                        % by the relative density
                        tank_out_term = Vz .* tanks_prev_profiles(t).(fn{k})(1:Nz) - diffusion_out .* Vz;
                        for t2 = 1:Nt % TP this is needed to separate pure advective fluxes from river inputs and outputs which
                            % are combined in obj.advective_terms
                            if size(obj.diffusion_matrix, 3) == 1 % i.e. not a timeseries of advection matrices
                                tank_out_term = tank_out_term - obj.advection_matrix(t, t2, 1) * tank_decline_profile_advection(t) .* Vz;
                            else
                                tank_out_term = tank_out_term - obj.advection_matrix(t, t2, i) * tank_decline_profile_advection(t) .* Vz;
                            end
                        end
                    elseif obj.river_inflow_switch > 1
                        tank_out_term = max(0, (Vz - (advection_out * Inflow.V  + diffusion_out .* Vz))) .* tanks_prev_profiles(t).(fn{k})(1:Nz);
                    end

                    % TP Modified this with a switch as otherwise putting mass in twice
                    if obj.river_inflow_switch == 1
                        profiles.(fn{k}) = diffusion_in + tank_out_term;
                    elseif obj.river_inflow_switch > 1
                        % Advection in + diffusion in + tank out term which is actually original mass in tank???
                        profiles.(fn{k}) = advection_in * Inflow.V + diffusion_in + tank_out_term; 
                    end
                    
                    % volume only
                    V_in = V_in_ad * Inflow.V + V_in_diff;
                    % volume only
                    V_out = advection_out * Inflow.V  + diffusion_out .* Vz;
                    
                    % back to concentreation by taking volume in and out in account:
                    profiles.(fn{k}) = profiles.(fn{k}) ./ (max(0, Vz - V_out) + V_in);

                    if k==1
                        obj.check_volume(1:Nz) = obj.check_volume(1:Nz) + V_in - V_out;
                    end
                else
                    for l=1:length(profiles.Chlz)
                        if isfield(Inflow, 'Chl') 
                            advection_inflow = Inflow.Chl{l};
                        else
                            advection_inflow = NaN;
                        end

                        if size(obj.advection_matrix, 3) == 1
                            [advection_in, advection_out, V_in_ad] = obj.advective_terms(advection_inflow, tanks_prev_profiles, Nz, t, 'Chlz', tank_decline_profile_advection, Nt, 1, l,obj.river_inflow_switch);
                        else
                            [advection_in, advection_out, V_in_ad] = obj.advective_terms(advection_inflow, tanks_prev_profiles, Nz, t, 'Chlz', tank_decline_profile_advection, Nt, i, l,obj.river_inflow_switch);
                        end

                        if size(obj.diffusion_matrix, 3) == 1
                            [diffusion_in, diffusion_out, V_in_diff] = obj.diffusive_terms(tanks_prev_profiles, Nz, t, 'Chlz', tank_decline_profile_diffusion, Nt, 1, l);
                        else
                            [diffusion_in, diffusion_out, V_in_diff] = obj.diffusive_terms(tanks_prev_profiles, Nz, t, 'Chlz', tank_decline_profile_diffusion, Nt, i, l);
                        end


                        tank_out_term = max(0, (Vz - (advection_out * Inflow.V  + diffusion_out .* Vz))) .* tanks_prev_profiles(t).Chlz{l}(1:Nz);
                        profiles.Chlz{l} = advection_in * Inflow.V + diffusion_in + tank_out_term;

                        V_in = V_in_ad * Inflow.V + V_in_diff;
                        V_out = advection_out * Inflow.V  + diffusion_out .* Vz;
                        % Taking volume in and out in account:
                        profiles.Chlz{l} = profiles.Chlz{l} ./ (max(0, Vz - V_out) + V_in);
                        % Tank volume only:
        %                profiles.(fn{k}) = profiles.(fn{k}) ./ Vz;

                        if k==1
                            obj.check_volume(1:Nz) = obj.check_volume(1:Nz) + V_in - V_out;
                        end
                    end
                end
            end
        end
                        
        function [in, out, V_in] = advective_terms(obj, Inflow, tanks_prev_profiles, Nz, t, fn, tank_decline_profile, Nt, i, l,Rivsw)
            %ADVECTIVE_TERMS calculates in, out and volume terms for the update of profiles by advection due to river inflow
            in = 0;
            out = 0;
            V_in = 0;
            if ~isnan(Inflow)
                for t2 = 1:Nt
                    if l ==0
                        in = in + obj.advection_matrix(t2, t, i) * tanks_prev_profiles(t2).(fn)(1:Nz) .* tank_decline_profile(t2, :)';
                    else
                        in = in + obj.advection_matrix(t2, t, i) * tanks_prev_profiles(t2).(fn){l}(1:Nz) .* tank_decline_profile(t2, :)';
                    end
                    out = out + obj.advection_matrix(t, t2, i) .* tank_decline_profile(t2, :)';
                    V_in = V_in + obj.advection_matrix(t2, t, i) * tank_decline_profile(t2, :)';
                end

                if Rivsw == 2 % TP river inputs if NOT using density-dependent inputs
                    in = in + Inflow * obj.river_in_coef(t) .* tank_decline_profile(t, :)';
                    V_in = V_in + obj.river_in_coef(t) * tank_decline_profile(t, :)';
                    out = out + obj.river_out_coef(t) * tank_decline_profile(t, :)';
                end
                
                if Rivsw == 3
                    %TP - inverting the advection profile but only for river inflows
                    inverse_tank_decline_profile = tank_decline_profile(t, Nz:-1:1)';
                    in = in + Inflow * obj.river_in_coef(t) .* inverse_tank_decline_profile;
                    V_in = V_in + obj.river_in_coef(t) * inverse_tank_decline_profile;
                    out = out + obj.river_out_coef(t) * inverse_tank_decline_profile;
                end
            end
        end
        
        function [in, out, V_in] = diffusive_terms(obj, tanks_prev_profiles, Nz, t, fn, tank_decline_profile, Nt, i, l)
            %DIFFUSIVE_TERMS calculates in, out and volume terms for the
            %update of profiles by eddy diffusion
            in = 0;
            out = 0;
            V_in = 0;
            for t2 = 1:Nt
                if l==0
                    in = in + obj.diffusion_matrix(t2, t, i) * obj.Vz_tanks(1:Nz, t2) .* tank_decline_profile(t2, :)' .* tanks_prev_profiles(t2).(fn)(1:Nz);
                else
                    in = in + obj.diffusion_matrix(t2, t, i) * obj.Vz_tanks(1:Nz, t2) .* tank_decline_profile(t2, :)' .* tanks_prev_profiles(t2).(fn){l}(1:Nz);
                end
                out = out + obj.diffusion_matrix(t, t2, i) .* tank_decline_profile(t2, :)';
                V_in = V_in + obj.diffusion_matrix(t2, t, i) * obj.Vz_tanks(1:Nz, t2) .* tank_decline_profile(t2, :)';
            end
        end
    end
end

function [profiles] = resize(profiles, Nz_tanks, Nz)
% RESIZE all profiles are resized to tank profile size Nz
    fn = fieldnames(profiles);
    for t=1:numel(Nz_tanks)
        if Nz_tanks(t) > Nz
            for k=1:numel(fn)
                if strcmp(fn{k},'Chlz') == 0
                    profiles(t).(fn{k}) = profiles(t).(fn{k})(1:Nz);
                else
                    for l=1:length(profiles(t).Chlz)
                        profiles(t).Chlz{l} = profiles(t).Chlz{l}(1:Nz);
                    end
                end
            end
        elseif Nz_tanks(t) < Nz
            for k=1:numel(fn)
                if strcmp(fn{k},'Chlz') == 0
                    profiles(t).(fn{k})(Nz) = 0;
                else
                    for l=1:length(profiles(t).Chlz)
                        profiles(t).Chlz{l}(Nz) = 0;
                    end
                end
           end
        end
    end
end

function [profiles] = decline_profiles_advection(profile, Nz_tanks)
%DECLINE_PROFILES calculates the decline profiles between two tanks
    profiles = zeros(length(Nz_tanks), length(Nz_tanks), max(Nz_tanks));
    for i=1:length(Nz_tanks)
        for j=1:length(Nz_tanks)
            min_Nz = min(Nz_tanks(i), Nz_tanks(j));
            profiles(i, j, 1:min_Nz) = profile(1:min_Nz) / sum(profile(1:min_Nz));
        end
    end
end

function [profiles] = decline_profiles_diffusion(profile, Nz_tanks)
%DECLINE_PROFILES calculates the decline profiles between two tanks
    profiles = zeros(length(Nz_tanks), length(Nz_tanks), max(Nz_tanks));
    for i=1:length(Nz_tanks)
        for j=1:length(Nz_tanks)
            min_Nz = min(Nz_tanks(i), Nz_tanks(j));
            profiles(i, j, 1:min_Nz) = profile(1:min_Nz);
        end
    end
end

function [Inflow] = get_inflow(I_sc, Inflw, i, Tf, river_inflow_switch, N_algae)
%GET_INFLOW calculates all inflow for time step i

    if river_inflow_switch ~= 0
        Inflow.V = I_sc.V * Inflw(i,1); % (scaled) inflow rate
        Inflow.T = I_sc.T + Inflw(i,2); %(adjusted) inflow temperature
        if (Inflow.T < Tf) %negative temperatures changed to Tf
            Inflow.T = Tf;
        end
        Inflow.POC = I_sc.POC * Inflw(i,3); %(scaled) inflow POC concentration
        Inflow.TP = I_sc.TP * Inflw(i,4); %(scaled) inflow TP concentration (incl. DOP & Chla)
        Inflow.DOP = I_sc.DOP * Inflw(i,5); %(scaled) inflow DOP concentration
        Inflow.DOC = I_sc.DOC * Inflw(i,6); %(scaled) inflow DOC concentration
        Inflow.DIC = I_sc.DIC * Inflw(i,7); %(scaled) inflow DIC concentration
        Inflow.O2 = I_sc.O * Inflw(i,8); %(scaled) inflow O2 concentration
        Inflow.NO3 = I_sc.NO3 * Inflw(i,9);
        Inflow.NH4 = I_sc.NH4 * Inflw(i,10);
        Inflow.SO4 = I_sc.SO4 * Inflw(i,11);
        Inflow.Fe2 = I_sc.Fe2 * Inflw(i,12);
        Inflow.Ca2 = I_sc.Ca2 * Inflw(i,13);
        Inflow.pH = I_sc.pH * Inflw(i,14);
        Inflow.CH4aq = I_sc.CH4aq * Inflw(i,15);
        Inflow.Fe3 = I_sc.Fe3 * Inflw(i,16);
        Inflow.Al3 = I_sc.Al3 * Inflw(i,17);
        Inflow.FeS = I_sc.FeS * Inflw(i,18);
        Inflow.CaCO3 = I_sc.CaCO3 * Inflw(i,19);
        Inflow.CH4g = I_sc.CH4g * Inflw(i,20);
        Inflow.POP = I_sc.POP * Inflw(i,21);
        % TP add in Silica inputs and CHl species inputs at end (1:N_algae)
        Inflow.Si = I_sc.Si * Inflw(i,22);
        for k=1:N_algae
            Inflow.Chl{k} = I_sc.Chl(k) * Inflw(i,22+k); %(scaled) inflow Chl a concentration
        end
        
        Inflow.P = (Inflow.TP - Inflow.DOP - Inflow.POP)/2;
        Inflow.P = Inflow.P .* (Inflow.P > 0);
        Inflow.PP = (Inflow.TP - Inflow.DOP - Inflow.POP)/2; %;
        Inflow.PP = Inflow.PP .* (Inflow.PP > 0);
        
    else
        Inflow.V = 0; % (scaled) inflow rate
        Inflow.T = NaN; %(adjusted) inflow temperature
        Inflow.C = NaN; %(scaled) inflow C concentration
        Inflow.POC = NaN; %(scaled) inflow S concentration
        Inflow.TP = NaN; %(scaled) inflow TP concentration (incl. DOP & Chla)
        Inflow.DOP = NaN; %(scaled) inflow DOP concentration
        Inflow.Chl = NaN; %(scaled) inflow Chl a concentration
        Inflow.DOC = NaN; %(scaled) inflow DOC concentration
        Inflow.DIC = NaN; %(scaled) inflow DIC concentration
        Inflow.O2 = NaN; %(scaled) inflow concentration
        Inflow.NO3 = NaN; %(scaled) inflow concentration
        Inflow.NH4 = NaN; %(scaled) inflow concentration
        Inflow.SO4 = NaN; %(scaled) inflow  concentration
        Inflow.Fe2 = NaN; %(scaled) inflow concentration
        Inflow.Ca2 = NaN; %(scaled) inflow concentration
        Inflow.pH = NaN; %(scaled) inflow concentration
        Inflow.CH4aq = NaN; %(scaled) inflow concentration
        Inflow.Fe3 = NaN; %(scaled) inflow concentration
        Inflow.Al3 = NaN; %(scaled) inflow concentration
        Inflow.FeS = NaN; %(scaled) inflow concentration
        Inflow.CaCO3 = NaN; %(scaled) inflow concentration
        Inflow.CH4g = NaN; %(scaled) inflow concentration
    end
end
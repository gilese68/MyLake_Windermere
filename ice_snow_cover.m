classdef ice_snow_cover
    %ICE_SNOW_COVER stores properties and methods related to ice and snow
    %cover
    %   A few properties are constant and defined in this class. Other
    %   properties are either user defined or calculated during the
    %   simulation. Methods updates properties during the simulation, e.g.
    %   ice and frazil thicknesses.
    properties (Constant)
        rho_fw = 1000;        %density of freshwater (kg m-3)
        rho_ice = 910;        %ice (incl. snow ice) density (kg m-3)
        rho_new_snow = 250;   %new-snow density (kg m-3)
        max_rho_snow = 450;   %maximum snow density (kg m-3)
        L_ice = 333500;       %latent heat of freezing (J kg-1)
        K_ice = 2.1;          %ice heat conduction coefficient (W m-1 K-1)
        C1 = 7.0;             %snow compaction coefficient #1
        C2 = 21.0;            %snow compaction coefficient #2
        Frazil2Ice_tresh = 0.03;  % treshold (m) where frazil is assumed to turn into a solid ice cover NEW!!!
    end
    properties
        tank_id         % tank number, for information
        rho_snow
        Tice
        XE_melt
        XE_surf
        Hi
        WEQs
        Hsi
        HFrazil
        IceIndicator
        IceSnowAttenuationCoefficient
        pp
        qq
        DoF
        DoM
        lambda_i        %PAR light attenuation coefficient for ice (m-1)
        lambda_s        %PAR light attenuation coefficient for snow (m-1)
    end
    
    methods
        function obj = ice_snow_cover(Ice0, t, lambda_i, lambda_s)
            %ICE_SNOW_COVER Construct an instance of this class
            %   Initialise properties and set user defined properties.
            % ice & snow parameter values            
            obj.rho_snow = obj.rho_new_snow;   %initial snow density (kg m-3)
            obj.Tice = NaN;                %ice surface temperature (initial value, deg C)
            obj.XE_melt = 0;               %energy flux that is left from last ice melting (initial value, W m-2)
            obj.XE_surf = 0;               %energy flux from water to ice (initial value,  J m-2 per day)

            %Initialisation of ice & snow variables
            obj.Hi = Ice0(1);               %total ice thickness (initial value, m)
            obj.WEQs = (obj.rho_snow / obj.rho_fw) * Ice0(2); %snow water equivalent  (initial value, m)
            obj.Hsi = 0;                %snow ice thickness (initial value = 0 m)
            obj.HFrazil = 0;              % (initial value, m) NEW!!!
            
            if ((obj.Hi <= 0) && (obj.WEQs > 0))
                error('Mismatch in initial ice and snow thicknesses')
            end

            if (obj.Hi<=0)
                obj.IceIndicator = 0;     %IceIndicator==0 means no ice cover
            else
                obj.IceIndicator = 1;
            end
            
            obj.pp = 1; %initial indexes for ice freezing/melting date arrays
            obj.qq = 1;
            obj.DoF = []; %initialize freezing dates
            obj.DoM = []; %initialize melting dates
            
            obj.tank_id = t; %tank number
            obj.lambda_i = lambda_i;
            obj.lambda_s = lambda_s;
        end
        
        function obj = set_attenuation_coefficient(obj)
            %SET_ATTENUATION_COEFFICIENT calculates the attenuation
            %coefficient depending on the ice and snow cover.
            if(obj.IceIndicator == 0)
                obj.IceSnowAttenuationCoefficient = 1; %no extra light attenuation due to snow and ice
            else    %extra light attenuation due to ice and snow
                obj.IceSnowAttenuationCoefficient = exp(-obj.lambda_i * obj.Hi) * exp(-obj.lambda_s * (obj.rho_fw / obj.rho_snow) * obj.WEQs);
            end
        end
        
        function [obj, dT] = frazil_ice_melting(obj, dT, Az, Cw, Vz)
            %FRAZING_ICE_MELTING calculates the frazil thickness
            % === NEW!!! frazil ice melting ===
            postemp = find(dT>0);
            if (isempty(postemp) == 0)
                %disp(['NOTE: positive night heat flux at T=' num2str(Tz(postemp),2)]) %NEW
                RelT = dT(postemp) ./ sum(dT(postemp));
                HFrazilnew = max(0, obj.HFrazil - sum(dT(postemp)) * 1 / ((Az(1) * obj.rho_ice * obj.L_ice) / (Cw * Vz(1)))); %
                sumdTnew = max(0, sum(dT(postemp))-(obj.HFrazil * Az(1) * obj.rho_ice * obj.L_ice) / (Cw * Vz(1)));
                dT(postemp) = RelT .* sumdTnew;
                obj.HFrazil = HFrazilnew;
            end
        end
        
        function [obj, Hi_new, dWEQs, dWEQnews] = ice_and_snow_formation(obj, i, Tf, Wt)
            %ICE_AND_SNOW_FORMATION updates the snow and ice thicknesses
            %with snow precipitation
            %Calculate ice surface temperature (Tice)
            if(obj.WEQs == 0) %if no snow
                alfa = 1 / (10 * obj.Hi);
                dHsi = 0;
            else
                K_snow = 2.22362 * (obj.rho_snow / 1000)^1.885; %Yen (1981)
                alfa = (obj.K_ice / K_snow) * (((obj.rho_fw / obj.rho_snow) * obj.WEQs) / obj.Hi);
                %Slush/snow ice formation (directly to ice)
                dHsi = max([0, obj.Hi * (obj.rho_ice / obj.rho_fw - 1) + obj.WEQs]);
                obj.Hsi = obj.Hsi + dHsi;
            end
            obj.Tice = (alfa * Tf + Wt(i,3)) / (1 + alfa);

            %Ice growth by Stefan's law
            Hi_new = sqrt((obj.Hi + dHsi)^2 + (2 * obj.K_ice / (obj.rho_ice * obj.L_ice)) * (24*60*60) * (Tf-obj.Tice));
            %snow fall
            dWEQnews = 0.001 * Wt(i,7); %mm->m
            dWEQs = dWEQnews - dHsi *(obj.rho_ice / obj.rho_fw); % new precipitation minus snow-to-snowice in snow water equivalent
        end
        
        function [obj, Hi_new, dWEQs, dWEQnews] = ice_and_snow_melting(obj, Tf, Qsw, Qlw, Qsl)
            %ICE_AND_SNOW_MELTING updates ice and snow thicknesses with
            %melting at the lake surface
            obj.Tice = Tf;    %ice surface at freezing point
            dWEQnews = 0; %No new snow
            if (obj.WEQs > 0)
                %snow melting in water equivalents
                dWEQs = -max([0, (60*60*24)*(((1 - obj.IceSnowAttenuationCoefficient) * Qsw) + Qlw + Qsl) / (obj.rho_fw * obj.L_ice)]);
                if ((obj.WEQs + dWEQs) < 0) %if more than all snow melts...
                    Hi_new = obj.Hi + (obj.WEQs + dWEQs) * (obj.rho_fw / obj.rho_ice); %...take the excess melting from ice thickness
                else
                    Hi_new = obj.Hi; %ice does not melt until snow is melted away
                end
            else
                %total ice melting
                dWEQs=0;
                Hi_new = obj.Hi - max([0, (60*60*24)*(((1 - obj.IceSnowAttenuationCoefficient) * Qsw) + Qlw + Qsl) / (obj.rho_ice * obj.L_ice)]);
                %snow ice part melting
                obj.Hsi = obj.Hsi - max([0, (60*60*24)*(((1 - obj.IceSnowAttenuationCoefficient) * Qsw) + Qlw + Qsl) / (obj.rho_ice * obj.L_ice)]);
                if (obj.Hsi <= 0)
                    obj.Hsi = 0;
                end
            end %if there is snow or not
        end
        
        function [obj] = update_ice_and_snow_thicknesses(obj, Hi_new, dWEQs)
            %UPDATE_ICE_AND_SNOW_THICKNESSES with melting due to heat flux
            %from water
            obj.Hi = Hi_new - (obj.XE_surf / (obj.rho_ice * obj.L_ice)); %new ice thickness (minus melting due to heat flux from water)
            obj.XE_surf = 0; %reset energy flux from water to ice (J m-2 per day)
            obj.WEQs = obj.WEQs + dWEQs; %new snow water equivalent

            if(obj.Hi < obj.Hsi)
                obj.Hsi = max(0, obj.Hi);    %to ensure that snow ice thickness does not exceed ice thickness
                %(if e.g. much ice melting much from bottom)
            end
        end
           
        function [obj] = update_snow_density(obj, dWEQnews, snow_compaction_switch, Wt, Tf, i)
            %UPDATE_SNOW_DENSITY with new snow density, and snow compaction
            %if switch is on.
            if(obj.WEQs <= 0)
                obj.WEQs = 0; %excess melt energy already transferred to ice above
                obj.rho_snow = obj.rho_new_snow;
            else
                %Update snow density as weighed average of old and new snow densities
                obj.rho_snow = obj.rho_snow * (obj.WEQs - dWEQnews) / obj.WEQs + obj.rho_new_snow * dWEQnews / obj.WEQs;
                if (snow_compaction_switch == 1)
                    %snow compaction
                    if (Wt(i,3) < Tf) %if air temperature is below freezing
                        rhos = 1e-3 * obj.rho_snow; %from kg/m3 to g/cm3
                        delta_rhos = 24 * rhos * obj.C1 * (0.5 * obj.WEQs) * exp(-obj.C2 * rhos) * exp(-0.08*(Tf - 0.5 * (obj.Tice + Wt(i,3))));
                        obj.rho_snow = min([obj.rho_snow + 1e+3 * delta_rhos, obj.max_rho_snow]);  %from g/cm3 back to kg/m3
                    else
                        obj.rho_snow = obj.max_rho_snow;
                    end
                end
            end
        end
        
        function [obj, Tz] = initial_freezing(obj, Tz, Tf, Vz, Cw, Az, i, M_start)
            %INITIAL FREEZING finds supercooled layers and initialises
            %freezing
            Supercooled = find(Tz < Tf);
            if (isempty(Supercooled)==0)
                %===NEW!!! (040707)
                if(Supercooled(1) ~= 1); disp('NOTE: non-surface subsurface supercooling'); end
                InitIceEnergy = sum((Tf - Tz(Supercooled)) .* Vz(Supercooled) * Cw);
                obj.HFrazil = obj.HFrazil + (InitIceEnergy / (obj.rho_ice * obj.L_ice)) / Az(1);
                Tz(Supercooled) = Tf;

                if ((obj.IceIndicator == 0) & (obj.HFrazil > obj.Frazil2Ice_tresh))
                    % NOTE to disable ice put zero hero; disable disp below
                    obj.IceIndicator = 1;
                    obj.Hi = obj.Hi + obj.HFrazil;
                    obj.HFrazil = 0;
                    obj.DoF(obj.qq) = i;
                    disp(['Tank ' num2str(obj.tank_id) ', Ice-on, ' datestr(datenum(M_start)+i-1)])
                    obj.qq = obj.qq + 1;
                end

                if (obj.IceIndicator == 1)
                    obj.Hi = obj.Hi + obj.HFrazil;
                    obj.HFrazil = 0;
                end
                Tz(1) = Tf; %set temperature of the first layer to freezing point
                %======================
            end
        end
        
        function obj = reset_ice_snow(obj, t, i, M_start)
            %RESET_ICE_SNOW ice has melted and variables are reset
            if(obj.Hi <= 0)
                obj.IceIndicator = 0;
                disp(['Tank ' num2str(t) ', Ice-off, ' datestr(datenum(M_start)+i-1)])
                obj.XE_melt = (-obj.Hi - (obj.WEQs * obj.rho_fw / obj.rho_ice)) * obj.rho_ice * obj.L_ice/(24*60*60);
                %(W m-2) snow part is in case ice has melted from bottom leaving some snow on top (reducing XE_melt)
                obj.Hi = 0;
                obj.WEQs = 0;
                obj.Tice = NaN;
                obj.DoM(obj.pp) = i;
                obj.pp = obj.pp+1;
            end
        end
    end
end


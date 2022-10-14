classdef sediments
    %SEDIMENTS regroups methods and properties related to sediments
    % Constant properties are stored in this class. Other user defined
    % properties are initialised at the beginning of the simulation.
    properties (Constant)
        K_sed = 0.035;      %thermal diffusivity of the sediment (m2 day-1)
        rho_sed = 2500;     %bulk density of the inorganic solids in sediment (kg m-3)
        rho_org = 1000;     %bulk density of the organic solids in sediment (kg m-3)
        cp_sed = 1000;      %specific heat capasity of the sediment (J kg-1 K-1) 
        N_sed = 26;         %total number of layers in the sediment column
    end
   
    properties
        F_sed_sld           %volume fraction of solids in sediment (= 1-porosity)
        Fmax_L_sed
        H_sed               %height of active sediment layer (m, wet mass)
        ksw                 %sediment pore water mass transfer coefficient (m/d)
    end
    
    methods
        function obj = sediments(Fmax_L, H_sed, F_sed_sld)
            %SEDIMENTS Construct an instance of this class
            %   Initialised user defined properties
            obj.Fmax_L_sed = Fmax_L;
            obj.H_sed = H_sed;
            obj.F_sed_sld = F_sed_sld;
            obj.ksw = 1e-3;
        end
        
        function [Qz_sed, Tzy_sed] = calculate_sediment_heat_flux(obj, sediment_heatflux_switch, Tzy_sed, Tz, dt, Nz, Az)
            %CALCULATE_SEDIMENT_HEAT_FLUX
            
            % Sediment vertical heat flux, Q_sed
            % (averaged over the whole top area of the layer, although actually coming only from the "sides")
            if (sediment_heatflux_switch == 1)
                % update top sediment temperatures
                dz_sf = 0.2; %fixed distance between the two topmost sediment layers (m)
                Tzy_sed(1,:) = Tz';
                Tzy_sed_upd = sedimentheat_v11(Tzy_sed, obj.K_sed, dt);
                Tzy_sed = Tzy_sed_upd;
                Qz_sed = obj.K_sed * obj.rho_sed * obj.cp_sed * (1/dz_sf) * (-diff([Az; 0]) ./ Az) .* (Tzy_sed(2,:)' - Tzy_sed(1,:)'); %(J day-1 m-2)
                %positive heat flux => from sediment to water
            else
                Qz_sed = zeros(Nz, 1);
            end
        end
    end
end


classdef switches
    %SWITCHES Store all switches value
    
    properties (Constant)
        snow_compaction_switch = 1;       %snow compaction: 0=no, 1=yes
        deposition_switch = 0;			%human impact, atm deposition , point source addition %% NEW_DOCOMO
        sediment_heatflux_switch = 1;     %heatflux from sediments: 0=no, 1=yes
        selfshading_switch = 1;           %light attenuation by chlorophyll a: 0=no, 1=yes
        tracer_switch = 1;                %simulate tracers:  0=no, 1=yes
        sediment_module = 0;            % sediment diagenetic module  %% NEW_DOCOMO
        wc_chemistry_module = 1;        % WC chemistry module: on/off
        wc_int_method = 0;              % WC chemistry module: method: 0 = Runge-Kutta 4th order; 1 = Buthcer's 5th Order;
        %fokema
        photobleaching = 0;               %photo bleaching: 0=TSA model, 1=FOKEMA model
        floculation_switch = 1;           % floculation according to Wachenfeldt 2008  %% NEW_DOCOMO
        resuspension_enabled = 0;         % Resuspension switch

        rate_estimator_switch = 0;        % save rates or not, additional cost of about 20% of computational time;
    end
    
end


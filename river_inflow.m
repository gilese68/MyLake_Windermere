function [tank] = river_inflow(Inflw, I_sc, Tf, i, ...
    ies80, zz, dz, tank, prev_profiles)

% Inflw key to columns 
% 1 = InflowQ	2 = InflowT	3 = POC	4 = InflowTP	5 = InflowDOP	6 = DOC	
% 7 = DIC	8 = O	9 = NO3 10 =	NH4	11 = SO4	12 = Fe2	13 = Ca2	
% 14 = pH	15 = CH4aq	16 =Fe3	17 = Al3 18 = FeS	19 = CaCO3	20 = CH4g
% 21 = POP	22 = Si


    tank.Inflow.V = I_sc.V * Inflw(i,1); % (scaled) inflow rate
    tank.Inflow.T = I_sc.T + Inflw(i,2); %(adjusted) inflow temperature
    if (tank.Inflow.T < Tf) %negative temperatures changed to Tf
        tank.Inflow.T = Tf;
    end
    tank.Inflow.POC = I_sc.POC * Inflw(i,3); %(scaled) inflow POC concentration
    tank.Inflow.TP = I_sc.TP * Inflw(i,4); %(scaled) inflow P concentration 
    tank.Inflow.DOP = I_sc.DOP * Inflw(i,5); %(scaled) inflow DOP concentration
    tank.Inflow.DOC = I_sc.DOC * Inflw(i,6); %(scaled) inflow DOC concentration
    tank.Inflow.DIC = I_sc.DIC * Inflw(i,7); %(scaled) inflow DIC concentration  
    tank.Inflow.O2 = I_sc.O * Inflw(i,8); %(scaled) inflow O2 concentration
    % inflow HS and H2S are neglected
    tank.Inflow.NO3 = I_sc.NO3 * Inflw(i,9);
    tank.Inflow.NH4 = I_sc.NH4 * Inflw(i,10);
    tank.Inflow.SO4 = I_sc.SO4 * Inflw(i,11);
    tank.Inflow.Fe2 = I_sc.Fe2 * Inflw(i,12);
    tank.Inflow.Ca2 = I_sc.Ca2 * Inflw(i,13);
    tank.Inflow.pH = I_sc.pH * Inflw(i,14);
    tank.Inflow.CH4aq = I_sc.CH4aq * Inflw(i,15);
    tank.Inflow.Fe3 = I_sc.Fe3 * Inflw(i,16);
    tank.Inflow.Al3 = I_sc.Al3 * Inflw(i,17);
    tank.Inflow.FeS = I_sc.FeS * Inflw(i,18);
    tank.Inflow.CaCO3 = I_sc.CaCO3 * Inflw(i,19);
    tank.Inflow.CH4g = I_sc.CH4g * Inflw(i,20);
    tank.Inflow.POP = I_sc.POP * Inflw(i,21);
    % TP add Si
    tank.Inflow.Si = I_sc.Si * Inflw(i,22); % %(scaled) inflow Si concentration
    % TP include loop round multiple algal types
    for ii = 1:size(tank.profiles.Chlz,2)
        tank.Inflow.Chl{ii} = I_sc.Chl(ii) * Inflw(i,22+ii); %(scaled) inflow Chl a concentration
    end
    
    tank.Inflow.Pz = (tank.Inflow.TP - tank.Inflow.DOP - tank.Inflow.POP)/2;
    tank.Inflow.Pz = tank.Inflow.Pz .* (tank.Inflow.Pz > 0);
    tank.Inflow.PP = (tank.Inflow.TP - tank.Inflow.DOP - tank.Inflow.POP)/2; %;
    tank.Inflow.PP = tank.Inflow.PP .* (tank.Inflow.PP > 0);
    
    if(tank.Inflow.V > 0)
        if (isnan(tank.Inflow.T))
            tank.lvlD = 0;
            tank.Inflow.T = prev_profiles.Tz(1, t);
        else
            tank.rho = polyval(ies80, max(0, prev_profiles.Tz)) + min(prev_profiles.Tz, 0);	% Density (kg/m3)
            rho_Inflow = polyval(ies80, max(0, tank.Inflow.T)) + min(tank.Inflow.T, 0);
            lvlG = find(tank.rho >= rho_Inflow);
            if (isempty(lvlG))
                lvlG = length(tank.rho);
            end
            tank.lvlD = zz(lvlG(1)); %level zz above which inflow is put
        end %if isnan...


        %Changes in properties due to inflow
        tank.profiles.Tz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Tz, tank.lvlD, tank.Inflow.V, tank.Inflow.T); %Temperature
        tank.profiles.POCz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.POCz, tank.lvlD, tank.Inflow.V, tank.Inflow.POC); % POC
        tank.profiles.DOPz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.DOPz, tank.lvlD, tank.Inflow.V, tank.Inflow.DOP); %Particulate organic P
             
        % TP include loop round multiple algal types
        for ii = 1:size(tank.profiles.Chlz,2)
          tank.profiles.Chlz{ii} = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Chlz{ii}, tank.lvlD, tank.Inflow.V, tank.Inflow.Chl(ii)); %Chlorophyll (group 1)
        end 
        
        tank.profiles.DOCz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.DOCz, tank.lvlD, tank.Inflow.V, tank.Inflow.DOC); %DOC
        tank.profiles.DICz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.DICz, tank.lvlD, tank.Inflow.V, tank.Inflow.DIC); %DIC
        tank.profiles.O2z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.O2z, tank.lvlD, tank.Inflow.V, tank.Inflow.O2); %O2
        tank.profiles.NO3z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.NO3z, tank.lvlD, tank.Inflow.V, tank.Inflow.NO3); %NO3
        tank.profiles.NH4z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.NH4z, tank.lvlD, tank.Inflow.V, tank.Inflow.NH4); %NH4
        tank.profiles.SO4z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.SO4z, tank.lvlD, tank.Inflow.V, tank.Inflow.SO4); %SO4
        tank.profiles.Fe2z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Fe2z, tank.lvlD, tank.Inflow.V, tank.Inflow.Fe2); %Fe2
        tank.profiles.Ca2z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Ca2z, tank.lvlD, tank.Inflow.V, tank.Inflow.Ca2); %Ca2
        tank.profiles.pHz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.pHz, tank.lvlD, tank.Inflow.V, tank.Inflow.pH); %pH
        tank.profiles.CH4aqz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.CH4aqz, tank.lvlD, tank.Inflow.V, tank.Inflow.CH4aq); %CH4aq
        tank.profiles.Fe3z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Fe3z, tank.lvlD, tank.Inflow.V, tank.Inflow.Fe3); %Fe3
        tank.profiles.Al3z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Al3z, tank.lvlD, tank.Inflow.V, tank.Inflow.Al3); %Al3
        tank.profiles.FeSz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.FeSz, tank.lvlD, tank.Inflow.V, tank.Inflow.FeS); %FeS
        tank.profiles.CaCO3z = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.CaCO3z, tank.lvlD, tank.Inflow.V, tank.Inflow.CaCO3); %CaCO3
        tank.profiles.CH4gz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.CH4gz, tank.lvlD, tank.Inflow.V, tank.Inflow.CH4g); %CH4g
        tank.profiles.POPz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.POPz, tank.lvlD, tank.Inflow.V, tank.Inflow.POP); %POP
        tank.profiles.PPz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.PPz, tank.lvlD, tank.Inflow.V, tank.Inflow.PP); %PP
        tank.profiles.Pz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Pz, tank.lvlD, tank.Inflow.V, tank.Inflow.Pz); %Pz
        tank.profiles.Siz = IOflow_v11(dz, zz, tank.Vz_tank, prev_profiles.Siz, tank.lvlD, tank.Inflow.V, tank.Inflow.Si); %Siz
    else
        tank.lvlD = NaN;
    end %if(Inflow.V>0)

end
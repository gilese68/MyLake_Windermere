function [profiles] = mixing_from_wind_energy(C_shelter, Az, Vz, tau, rho, dt, dz, zz, g, profiles, ies80)

    TKE = C_shelter * Az(1) * sqrt(tau^3/rho(1)) * (24*60*60*dt); %Turbulent kinetic energy (J day-1) over the whole lake

    %Wind mixing
    WmixIndicator = 1;
    Bef_wind = sum(diff(rho) == 0); %just a watch variable
    
    while (WmixIndicator == 1)
        d_rho = diff(rho);
        inx = find(d_rho>0);
        if (isempty(inx) == 0) %if water column not already fully mixed
            zb = inx(1);
            MLD = dz * zb; %mixed layer depth
            dD = d_rho(zb); %density difference
            Zg = sum( Az(1:zb+1) .* zz(1:zb+1) ) / sum(Az(1:zb+1)); %Depth of center of mass of mixed layer
            V_weight = Vz(zb+1) * sum(Vz(1:zb)) / (Vz(zb+1) + sum(Vz(1:zb)));
            POE = (dD * g * V_weight * (MLD + dz/2 - Zg));
            KP_ratio = TKE / POE;
            
            fn = {'Tz', 'POCz', 'Pz', 'Chlz', 'PPz', 'DOPz', 'DOCz', ...
                'DICz', 'O2z', 'NO3z', 'NH4z', 'SO4z', 'HSz', 'H2Sz', 'Fe2z', ...
                'Ca2z', 'pHz', 'CH4aqz', 'Fe3z', 'Al3z', 'FeSz', 'CaCO3z', ...
                'CH4gz', 'POPz','Siz'};
            if (KP_ratio >= 1)

                for k = 1:numel(fn)
                    if strcmp(fn{k}, 'Chlz') == 0
                        mix = sum(Vz(1:zb+1) .* profiles.(fn{k})(1:zb+1)) / sum(Vz(1:zb+1));
                        profiles.(fn{k})(1:zb+1) = mix;
                    else
                        for l=1:length(profiles.Chlz)
                            mix = sum(Vz(1:zb+1) .* profiles.Chlz{l}(1:zb+1)) / sum(Vz(1:zb+1));
                            profiles.Chlz{l}(1:zb+1) = mix;
                        end
                    end
                end

                rho = polyval(ies80, max(0, profiles.Tz)) + min(profiles.Tz,0);
                TKE = TKE - POE;

            else %if KP_ratio < 1, then mix with the remaining TKE part of the underlying layer

                for k = 1:numel(fn)
                    if strcmp(fn{k}, 'Chlz') == 0
                        mix = sum([Vz(1:zb); KP_ratio * Vz(zb+1)] .* profiles.(fn{k})(1:zb+1)) / sum([Vz(1:zb); KP_ratio * Vz(zb+1)]);
                        profiles.(fn{k})(1:zb) = mix;
                        profiles.(fn{k})(zb+1) = KP_ratio * mix + (1-KP_ratio) * profiles.(fn{k})(zb+1);
                    else
                        for l=1:length(profiles.Chlz)
                            mix = sum([Vz(1:zb); KP_ratio * Vz(zb+1)] .* profiles.Chlz{l}(1:zb+1)) / sum([Vz(1:zb); KP_ratio * Vz(zb+1)]);
                            profiles.Chlz{l}(1:zb) = mix;
                            profiles.Chlz{l}(zb+1) = KP_ratio * mix + (1-KP_ratio) * profiles.Chlz{l}(zb+1);
                        end
                    end
                end

                rho = polyval(ies80, max(0, profiles.Tz)) + min(profiles.Tz, 0);
                TKE = 0;
                WmixIndicator = 0;
            end %if (KP_ratio>=1)
        else
            WmixIndicator = 0;
        end %if water column (not) already mixed
    end %while

    Aft_wind = sum(diff(rho)==0); %just a watch variable
% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2004

% VERSION 1.2.1a, based on convection_v12 (with three modified lines of code, marked with NEW!!!)

% Convection module
% Code checked by TSA, 07.03.05
% Last modified by TSA, 17.07.07
% Modified by PK 30.12.2010 (DIC) & 14.02.2011 (O2)

function [profiles] = ...
    convection_v2(profiles, Tprof_prev, Vz, Cw, f_par, lambdaz_wtot_avg, zz, swa_b0, tracer_switch, springautumn)

% Inputs (with extension "_in") and Outputs:
%       Tz   : Temperature profile
%       Cz   : Tracer profile
%       Sz   : Suspended inorg. matter profile
%       Pz   : Dissolved inorg. P profile
%       Chlz : Chlorophyll a profile
%       PPz  : Phosphorus bound to inorganic particles profile
%       DOPz  : Dissolved organic phosphorus profile
%       DOCz  : Particulate inorganic phosphorus profile
%       DICz  : Dissolved inorganic carbon profile (PK)

% Inputs:
%       Tprof_prev   : Temperature profile from previous timestep
%       etc.

% These variables are still global and not transferred by functions
global ies80;

Trhomax = 3.98; %temperature of maximum water density (deg C)
Nz = length(zz); %total number of layers in the water column
dz = zz(2) - zz(1); %model grid step

% Convective mixing adjustment
% Mix successive layers until stable density profile
rho = polyval(ies80, max(0, profiles.Tz)) + min(profiles.Tz, 0);	% Density (kg/m3)
d_rho = [diff(rho); 1]; %d_rho = how much a layer is lighter than layer below; last cell in "d_rho" is always positive (sediment)

while any(d_rho < 0)
    blnUnstb_layers = (d_rho <= 0); %1=layer is heavier or equal than layer below, 0=layer is lighter than layer below
    A_Unstb = find(diff([0; blnUnstb_layers]) == 1); %layer index(es) where unstable/neutral column(s) start(s)
    B_Unstb = find(diff([0; blnUnstb_layers]) == -1) - 1;%layer index(es) where unstable/neutral column(s) end(s)

    for n = 1:length(A_Unstb)
        j = [A_Unstb(n):B_Unstb(n)+1];
        
        % TODO: [LC] move outside of convection
        fn = {'Tz', 'POCz', 'Pz', 'PPz', 'DOPz', 'DOCz', ...
            'DICz', 'O2z', 'NO3z', 'NH4z', 'SO4z', 'HSz', 'H2Sz', 'Fe2z', ...
            'Ca2z', 'pHz', 'CH4aqz', 'Fe3z', 'Al3z', 'FeSz', 'CaCO3z', ...
            'CH4gz', 'POPz','Siz'};
        for k = 1:numel(fn)
            mix = sum(profiles.(fn{k})(j) .* Vz(j)) / sum(Vz(j));
            profiles.(fn{k})(j) = mix * ones(size(profiles.(fn{k})(j)));
        end
        
        % Algae profiles
        for k = 1:length(profiles.Chlz)
            mix = sum(profiles.Chlz{k}(j) .* Vz(j)) / sum(Vz(j));
            profiles.Chlz{k}(j) = mix * ones(size(profiles.Chlz{k}(j)));
        end

    end

    rho = polyval(ies80,max(0,profiles.Tz)) + min(profiles.Tz,0);
    d_rho = [diff(rho); 1];
end

if (springautumn==1)
    % Spring/autumn turnover
    % don't allow temperature jumps over temperature of maximum density
    if(((Tprof_prev(1) > Trhomax) & (profiles.Tz(1) < Trhomax)) | ((Tprof_prev(1) < Trhomax) & (profiles.Tz(1) > Trhomax))) %NEW!!! (Tprof, ">" instead of ">=")!

        jumpinx = find(((Tprof_prev > Trhomax) & (profiles.Tz < Trhomax)) | ((Tprof_prev < Trhomax) & (profiles.Tz > Trhomax))); %NEW!!!!!
        if (sum(jumpinx == 1) == 0); disp('NOTE: Non-surface jumps over temperature of maximum density'); end %NEW!!!!!

        intSign = sign(Trhomax - profiles.Tz(1)); %plus in autumn turnover, minus in spring turnover
        XE_turn = cumsum((profiles.Tz - Trhomax) .* Vz * Cw * intSign); %always starts negative
        Dummy = find(XE_turn > 0);
        
        if(isempty(Dummy) == 1)
            profiles.Tz(:) = Trhomax;
            if (intSign == 1)
                profiles.Tz(1) = profiles.Tz(1) + intSign * XE_turn(end) / (Vz(1) * Cw); %put overshoot on top layer
            else
                profiles.Tz = profiles.Tz + (-diff(intSign * XE_turn(end) * (f_par * exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]) + ...
                    (1-f_par) * exp([0; -swa_b0*ones(Nz,1)] .* [zz; zz(end)+dz])) ) ./(Vz*Cw));
                %distribute overshoot as shortwave energy
            end
        else
            profiles.Tz(1:Dummy(1)-1) = Trhomax;
            profiles.Tz(Dummy(1)) = Trhomax + intSign * XE_turn(Dummy(1)) / (Vz(Dummy(1)) * Cw);
        end
    end
end %springautumn

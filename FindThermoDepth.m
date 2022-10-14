function [ thermoD,thermoInd,drho_dz,SthermoD,SthermoInd ] = ...
    FindThermoDepth( rhoVar,depths,Smin )
%----Author: Jordan S Read 2009 ----
% updated 27 march 2010

warning('off', 'signal:findpeaks:noPeaks');
warning('off', 'signal:findpeaks:largeMinPeakHeight');

seasonal = boolean(0);
if nargout > 3
    seasonal = 1;
end
if nargin < 3
    Smin = 0.1;
end
dRhoPerc = 0.15; %min percentage max for unique thermocline step
numDepths = length(depths);
drho_dz = NaN(1,numDepths-1);

for i = 1:numDepths-1
    drho_dz(i) = (rhoVar(i+1)-rhoVar(i))/...
        (depths(i+1)-depths(i));
end
if seasonal
    %look for two distinct maximum slopes, lower one assumed to be seasonal
    [mDrhoZ,thermoInd] = max(drho_dz);          %find max slope
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (drho_dz(thermoInd+1)-drho_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (drho_dz(thermoInd)-drho_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
    end
    dRhoCut = max([dRhoPerc*mDrhoZ Smin]);
    [pks,locs] = findpeaks(drho_dz,'minpeakheight',dRhoCut);
    if isempty(pks)
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    else
        mDrhoZ = pks(length(pks));
        SthermoInd = locs(length(pks));
        if SthermoInd > thermoInd+1
            SthermoD = mean([depths(SthermoInd)...
                depths(SthermoInd+1)]);
            if SthermoInd > 1 && SthermoInd < numDepths-1
                Sdn = -(depths(SthermoInd+1)-depths(SthermoInd))/...
                    (drho_dz(SthermoInd+1)-drho_dz(SthermoInd));
                Sup = (depths(SthermoInd)-depths(SthermoInd-1))/...
                    (drho_dz(SthermoInd)-drho_dz(SthermoInd-1));
                upD  = depths(SthermoInd);
                dnD  = depths(SthermoInd+1);
                SthermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));
            end
        else
            SthermoD = thermoD;
            SthermoInd = thermoInd;
        end
    end
    if SthermoD < thermoD;
        SthermoD = thermoD;
        SthermoInd = thermoInd;
    end
    
else
    [mDrhoZ,thermoInd] = max(drho_dz);          %find max slope
        thermoD = mean([depths(thermoInd)...
        depths(thermoInd+1)]);                  %depth of max slope
    if thermoInd > 1 && thermoInd < numDepths-1 %if within range, 
        Sdn = -(depths(thermoInd+1)-depths(thermoInd))/...
            (drho_dz(thermoInd+1)-drho_dz(thermoInd));
        Sup = (depths(thermoInd)-depths(thermoInd-1))/...
            (drho_dz(thermoInd)-drho_dz(thermoInd-1));
        upD  = depths(thermoInd);
        dnD  = depths(thermoInd+1);
        thermoD = dnD*(Sdn/(Sdn+Sup))+upD*(Sup/(Sdn+Sup));% Interpolate mixed depth
    end
end
warning('on', 'signal:findpeaks:noPeaks');
warning('on','signal:findpeaks:largeMinPeakHeight');

end

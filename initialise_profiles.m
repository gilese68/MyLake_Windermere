function profiles = initialise_profiles(In, zz, dz, N_alg)
    %INITIALISE_PROFILES set all profiles to user defined initial profiles
    
    fn = {'Tz', 'POCz', 'TPz', 'DOPz', 'DOCz', ...
    'DICz', 'O2z', 'TPz_sed', 'Chlz_sed', 'FIM', 'NO3z', 'NH4z', 'SO4z', 'HSz', 'H2Sz', 'Fe2z', ...
    'Ca2z', 'pHz', 'CH4aqz', 'Fe3z', 'Al3z', 'FeSz', 'CaCO3z', ...
    'CH4gz', 'POPz', 'Siz'};
    for k = 1:numel(fn)
        profiles.(fn{k}) = interp1(In.Z, In.(fn{k}), zz+dz/2);
    end

    profiles.Pz = (profiles.TPz-profiles.DOPz-profiles.POPz);
    profiles.PPz = (profiles.TPz-profiles.DOPz-profiles.POPz); % (mg m-3) NEW!!!
    profiles.Pz = profiles.Pz .* (profiles.Pz>0);
    profiles.PPz = profiles.PPz .* (profiles.PPz>0);
    
    % Algae profiles
    for k = 1:N_alg
        profiles.Chlz{k} = interp1(In.Z, In.Chlz(:,k), zz+dz/2);
    end
    
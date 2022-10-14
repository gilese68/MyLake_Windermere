classdef vertical_diffusion
    %VERTICAL_DIFFUSION Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        % diffusion parameterization exponents
        Kz_b1 = 0.43
        Kz_b1_ice = 0.43
    end
    
    properties
        Kz_K1           % open water diffusion parameter (-)
        Kz_K1_ice       % under ice diffusion parameter (-)
        Kz_N0           % min. stability frequency (s-2)
    end
    
    methods
        function obj = vertical_diffusion(Kz_K1, Kz_K1_ice, Kz_N0)
            %VERTICAL_DIFFUSION Construct an instance of this class
            %   Detailed explanation goes here
            obj.Kz_K1 = Kz_K1;
            obj.Kz_K1_ice = Kz_K1_ice;
            obj.Kz_N0 = Kz_N0;           
        end
        
        function [Fi, Kz] = calculate_turbulent_vertical_diffusion(obj, Tz, zz, ice_snow, Vz, Az, dz, dt, ies80, g)
            % Vertical turbulent diffusion
            rho = polyval(ies80, max(0, Tz(:))) + min(Tz(:), 0);  % Water density (kg m-3)
            % Note: in equations of rho it is assumed that every supercooled degree lowers density by
            % 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")

            N2 = g * (diff(log(rho)) ./ diff(zz));	% Brunt-Vaisala frequency (s-2) for level (zz+1) %NOTE: why log here??
            if (ice_snow.IceIndicator == 0)
                Kz = obj.Kz_K1 * max(obj.Kz_N0, N2).^(-obj.Kz_b1);	% Vertical diffusion coeff. in ice free season (m2 day-1)
                % for level (zz+1)
            else
                Kz = obj.Kz_K1_ice * max(obj.Kz_N0, N2).^(-obj.Kz_b1_ice); % Vertical diffusion coeff. under ice cover (m2 day-1)
                % for level (zz+1)
            end

            Fi = tridiag_DIF_v11([NaN; Kz], Vz, Az, dz, dt); %Tridiagonal matrix for general diffusion
        end
    end
end

% Below are the two functions for calculating tridiagonal matrix Fi for solving the
% 1) diffusion equation (tridiag_DIF_v11), and
% 2) advection-diffusion equation (tridiag_HAD_v11) by fully implicit hybrid exponential numerical scheme,
% based on Dhamotharan et al. 1981,
%'Unsteady one-dimensional settling of suspended sediments', Water Resources Research 17(4), 1125-1132
% code checked by TSA, 16.03.2004


%Inputs:
% Kz    diffusion coefficient at layer interfaces (plus surface) N (N,1)
% U     vertical settling velocity (scalar)
% Vz    layer volumes (N,1)
% Az    layer interface areas (N,1)
% dz    grid size
% dt    time step

%Output:
% Fi    tridiagonal matrix for solving new profile Cz

% az = (dt/dz) * [0; Kz] .* (Az ./ Vz);
% bz = (dt/dz) * [Kz; 0] .* ([Az(2:end); 0] ./ Vz);
% Gi = [-bz (1 + az + bz) -az];

%=== DIFFUSIVE EQUATION ===
function [Fi] = tridiag_DIF_v11(Kz,Vz,Az,dz,dt)

Nz=length(Vz); %number of grid points/layers

% Linearized heat conservation equation matrix (diffusion only)
az = (dt/dz) * Kz .* (Az ./ Vz);                                        %coefficient for i-1
cz = (dt/dz) * [Kz(2:end); NaN] .* ([Az(2:end); NaN] ./ Vz);            %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i+1
%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged
bz(1)= 1 + az(1) + cz(1);


%Boundary conditions, bottom

%az(end) remains unchanged
cz(end) = 0;
bz(end) = 1 + az(end) + cz(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function
end

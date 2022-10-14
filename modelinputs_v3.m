% === MyLake model, Version 3
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
% modified by JBA consulting and Lancaster University 2020
% Module for reading input data and parameters

function [In, tt, Ice0, Wt, Inflw, Phys_par, Bio_par, Tanks_par] ...
    = modelinputs_v3(M_start,M_stop,init_filename,...
    input_filename,param_filename,dt,tanks_params)

% Inputs:
%       M_start : Model start date [year, month, day]
%       M_stop : Model stop date [year, month, day]
%       + Input filenames and sheetnames
%		dt		: Time step (= 1 day)
% Outputs:
%		tt		: Solution time domain (day)
%       In_Z    : Depths read from initial profiles file (m)
%       In_Az   : Areas read from initial profiles file (m2)
%       In_Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In_Cz   : Initial tracer profile read from initial profiles file (-)
%       In_POCz   : Initial POC (mg m-3)
%       In_TPz  : Initial total P profile read from initial profiles file (mg m-3)
%       In_DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In_Chlz : Initial chlorophyll a profile read from initial profiles file (mg m-3)
%       In_DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In_TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_Chlz_sed : Initial chlorophyll a profile in the sediment compartments read from initial profiles file (mg m-3)
% ...
%       In_FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 15 parameters  that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (15 * 2)
%       Bio_par_names   : Names for Bio_par

global ies80;

% == Read model parameter file

data = load(param_filename);
% read 46 lines
par = data(:, 2:end);
Phys_par = par(1:37,:);
Bio_par = par(38:61,:);

% Vertical settling velocities
U = Bio_par(8:9,:); %for sedimenting velocities
if any(U<0)
    error('Given settling velocity must be positive')
end

%% Tank parameters
Tanks_par.Nt = tanks_params{1,1}; % Number of tanks
if tanks_params{2,1}(end-3:end)=='.mat'
    load(tanks_params{2,1});
    Tanks_par.tank_advection_matrix = tank_advection_matrix;
else
    Tanks_par.tank_advection_matrix = xlsread(tanks_params{2,1});
end
if tanks_params{3,1}(end-3:end)=='.mat'
    load(tanks_params{3,1});
    Tanks_par.tank_diffusion_matrix = tank_diffusion_matrix;
else
    Tanks_par.tank_diffusion_matrix = xlsread(tanks_params{3,1});
end
if ~isempty(tanks_params{4,1})
    Tanks_par.tank_areas = xlsread(tanks_params{4,1});
else
    Tanks_par.tank_areas = NaN;
end
Tanks_par.river_inflow_switch = tanks_params{5,1};
Tanks_par.river_inflow_tank_id = tanks_params{6,1};
Tanks_par.river_in_coef = tanks_params{7,1};
Tanks_par.river_out_coef = tanks_params{8,1};
Tanks_par.decline_parameter_advection = tanks_params{9,1};
Tanks_par.decline_parameter_diffusion = tanks_params{10,1};
Tanks_par.init_filename_list = tanks_params{11,1};


%% == Read morphometric and INITIAL DEPTH PROFILES

% Use init_filename if list from Tanks_par is empty
if isempty(Tanks_par.init_filename_list)
    init_filename_list = init_filename;
else
    init_filename_list = Tanks_par.init_filename_list;
end

if length(init_filename_list) ~= Tanks_par.Nt && length(init_filename_list) ~= 1
    error('Number of profile files must be equal to 1 or the number of tanks.')
end

In = struct;
Ice0 = zeros(length(init_filename_list), 2);
for i=1:length(init_filename_list)
    init_filename = init_filename_list(i);
    InitMx = dlmread(init_filename{1}(~isspace(init_filename{1})), '\t', 2, 0);
    %% changed below 3:end to 1:end
    In(i).Z = InitMx(1:end,1);
    In(i).Az = InitMx(1:end,2);
    In(i).Tz = InitMx(1:end,3);
    In(i).POCz = InitMx(1:end,4);
    In(i).TPz = InitMx(1:end,5);
    In(i).DOPz = InitMx(1:end,6);
    In(i).DOCz = InitMx(1:end,7);
    In(i).TPz_sed = InitMx(1:end,8);
    In(i).Chlz_sed = InitMx(1:end,9);
    In(i).FIM = InitMx(1:end,10);
    Ice0(i, :) = InitMx(1,11:12);
    In(i).O2z = InitMx(1:end,13);
    In(i).DICz = InitMx(1:end,14);
    In(i).NO3z = InitMx(1:end,15);
    In(i).NH4z = InitMx(1:end,16);
    In(i).SO4z = InitMx(1:end,17);
    In(i).HSz = InitMx(1:end,18);
    In(i).H2Sz = InitMx(1:end,19);
    In(i).Fe2z = InitMx(1:end,20);
    In(i).Ca2z = InitMx(1:end,21);
    In(i).pHz = InitMx(1:end,22);
    In(i).CH4aqz = InitMx(1:end,23);
    In(i).Fe3z = InitMx(1:end,24);
    In(i).Al3z = InitMx(1:end,25);
    In(i).FeSz = InitMx(1:end,26);
    In(i).CaCO3z = InitMx(1:end,27);
    In(i).CH4gz = InitMx(1:end,28);
    In(i).POPz = InitMx(1:end,29);
    In(i).Siz = InitMx(1:end,30); % New silica profile LC: 2020
    % Now from 1: to N species appended to array
    In(i).Chlz = InitMx(1:end,31:end);
    
end

tt = [datenum(M_start):dt:datenum(M_stop)]';		% Solution time domain

%%%%%%%%%%%%%%%%% Read INFLOW TIMESERIES DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
InputMx = dlmread(input_filename, ',', 2, 0);

%% changed below 3:end to 1:end 2011-09-29
In_Date=InputMx(1:end,1:3);
In_Met=InputMx(1:end,4:10);
In_Inflow=InputMx(1:end,11:end);

tmet=datenum(In_Date);

dum=100*((tmet(end)-tmet(1)+1)-length(tmet))/(tmet(end)-tmet(1)+1);
% $$$ disp(['Percent missing dates in meteorology and inflow data: ']);
% $$$ disp([num2str(dum) ' %']);

dum=100*sum(isnan(In_Met))./length(tmet);
% $$$ disp(['Percent missing values in meteorology data (values correspond to columns 4-10 in input file): ']);
% $$$ disp([num2str(dum) ' %']);

dum=100*sum(isnan(In_Inflow))./length(tmet);
% $$$ disp(['Percent missing values in inflow data (values correspond to columns 11-17 in input file): ']);
% $$$ disp([num2str(dum) ' %']);
% $$$ disp(' ')

clear Wt
for i=1:7 %Interpolate over missing values and dates
    nonnans = find(isnan(In_Met(:,i))==0);
    if(isempty(nonnans)) % if the whole column is NaNs then preserve it
        Wt(:,i) = NaN*ones(length(tt(:)),1);
    else
        repaired = interp1(nonnans,In_Met(nonnans,i),[1:length(In_Met(:,i))]);
        %  sometimes tmet non-unique. to check it:
        % u=unique(tmet);
        %         n=histc(tmet,u);
        % u(n>1)
        % datestr(u(n>1))
        % ans =
        % 01-Mar-2100
        
        Wt(:,i) = interp1(tmet, repaired, tt(:));
    end
end
% Wt(:,1)  Global radiation (MJ/(m^2 day))
% Wt(:,2)  Cloud cover (-)
% Wt(:,3)  Air temperature (deg. C, at 2 m height)
% Wt(:,4)  Relative humidity (%, at 2 m height)
% Wt(:,5)  Air pressure (mbar)
% Wt(:,6)  Wind speed (m/s at 10 m height)
% Wt(:,7)  Precipitation (mm/day)

clear Inflw
for i=1:size(In_Inflow,2) %Interpolate over missing values and dates
    nonnans = find(isnan(In_Inflow(:,i))==0);
    if(isempty(nonnans)) % if the whole column is NaNs then preserve it
        Inflw(:,i) = NaN*ones(length(tt(:)),1);
    else
        repaired = interp1(nonnans,In_Inflow(nonnans,i),[1:length(In_Inflow(:,i))]);
        Inflw(:,i) = interp1(tmet, repaired, tt(:));
    end
end

% International Equation of State 1980
% 5-order polynomial for density as function of temperature
ies80 = [6.536332e-9,-1.120083e-6,1.001685e-4,-9.09529e-3,6.793952e-2,999.842594];

% Default turbulence and wind shelter parameterization (Hondzo and Stefan, 1993; Ellis et al., 1991)
if(isnan(Phys_par(2)))
    Phys_par(2) = 0.00706*(In_Az(1)/1e6)^0.56; % default diffusion coeff. parameterisation
end

if(isnan(Phys_par(3)))
    Phys_par(3) = 8.98e-4;		%default value for diffusion coeff. in ice-covered water
end

if(isnan(Phys_par(4)))
    Phys_par(4) = 7e-5;			% default minimum allowed stability frequency, N2 > N0 <=> Kz < Kmax (1/s2)
end

if(isnan(Phys_par(5)))
    Phys_par(5) =  1-exp(-0.3*In_Az(1)/1e6);			% default wind sheltering parameterisation
end

function [ sediment_concentrations, sediment_params, sediment_matrix_templates] = sediment_init( pH, max_depth, temperature, rate_estimator_switch)
  % Input:
  % bottom pH of the lake, temperature at SWI.

  % Output:
  % init sediment_concentrations - concentrations read from file
  % sediment_params - all sediment params
  % sediment_matrix_templates - all matrix templates for sediment module
  % species_sediment - which species to simulate
    global sed_par_file sediment_params


    sediment_params = params(max_depth, temperature);
    pH = 8;
    sediment_concentrations = init_concentrations(pH); % NOTE: init pH in sediment = 8;
    sediment_matrix_templates = templates();
    sediment_params.rate_estimator_switch = rate_estimator_switch;
end

function [sediment_params] = params(max_depth, temperature)
    global sed_par_file
    f=fopen(sed_par_file);
    % f=fopen('calibration_k_values.txt');
    data = textscan(f,'%s%f', 75,'Delimiter', '\t');
    fclose(f);

    % chemical constants from file
    sediment_params.k_Chl = data{2}(1);
    sediment_params.k_POP = data{2}(2);
    sediment_params.k_POC = data{2}(3);
    sediment_params.k_DOP = data{2}(4);
    sediment_params.k_DOC = data{2}(5);
    sediment_params.Km_O2 = data{2}(6);
    sediment_params.Km_NO3 = data{2}(7);
    sediment_params.Km_FeOH3 = data{2}(8);
    sediment_params.Km_FeOOH = data{2}(9);
    sediment_params.Km_SO4 = data{2}(10);
    sediment_params.Km_oxao = data{2}(11);
    sediment_params.Km_amao = data{2}(12);
    sediment_params.Kin_O2 = data{2}(13);
    sediment_params.Kin_NO3 = data{2}(14);
    sediment_params.Kin_FeOH3 = data{2}(15);
    sediment_params.Kin_FeOOH = data{2}(16);
    sediment_params.k_amox = data{2}(17);
    sediment_params.k_Feox = data{2}(18);
    sediment_params.k_Sdis = data{2}(19);
    sediment_params.k_Spre = data{2}(20);
    sediment_params.k_fes2pre = data{2}(21);
    sediment_params.k_alum = data{2}(22);
    sediment_params.k_pdesorb_a = data{2}(23);
    sediment_params.k_pdesorb_b = data{2}(24);
    sediment_params.k_fesox = data{2}(25);
    sediment_params.k_fes2ox = data{2}(26);
    sediment_params.k_tS_Fe = data{2}(27);
    sediment_params.K_FeS = data{2}(28);
    sediment_params.k_fe_dis = data{2}(29);
    sediment_params.k_fe_pre = data{2}(30);
    sediment_params.k_apa_pre = data{2}(31);
    sediment_params.k_apa_dis = data{2}(32);
    sediment_params.K_apa = data{2}(33);
    sediment_params.k_CCpre = data{2}(34);
    sediment_params.k_CCdis = data{2}(35);
    sediment_params.K_CC = data{2}(36);
    sediment_params.k_FCpre = data{2}(37);
    sediment_params.k_FCdis = data{2}(38);
    sediment_params.K_FC = data{2}(39);
    sediment_params.k_viv_pre = data{2}(40);
    sediment_params.k_viv_dis = data{2}(41);
    sediment_params.K_viv = data{2}(42);
    sediment_params.k_oms = data{2}(43);
    sediment_params.k_tsox = data{2}(44);
    sediment_params.k_fespre = data{2}(45);
    sediment_params.k_ch4_o2 = data{2}(46);
    sediment_params.k_ch4_so4 = data{2}(47);
    sediment_params.Kh_CH4 = data{2}(48);
    sediment_params.k_ch4_dis = data{2}(49);
    sediment_params.w_CH4g = data{2}(50);
    sediment_params.Kh_CO2 = data{2}(51);
    sediment_params.accel = data{2}(52);
    sediment_params.Kd_fe2 = data{2}(53);
    sediment_params.k_pdesorb_c = data{2}(54);
    phi_in = data{2}(55);
    phi_f  = data{2}(56);
    X_b   = data{2}(57);
    tortuosity = data{2}(58);
    sediment_params.w = data{2}(59);      % time-dependent burial rate w = 0.1
    sediment_params.n = data{2}(60);;  % points in spatial grid
    sediment_params.depth = data{2}(61);  % sediment depth
    alpha0 = data{2}(62);
    sediment_params.Cx1 = data{2}(63);
    sediment_params.Ny1 = data{2}(64);
    sediment_params.Pz1 = data{2}(65);
    sediment_params.Cx2 = data{2}(66);
    sediment_params.Ny2 = data{2}(67);
    sediment_params.Pz2 = data{2}(68);
    sediment_params.Cx3 = data{2}(69);
    sediment_params.Ny3 = data{2}(70);
    sediment_params.Pz3 = data{2}(71);
    sediment_params.effective_depth = data{2}(72);
    sediment_params.n_of_time_steps_during_1_dt_of_myLake = data{2}(73);  % time step
    sediment_params.pH_algorithm = data{2}(74);
    sediment_params.SO4_bottom_fx = data{2}(75);


    % Estimation of params:
    temperature = 4.8497; % mean value of observed bottom temp; NOTE: We need to make better estimation of temp for diff coef.
    abs_temp = temperature+273.15; % [K]
    P_atm = 1.01; % [Bar]
    P_water = 998 * 9.8 * max_depth/10^5; % [Bar]
    pressure = P_atm + P_water; % [Bar]
    salinity = 0;
    viscosity = viscosity_fun(temperature,pressure,salinity);


    % pH module. NOTE: estimates pH once a day
    % !!!!!!! Recommend to use #0
    % Specify pH algorithm:
    % 0. Disabled (simulation time: ~80 sec per year)
    % 1. Phreeqc  adds computational time. (simulation time: ~ 120 sec per year) - (under test)
    % 2. New algorithm by Markelov (simulation time: ~ 180 sec per year) (under test)



    aq_system.carb_acid = acid([6.52, 10.56], 0, 0);
    aq_system.amonia = acid([9.2503], 1, 0);
    aq_system.sulf = acid([6.8861], 0, 0);
    aq_system.ca = neutral(2, 0);
    aq_system.fe2 = neutral(2, 0);
    aq_system.no3 = neutral(-1, 0);
    aq_system.so4 = neutral(-2, 0);
    aq_system.p_acid= acid([2.148, 7.198, 12.319], 0, 0);
    sediment_params.aq_system = aq_system;


    sediment_params.years = 1/365;  % 1 day #35
    sediment_params.ts = 1/365/(sediment_params.n_of_time_steps_during_1_dt_of_myLake);  % time step
    x = linspace(0, sediment_params.depth, sediment_params.n);
    sediment_params.x = x;      % x-axis
    sediment_params.phi =  ones(size(x'))*0.92 ;%( phi_in - phi_f ) * exp( -x' / X_b ) + phi_f;;      % porosity
    sediment_params.tortuosity = tortuosity;  % tortuosity
    alfax = alpha0*exp(-0.25*x);
    sediment_params.alfax = alfax';   % bioirrigation
    sediment_params.viscosity = viscosity;
    sediment_params.temperature = temperature;
    sediment_params.pressure = pressure;
    sediment_params.salinity = salinity;

    %===================== % effective molecular diffusion
    % Linear regression of Diffusion coefficients for cations and anions (Boudreau, 1997):
    sediment_params.D_H    = lr_ion_diffusion(54.4, 1.555, temperature);
    sediment_params.D_OH   = lr_ion_diffusion(25.9, 1.094, temperature);
    sediment_params.D_HCO3 = lr_ion_diffusion(5.06, 0.275, temperature);
    sediment_params.D_CO3  = lr_ion_diffusion(4.33, 0.199, temperature);
    sediment_params.D_NO3  = lr_ion_diffusion(9.5,  0.388, temperature);
    sediment_params.D_SO4  = lr_ion_diffusion(4.88, 0.232, temperature);
    sediment_params.D_NH4  = lr_ion_diffusion(9.5,  0.413, temperature);
    sediment_params.D_Fe2  = lr_ion_diffusion(3.31, 0.15,  temperature);
    sediment_params.D_PO4  = lr_ion_diffusion(2.62, 0.143, temperature);
    sediment_params.D_Ca2  = lr_ion_diffusion(3.6,  0.179, temperature);
    sediment_params.D_HS   = lr_ion_diffusion(10.4, 0.273, temperature);

    % Empirical correlation of Wilke and Chang (1955) as corrected by Hayduk and Laudie (1974)
    sediment_params.D_NH3 = hayduk_laudie_diffusion(viscosity, abs_temp, 24.5);
    sediment_params.D_O2  = hayduk_laudie_diffusion(viscosity, abs_temp, 27.9);
    sediment_params.D_CO2 = hayduk_laudie_diffusion(viscosity, abs_temp, 37.3);
    sediment_params.D_CH4aq = hayduk_laudie_diffusion(viscosity, abs_temp, 37.7);
    sediment_params.D_CH4g = hayduk_laudie_diffusion(viscosity, abs_temp, 37.7);
    sediment_params.D_CO2g = hayduk_laudie_diffusion(viscosity, abs_temp, 37.7);

    % Diffusion coefficient based on Einstein relation:


    % User specified diffusion coefficients and other params (if there is no values found above):
    sediment_params.D_H2S = 284;
    sediment_params.D_HS  = 284;
    sediment_params.D_S0  = 100;
    sediment_params.D_DOP  = 85.14; %  0.27 ?? 10-5 cm2 s-1 taken from Diffusion processes of soluble organic substances in soil and their effect on ecological processes, Roland
    sediment_params.D_DOC  = 85.14; %  0.27 ?? 10-5 cm2 s-1 taken from Diffusion processes of soluble organic substances in soil and their effect on ecological processes, Roland
    sediment_params.Db    = 5;

    sediment_params.CH4_solubility = sediment_params.Kh_CH4 * pressure * 1e3; % in mmol/L
    sediment_params.CO2_solubility = sediment_params.Kh_CO2 * pressure * 1e3; % in mmol/L


end

function [sediment_concentrations ] = init_concentrations(pH)
    %% Init concentrations of sediment species
    global sediment_params

    z_max = sediment_params.depth;
    n = sediment_params.n;
    dz  = z_max/(n-1);
    read_file = true;
    if read_file
        Init = dlmread('IO/sediment_initial_concentrations.txt', '\t', 1, 0);
        zz = [0:dz:z_max]'; %solution depth domain
        in_z = Init(1:end, 1);
        sediment_concentrations.POP = interp1(in_z, Init(1:end, 2), zz);
        sediment_concentrations.POC = interp1(in_z, Init(1:end, 3), zz);
        sediment_concentrations.DOP = interp1(in_z, Init(1:end, 4), zz);
        sediment_concentrations.DOC = interp1(in_z, Init(1:end, 5), zz);
        sediment_concentrations.O2 = interp1(in_z, Init(1:end, 6), zz);
        sediment_concentrations.NO3 = interp1(in_z, Init(1:end, 7), zz);
        sediment_concentrations.NH4 = interp1(in_z, Init(1:end, 8), zz);
        sediment_concentrations.NH3 = interp1(in_z, Init(1:end, 9), zz);
        sediment_concentrations.FeOH3 = interp1(in_z, Init(1:end, 10), zz);
        sediment_concentrations.FeOOH = interp1(in_z, Init(1:end, 11), zz);
        sediment_concentrations.Fe2 = interp1(in_z, Init(1:end, 12), zz);
        sediment_concentrations.SO4 = interp1(in_z, Init(1:end, 13), zz);
        sediment_concentrations.H2S = interp1(in_z, Init(1:end, 14), zz);
        sediment_concentrations.HS = interp1(in_z, Init(1:end, 15), zz);
        sediment_concentrations.PO4 = interp1(in_z, Init(1:end, 16), zz);
        sediment_concentrations.PO4adsa = interp1(in_z, Init(1:end, 17), zz);
        sediment_concentrations.PO4adsb = interp1(in_z, Init(1:end, 18), zz);
        sediment_concentrations.S0 = interp1(in_z, Init(1:end, 19), zz);
        sediment_concentrations.S8 = interp1(in_z, Init(1:end, 20), zz);
        sediment_concentrations.FeS = interp1(in_z, Init(1:end, 21), zz);
        sediment_concentrations.FeS2 = interp1(in_z, Init(1:end, 22), zz);
        sediment_concentrations.AlOH3 = interp1(in_z, Init(1:end, 23), zz);
        sediment_concentrations.Ca2 = interp1(in_z, Init(1:end, 24), zz);
        sediment_concentrations.Ca3PO42 = interp1(in_z, Init(1:end, 25), zz);
        sediment_concentrations.OMS = interp1(in_z, Init(1:end, 26), zz);
        sediment_concentrations.H3O   = interp1(in_z, Init(1:end, 27), zz);
        sediment_concentrations.CaCO3 = interp1(in_z, Init(1:end, 28), zz);
        sediment_concentrations.CO2 = interp1(in_z, Init(1:end, 29), zz);
        sediment_concentrations.CO3 = interp1(in_z, Init(1:end, 30), zz);
        sediment_concentrations.HCO3 = interp1(in_z, Init(1:end, 31), zz);
        sediment_concentrations.CO2g = interp1(in_z, Init(1:end, 32), zz);
        sediment_concentrations.Chl = interp1(in_z, Init(1:end, 33), zz);
        sediment_concentrations.CH4aq = interp1(in_z, Init(1:end, 34), zz);
        sediment_concentrations.CH4g = interp1(in_z, Init(1:end, 35), zz);
        sediment_concentrations.FeCO3 = interp1(in_z, Init(1:end, 36), zz);
        sediment_concentrations.Fe3PO42 = interp1(in_z, Init(1:end, 37), zz);
        sediment_concentrations.PO4adsc = interp1(in_z, Init(1:end, 38), zz);
    else
        sediment_concentrations.POP     = ones(n,1) * 0;
        sediment_concentrations.POC     = ones(n,1) * 0;
        sediment_concentrations.DOP    = ones(n,1) * 0;
        sediment_concentrations.DOC    = ones(n,1) * 0;
        sediment_concentrations.O2      = ones(n,1) * 0;
        sediment_concentrations.NO3     = ones(n,1) * 0;
        sediment_concentrations.NH4     = ones(n,1) * 0;
        sediment_concentrations.NH3     = ones(n,1) * 0;
        sediment_concentrations.FeOH3   = ones(n,1) * 0;
        sediment_concentrations.FeOOH   = ones(n,1) * 0;
        sediment_concentrations.Fe2     = ones(n,1) * 0;
        sediment_concentrations.SO4     = ones(n,1) * 0;
        sediment_concentrations.H2S     = ones(n,1) * 0;
        sediment_concentrations.HS      = ones(n,1) * 0;
        sediment_concentrations.PO4     = ones(n,1) * 0;
        sediment_concentrations.PO4adsa = ones(n,1) * 0;
        sediment_concentrations.PO4adsb = ones(n,1) * 0;
        sediment_concentrations.S0      = ones(n,1) * 0;
        sediment_concentrations.S8      = ones(n,1) * 0;
        sediment_concentrations.FeS     = ones(n,1) * 0;
        sediment_concentrations.FeS2    = ones(n,1) * 0;
        sediment_concentrations.AlOH3   = ones(n,1) * 0;
        sediment_concentrations.Ca2     = ones(n,1) * 0;
        sediment_concentrations.Ca3PO42 = ones(n,1) * 0;
        sediment_concentrations.OMS     = ones(n,1) * 0;
        sediment_concentrations.H3O    = ones(n,1) * 10^-pH*1e3;
        sediment_concentrations.CaCO3   = ones(n,1) * 0;
        sediment_concentrations.CO2     = ones(n,1) * 0;
        sediment_concentrations.CO3     = ones(n,1) * 0;
        sediment_concentrations.HCO3    = ones(n,1) * 0;
        sediment_concentrations.CO2g   = ones(n,1) * 0;
        sediment_concentrations.Chl   = ones(n,1) * 0;
        sediment_concentrations.CH4aq   = ones(n,1) * 0;
        sediment_concentrations.CH4g   = ones(n,1) * 0;
        sediment_concentrations.FeCO3   = ones(n,1) * 0;
        sediment_concentrations.Fe3PO42   = ones(n,1) * 0;
        sediment_concentrations.PO4adsc   = ones(n,1) * 0;

    end
end

function [D] = einstein_diffusion(D_ref, abs_temp, viscosity)
  %% einstein_diffusion: Einstein's relation of diffusion coefficient derived from Einstein's formula D/D0
  % input:
  T_ref = 273.15 + 25; % [K]
  viscosity_ref = 1.002*10^-2; % [g s-1 cm-1]
  % D_ref - reference diffusion at 25'C [cm^2 year-1]
  % viscosity_ref - reference viscosity at 25'C
  % viscosity - current viscosity
  % abs_temp - current temperature in K

  % Output:
  % D - diffusion coefficient [cm^2 year-1]

  D = D_ref*(abs_temp/T_ref)*(viscosity_ref/viscosity);
end

function [D] = hayduk_laudie_diffusion(viscosity, abs_temp, molar_volume)
  %% hayduk_laudie_diffusion: empirical correlation of Wilke and Chang (1955) as corrected by Hayduk and Laudie (1974):

  % Input:
  % viscosity
  % abs_temp - temperature in K
  % molar_volume - is the molar volume of the nonelectrolyte (at the normal boiling temperature of that solute). [experimental]
  % 3.156*10^7 [seconds] ---> [year]

  % Output:
  % Diffusion coeff. in [cm^2 year-1]
  D = 4.71978 * 10^-9 * abs_temp / ( viscosity * molar_volume^0.6 ) * 3.156 * 10^7;
end

function [D] = lr_ion_diffusion(m0, m1, t)
  %%  lr_ion_diffusion: Linear Regressions??? of the Infinite- Dilution Diffusion Coefficients Do for Anions and Cations against Temperature (Boudreau, B.P., 1997. Diagenetic Models and Their Implementation)

  % Output:
  % Diffusion coeff. in [cm^2 year-1]
  % 3.156*10^7 [seconds] ---> [year]

  D = ( m0 + m1*t ) * 10^-6 * 3.156 * 10^7;
end

function [u] = viscosity_fun(temperature,pressure,salinity)
  %% viscosity: Values of the dynamic viscosity can be calculated from an empirical equation
  % developed by Matthaus (as quoted in Kukulka et al., 1987), which is claimed to be accurate to within 0.7% for the temperature, t, salinity, S, and pressure, P, ranges of 0
  % ??? C ??? 30, 0 ??? S ??? 36, and 1 to 1000 bars, respectively (Boudreau, 1997):

  % Input units
  t = temperature; % [Celsius]
  p = pressure; % [Bar]
  s = salinity; % [ppt] NOTE:not sure about units here!

  % Output:
  % viscosity [g cm-1 s-1] = Poise units
  u = (1.7910 - (6.144*10^-2)*t + (1.4510*10^-3)*t^2 - (1.6826*10^-5)*t^3 - (1.5290*10^-4)*p + (8.3885*10^-8)*p^2 + (2.4727*10^-3)*s + t*(6.0574*10^-6 *p - 2.6760*10^-9*p^2) + s*(4.8429*10^-5*t - 4.7172*10^-6*t^2 + 7.5986*10^-8*t^3))*10^-2;

end

function [sediment_matrix_templates] = templates()
    global sediment_params
    n  = sediment_params.n;
    depth = sediment_params.depth;
    x = linspace(0,depth,n);
    dx = x(2)-x(1);
    t = 0:sediment_params.ts:sediment_params.years;
    dt    = t(2)-t(1);
    v = sediment_params.w;
    phi  = sediment_params.phi;
    tortuosity = sediment_params.tortuosity;
    Db    = sediment_params.Db;

    % formation of templates:
    % Solid template the same for all solid species due to diffusion and advection coef the same for all.
    [Solid_AL, Solid_AR, solid_flux_coef] = cn_template_neumann(Db, v, phi, dx, dt, n);

    % solute templates:
    [O2_AL, O2_AR]        = cn_template_dirichlet(sediment_params.D_O2 + Db, tortuosity, v, phi, dx, dt, n);
    [NO3_AL, NO3_AR]        = cn_template_dirichlet(sediment_params.D_NO3 + Db, tortuosity, v, phi, dx, dt, n);
    [NH4_AL, NH4_AR]        = cn_template_dirichlet(sediment_params.D_NH4 + Db, tortuosity, v, phi, dx, dt, n);
    [Fe2_AL, Fe2_AR]        = cn_template_dirichlet(sediment_params.D_Fe2 + Db, tortuosity, v, phi, dx, dt, n);
    [H2S_AL, H2S_AR]        = cn_template_dirichlet(sediment_params.D_H2S + Db, tortuosity, v, phi, dx, dt, n);
    [S0_AL, S0_AR]        = cn_template_dirichlet(sediment_params.D_S0 + Db, tortuosity, v, phi, dx, dt, n);
    [PO4_AL, PO4_AR]        = cn_template_dirichlet(sediment_params.D_PO4 + Db, tortuosity, v, phi, dx, dt, n);
    [Ca2_AL, Ca2_AR]        = cn_template_dirichlet(sediment_params.D_Ca2 + Db, tortuosity, v, phi, dx, dt, n);
    [HS_AL, HS_AR]        = cn_template_dirichlet(sediment_params.D_HS + Db, tortuosity, v, phi, dx, dt, n);
    [OH_AL, OH_AR]        = cn_template_dirichlet(sediment_params.D_OH + Db, tortuosity, v, phi, dx, dt, n);
    [CO2_AL, CO2_AR]        = cn_template_dirichlet(sediment_params.D_CO2 + Db, tortuosity, v, phi, dx, dt, n);
    [CO3_AL, CO3_AR]        = cn_template_dirichlet(sediment_params.D_CO3 + Db, tortuosity, v, phi, dx, dt, n);
    [HCO3_AL, HCO3_AR]        = cn_template_dirichlet(sediment_params.D_HCO3 + Db, tortuosity, v, phi, dx, dt, n);
    [NH3_AL, NH3_AR]        = cn_template_dirichlet(sediment_params.D_NH3 + Db, tortuosity, v, phi, dx, dt, n);
    [CO2g_AL, CO2g_AR]        = cn_template_dirichlet(sediment_params.D_CO2g + Db, tortuosity, v, phi, dx, dt, n);
    [DOP_AL, DOP_AR]        = cn_template_dirichlet(sediment_params.D_DOP + Db, tortuosity, v, phi, dx, dt, n);
    [DOC_AL, DOC_AR]        = cn_template_dirichlet(sediment_params.D_DOC + Db, tortuosity, v, phi, dx, dt, n);
    [CH4aq_AL, CH4aq_AR]        = cn_template_dirichlet(sediment_params.D_CH4aq + Db, tortuosity, v, phi, dx, dt, n);
    [CH4g_AL, CH4g_AR]        = cn_template_dirichlet(sediment_params.D_CH4g + Db, tortuosity, sediment_params.w_CH4g, phi, dx, dt, n);  % NOTE: Rising velocity of the gas

    % custom boundary condittions for SO4. In Vansjo there is a source of SO4 at the bottom.
    [SO4_AL, SO4_AR, SO4_flux_coef]        = cn_template_custom(sediment_params.D_SO4 + Db, tortuosity, v, phi, dx, dt, n);

    sediment_params.solid_flux_coef = solid_flux_coef;
    sediment_params.SO4_flux_coef = SO4_flux_coef;

    sediment_matrix_templates = {...

        Solid_AL, Solid_AR, 'Solid';  % 1
        O2_AL, O2_AR,    'O2'; % 2
        NO3_AL, NO3_AR, 'NO3'; % 3
        SO4_AL, SO4_AR, 'SO4'; % 4
        NH4_AL, NH4_AR, 'NH4'; % 5
        Fe2_AL, Fe2_AR, 'Fe2'; % 6
        H2S_AL, H2S_AR, 'H2S'; % 7
        S0_AL, S0_AR, 'S0'; % 8
        PO4_AL, PO4_AR, 'PO4'; % 9
        Ca2_AL, Ca2_AR, 'Ca2'; % 10
        HS_AL, HS_AR, 'HS'; % 11
        CO2_AL, CO2_AR, 'CO2'; % 12
        CO3_AL, CO3_AR, 'CO3'; % 13
        HCO3_AL, HCO3_AR, 'HCO3'; % 14
        NH3_AL, NH3_AR, 'NH3'; % 15
        CO2g_AL, CO2g_AR, 'CO2g'; %16
        DOP_AL, DOP_AR, 'DOP'; %17
        DOC_AL, DOC_AR, 'DOC'; %18
        CH4aq_AL, CH4aq_AR, 'DOC'; %19
        CH4g_AL, CH4g_AR, 'DOC'; %20
    };

end


function [AL, AR] = cn_template_dirichlet(D_m, tortuosity, v, phi, dx, dt, n)
  %MATRICES Formation of matrices for species
  % ======================================================================

    D = D_m / tortuosity^2;
    s = phi * D * dt / dx / dx;
    q = phi * v * dt / dx;


    AL      = spdiags([-s/2-q/4 phi+s -s/2+q/4],[-1 0 1],n,n);
    AL(1,1) = phi(1);
    AL(1,2) = 0;
    AL(n,n) = phi(n)+s(n);
    AL(n,n-1) = -s(n);

    AR      = spdiags([ s/2+q/4 phi-s s/2-q/4],[-1 0 1],n,n);
    AR(1,1) = phi(1);
    AR(1,2) = 0;
    AR(n,n) = phi(n)-s(n);
    AR(n,n-1) = s(n);
end

function [AL, AR, flux_coef] = cn_template_neumann(D, v, phi, dx, dt, n)
  %MATRICES Formation of matrices for species
  % ======================================================================

    s = (1-phi) * D * dt / dx / dx;
    q = (1-phi) * v * dt / dx;

    AL      = spdiags([-s/2-q/4 (1-phi+s) -s/2+q/4],[-1 0 1],n,n);
    AL(1,1) = 1-phi(1)+s(1);
    AL(1,2) = -s(1);
    AL(n,n) = 1-phi(n)+s(n);
    AL(n,n-1) = -s(n);

    AR      = spdiags([ s/2+q/4 (1-phi-s)  s/2-q/4],[-1 0 1],n,n);
    AR(1,1) = 1-phi(1)-s(1);
    AR(1,2) = +s(1);
    AR(n,n) = 1-phi(n)-s(n);
    AR(n,n-1) = s(n);

    flux_coef = dx / D * (2*s(1) - q(1)) / (1-phi(1)); % NOTE: Tested, it is correct
end

function [AL, AR, flux_coef] = cn_template_custom(D_m, tortuosity, v, phi, dx, dt, n)
  %MATRICES Formation of matrices for species with custom BC
  % in this example it is Sulfate in Vansjo model
  % Dirichlet at the top and Neumann at the bottom
  % ======================================================================

    D = D_m / tortuosity^2;
    s = phi * D * dt / dx / dx;
    q = phi * v * dt / dx;

    AL      = spdiags([-s/2-q/4 phi+s -s/2+q/4],[-1 0 1],n,n);
    AL(1,1) = phi(1);
    AL(1,2) = 0;
    AL(n,n) = phi(n)+s(n);
    AL(n,n-1) = -s(n);

    AR      = spdiags([ s/2+q/4 phi-s s/2-q/4],[-1 0 1],n,n);
    AR(1,1) = phi(1);
    AR(1,2) = 0;
    AR(n,n) = phi(n)-s(n);
    AR(n,n-1) = s(n);

    flux_coef = dx / D * (2*s(end) - q(end)) / phi(end);
end

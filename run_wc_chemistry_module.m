function [profiles, C_new, wc_rates_av, mylake_temp_results] = run_wc_chemistry_module(profiles, DayFrac, lambdaz_wtot, H_sw_z, mylake_params, sediment_params, dt, wc_int_method, N_algae)

% Conversions
profiles.O2z     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.O2z, 31998.8);
profiles.DOCz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.DOCz, 12010.7);
profiles.NO3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.NO3z, 62004);
profiles.Fe3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Fe3z, 106867.0); %Fe(OH)3
profiles.SO4z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.SO4z, 96062);
profiles.NH4z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.NH4z, 18038);
profiles.Fe2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Fe2z, 55845);
profiles.H2Sz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.H2Sz, 34080.9);
profiles.HSz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.HSz, 33072.9);
profiles.Pz      = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Pz, 30973.762); % PO4
profiles.Al3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Al3z, 78003.6); % Al(OH)3
profiles.PPz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.PPz, 30973.762); % PO4-s
profiles.Ca2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Ca2z, 40078.0); % Ca2+
profiles.CO2z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CO2z, 44009.5);
profiles.DOPz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.DOPz, 30973.762); %DOPz
profiles.POCz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.POCz,  12010.7);
profiles.POPz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.POPz,  30973.762);
profiles.CH4aqz  = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CH4aqz,  16042.5);
profiles.CH4gz   = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CH4gz,  16042.5);
profiles.HCO3z   = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.HCO3z,  61016.8);
profiles.CO3z    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CO3z,  60008.9);
profiles.CaCO3z  = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CaCO3z,  100086.9); % Molar mass of CaCO3 = 100.0869 moles/gram
profiles.FeSz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.FeSz,  87910.0);  
profiles.Siz    = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Siz, 60840); % Molar mass of SiO2 = 60.084 moles/gram
for k=1:N_algae
    ASM(:,k) = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Chlz{k}, 30973.762); % TP only to pass to C0 conversion function
    profiles.Chlz{k} = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Chlz{k}, 30973.762); 
end

% Passing some MyLake results for chem module
mylake_temp_results.Tz = profiles.Tz;
mylake_temp_results.DayFrac = DayFrac;
mylake_temp_results.lambdaz_wtot = lambdaz_wtot;
mylake_temp_results.H_sw_z = H_sw_z;
mylake_temp_results.pHz = profiles.pHz;

%C0 = [profiles.O2z, profiles.DOCz, profiles.NO3z, profiles.Fe3z, profiles.SO4z, profiles.NH4z, profiles.Fe2z, profiles.H2Sz, profiles.HSz, profiles.Pz, profiles.Al3z, profiles.PPz, profiles.Ca2z, profiles.CO2z, profiles.DOPz, profiles.POCz, profiles.POPz, profiles.CH4aqz, profiles.CH4gz, profiles.HCO3z, profiles.CO3z, profiles.CaCO3z, profiles.FeSz, profiles.Siz, profiles.Chlz{1}, profiles.Chlz{2}];
C0 = [profiles.O2z, profiles.DOCz, profiles.NO3z, profiles.Fe3z, profiles.SO4z, profiles.NH4z, profiles.Fe2z, profiles.H2Sz, profiles.HSz, profiles.Pz, profiles.Al3z, profiles.PPz, profiles.Ca2z, profiles.CO2z, profiles.DOPz, profiles.POCz, profiles.POPz, profiles.CH4aqz, profiles.CH4gz, profiles.HCO3z, profiles.CO3z, profiles.CaCO3z, profiles.FeSz, profiles.Siz, ASM];
[C_new, wc_rates_av] = wc_chemical_reactions_module(mylake_params, sediment_params, mylake_temp_results, C0, dt, sediment_params.n_of_time_steps_during_1_dt_of_myLake, wc_int_method);

profiles.O2z  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,1), 31998.8);
profiles.DOCz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,2), 12010.7);
profiles.NO3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,3), 62004);
profiles.Fe3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,4), 106867.0); %Fe(OH)3
profiles.SO4z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,5), 96062);
profiles.NH4z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,6), 18038);
profiles.Fe2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,7), 55845);
profiles.H2Sz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,8), 34080.9);
profiles.HSz  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,9), 33072.9);
profiles.Pz   = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,10), 30973.762); % PO4
profiles.Al3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,11), 78003.6);
profiles.PPz  = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,12), 30973.762);
profiles.Ca2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,13), 40078.0); % Ca2+
profiles.CO2z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,14), 44009.5);
profiles.DOPz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,15), 30973.762); %DOPz
profiles.POCz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,16), 12010.7);
profiles.POPz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,17), 30973.762);
profiles.CH4aqz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,18), 16042.5);
profiles.CH4gz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,19), 16042.5);
profiles.HCO3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,20), 61016.8);
profiles.CO3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,21), 60008.9);
profiles.CaCO3z = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,22), 100086.9);
profiles.FeSz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,23), 87910.0);
profiles.Siz = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,24), 60840);
% TP algae
for k=1:N_algae
    profiles.Chlz{k}    = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C_new(:,24+k), 30973.762);
end

function C = convert_mg_per_qubic_m_to_umol_per_qubic_cm(C,M_C)
C = C./M_C;
%end of function

function C = convert_umol_per_qubic_cm_to_mg_per_qubic_m(C,M_C)
C = C.*M_C;
%end of function

function [C_new, rates] = wc_chemical_reactions_module(mylake_params, sediment_params, mylake_temp_results, C0, dt_mylake, n_of_time_steps_during_1_dt_of_myLake, method)
% ts - how many time steps during 1 day
% dt - time step in chemical module [years]
% dt_mylake - time step in MyLake [days]
% n_of_time_steps_during_1_dt_of_myLake - amount of steps during of 1 dt of MyLake;


dt = (dt_mylake/365) / n_of_time_steps_during_1_dt_of_myLake;

if method == 0
    [C_new, rates] = rk4(mylake_params, sediment_params, mylake_temp_results, C0, dt, n_of_time_steps_during_1_dt_of_myLake);
elseif method == 1
    [C_new, rates] = butcher5(mylake_params, sediment_params, mylake_temp_results, C0, dt, n_of_time_steps_during_1_dt_of_myLake);
end
C_new = (C_new>0).*C_new;

%end of function


%% rk4: Runge-Kutta 4th order integration
function [C_new, rates_av] = rk4(mylake_params, sediment_params, mylake_temp_results, C0, dt,n)
% ts - time step in chemical module [years]
% dt - time step in MyLake [days]
% (1/365/ts) = is how many steps during 1 dt of Mylake
% (dt/365) = is conversion of [days] to [years]

for i = 1:n
    [dcdt_1, r_1] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+0.5.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+0.5.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0+k_3, dt);
    k_4 = dt.*dcdt_4;
    C_new = C0 + (k_1+2.*k_2+2.*k_3+k_4)/6;
    C0 = C_new;
    
    % average rate
    if mylake_params.rate_estimator_switch
        fields = fieldnames(r_1);
        for fld_idx = 1:numel(fields)
            r.(fields{fld_idx}) = (r_1.(fields{fld_idx}) + 2*r_2.(fields{fld_idx}) + 2*r_3.(fields{fld_idx}) + r_4.(fields{fld_idx}))/6;
        end
        
        rates(i) = r;
    end
end

if mylake_params.rate_estimator_switch
    fields = fieldnames(rates);
    for i = 1:numel(fields)
        rates_av.(fields{i}) = 0;
        for j=1:n-1
            rates_av.(fields{i}) = rates_av.(fields{i}) + rates(j).(fields{i});
        end
        rates_av.(fields{i}) = rates_av.(fields{i})/(n-1);
    end
else
    rates_av = false;
end

%% butcher5: Butcher's Fifth-Order Runge-Kutta
function [C_new, rates_av] = butcher5(mylake_params, sediment_params, mylake_temp_results, C0, dt,n)

for i = 1:n
    [dcdt_1, r_1] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 1/4.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 1/8.*k_1 + 1/8.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 - 1/2.*k_2 + k_3, dt);
    k_4 = dt.*dcdt_4;
    [dcdt_5, r_5] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 + 3/16.*k_1 + 9/16.*k_4, dt);
    k_5 = dt.*dcdt_5;
    [dcdt_6, r_6] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C0 - 3/7.*k_1 + 2/7.*k_2 + 12/7.*k_3 - 12/7.*k_4 + 8/7.*k_5, dt);
    k_6 = dt.*dcdt_6;
    C_new = C0 + (7.*k_1 + 32.*k_3 + 12.*k_4 + 32.*k_5 + 7.*k_6)/90;
    C0 = C_new;
    
    if mylake_params.rate_estimator_switch
        fields = fieldnames(r_1);
        for fld_idx = 1:numel(fields)
            r.(fields{fld_idx}) = (7*r_1.(fields{fld_idx}) + 32*r_3.(fields{fld_idx}) + 12*r_4.(fields{fld_idx}) + 32*r_5.(fields{fld_idx}) + 7*r_6.(fields{fld_idx}))/90;
        end
        rates(i) = r;
    end
end

if mylake_params.rate_estimator_switch
    fields = fieldnames(rates);
    for i = 1:numel(fields)
        rates_av.(fields{i}) = 0;
        for j=1:ts-1
            rates_av.(fields{i}) = rates_av.(fields{i}) + rates(j).(fields{i});
        end
        rates_av.(fields{i}) = rates_av.(fields{i})/(ts-1);
    end
else
    rates_av = false;
end

function [dcdt, r] = wc_rates(mylake_params, sediment_params, mylake_temp_results, C, dt)
% parameters for water-column chemistry
% NOTE: the rates are the same as in sediments (per year!! not day). Units are per "year" due to time step is in year units too;

if any(isnan(C))
    error('NaN')
end

dcdt=zeros(size(C));

O2 = C(:,1) .* (C(:,1)>0);
DOC = C(:,2) .* (C(:,2)>0);
NO3 = C(:,3) .* (C(:,3)>0);
Fe3 = C(:,4) .* (C(:,4)>0);
SO4 = C(:,5) .* (C(:,5)>0);
NH4 = C(:,6) .* (C(:,6)>0);
Fe2 = C(:,7) .* (C(:,7)>0);
H2S = C(:,8) .* (C(:,8)>0);
HS = C(:,9) .* (C(:,9)>0);
Pz = C(:,10) .* (C(:,10)>0);
Al3 = C(:,11) .* (C(:,11)>0);
PP = C(:,12) .* (C(:,12)>0);
Ca2 = C(:,13) .* (C(:,13)>0);
CO2 = C(:,14) .* (C(:,14)>0);
DOP = C(:,15) .* (C(:,15)>0);
POC = C(:,16) .* (C(:,16)>0);
POP = C(:,17) .* (C(:,17)>0);
CH4aq = C(:,18) .* (C(:,18)>0);
CH4g = C(:,19) .* (C(:,19)>0);
HCO3 = C(:,20) .* (C(:,20)>0);
CO3 = C(:,21) .* (C(:,21)>0);
CaCO3 = C(:,22) .* (C(:,22)>0);
FeS = C(:,23) .* (C(:,23)>0);
Si = C(:,24).* (C(:,24)>0);
% TP generalise for multiple species

for ii = 1:size(mylake_params.m_twty,2)
    Chl{ii} = C(:,ii+24) .* (C(:,ii+24)>0);
end

k_Chl =  sediment_params.k_Chl;
k_POP =  sediment_params.k_POP;
k_POC = sediment_params.k_POC;
k_DOP = sediment_params.k_DOP;
k_DOC = sediment_params.k_DOC;
Km_O2 = sediment_params.Km_O2;
Km_NO3 = sediment_params.Km_NO3;
Km_FeOH3 = sediment_params.Km_FeOH3;
Km_FeOOH = sediment_params.Km_FeOOH;
Km_SO4 = sediment_params.Km_SO4;
Km_oxao = sediment_params.Km_oxao;
Km_amao = sediment_params.Km_amao;
Kin_O2 = sediment_params.Kin_O2;
Kin_NO3  = sediment_params.Kin_NO3;
Kin_FeOH3 = sediment_params.Kin_FeOH3;
Kin_FeOOH = sediment_params.Kin_FeOOH;
k_amox = sediment_params.k_amox;
k_Feox = sediment_params.k_Feox;
k_Sdis = sediment_params.k_Sdis;
k_Spre = sediment_params.k_Spre;
k_fes2pre = sediment_params.k_fes2pre;
k_pdesorb_c = sediment_params.k_pdesorb_c;
k_pdesorb_a = sediment_params.k_pdesorb_a;
k_pdesorb_b = sediment_params.k_pdesorb_b;
% k_alum = sediment_params.k_alum;
k_fesox   = sediment_params.k_fesox;
k_fes2ox   = sediment_params.k_fes2ox;
k_tS_Fe = sediment_params.k_tS_Fe;
K_FeS = sediment_params.K_FeS;
k_fe_dis = sediment_params.k_fe_dis;
k_fe_pre = sediment_params.k_fe_pre;
k_apa_pre  = sediment_params.k_apa_pre;
k_apa_dis  = sediment_params.k_apa_dis;
K_apa = sediment_params.K_apa;
k_CCpre = sediment_params.k_CCpre;
k_CCdis = sediment_params.k_CCdis;
K_CC = sediment_params.K_CC;
k_FCpre = sediment_params.k_FCpre;
k_FCdis = sediment_params.k_FCdis;
K_FC = sediment_params.K_FC;
k_viv_pre = sediment_params.k_viv_pre;
k_viv_dis = sediment_params.k_viv_dis;
K_viv = sediment_params.K_viv;
k_oms = sediment_params.k_oms;
k_tsox = sediment_params.k_tsox;
k_fespre = sediment_params.k_fespre;
k_ch4_o2 = sediment_params.k_ch4_o2;
k_ch4_so4 = sediment_params.k_ch4_so4;
Kh_CH4 = sediment_params.Kh_CH4;
Kh_CO2 = sediment_params.Kh_CO2;
k_ch4_dis = sediment_params.k_ch4_dis;
accel = sediment_params.accel;
Cx1   = sediment_params.Cx1;
Ny1   = sediment_params.Ny1;
Pz1   = sediment_params.Pz1;
Cx2   = sediment_params.Cx2;
Ny2   = sediment_params.Ny2;
Pz2   = sediment_params.Pz2;
Cx3   = sediment_params.Cx3;
Ny3   = sediment_params.Ny3;
Pz3   = sediment_params.Pz3;

z = mylake_params.zz;

T_ref = mylake_params.T_ref;
Q10 = mylake_params.Q10;
dop_twty = mylake_params.dop_twty;
% Algal parameters
g_twty = mylake_params.g_twty;
m_twty = mylake_params.m_twty;
P_half = mylake_params.P_half;
N_half = mylake_params.N_half;
Si_half = mylake_params.Si_half;

theta_m = mylake_params.theta_m;
dz = mylake_params.dz;
floculation_switch = mylake_params.floculation_switch;

Tz = mylake_temp_results.Tz;
DayFrac = mylake_temp_results.DayFrac;
lambdaz_wtot = mylake_temp_results.lambdaz_wtot;
H_sw_z = mylake_temp_results.H_sw_z;

k_POP         = sediment_params.k_POP * mylake_params.wc_factor;
k_POC         = sediment_params.k_POC * mylake_params.wc_factor;
k_DOP        = sediment_params.k_DOP * mylake_params.wc_factor;
k_DOC        = sediment_params.k_DOC * mylake_params.wc_factor;
k_POP_q10       = k_POP .* Q10.^((Tz-T_ref)/10);
k_POC_q10       = k_POC .* Q10.^((Tz-T_ref)/10);
k_DOP_q10       = k_DOP .* Q10.^((Tz-T_ref)/10);
k_DOC_q10       = k_DOC .* Q10.^((Tz-T_ref)/10);

% MyLake "old" chemistry:
% Conversion of units to "per year" and umoles
dop_twty_y = dop_twty*365; % In MyLake parameters with set to 0;

%Algal parameters:
g_twty_y = g_twty*365;
m_twty_y = m_twty*365;
P_half_molar = P_half/30973.762;
N_half_molar = N_half/14006.7;
Si_half_molar = Si_half/60840;
% DOP: ( we treed this as bio reaction, see rate "Rc")
R_dDOP = 0; % dop_twty_y .* DOP .* theta_m.^(Tz-20);  %Mineralisation to P

%Chl:
Nitrogen = NO3 + NH4 + 1e-16;
N_frac = NO3./Nitrogen;
N_frac = N_frac .* (N_frac > 0);

Growth_bioz = zeros(length(Pz), length(g_twty_y));
Loss_bioz = zeros(length(Pz), length(g_twty_y));
R_bioz = zeros(length(Pz), length(g_twty_y));
for alg=1:length(g_twty_y)
    
    if mylake_params.N_limited(alg) == 1 & mylake_params.Si_limited(alg)  == 0 % N-Limited and P-limited
        %Growth_bioz(:, alg)=g_twty_y(alg)*theta_m.^(Tz-20) .* (Pz./(P_half_molar(alg)+Pz)) .* (Nitrogen./(N_half_molar(alg)+Nitrogen)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z(:, alg); 0]);
        Growth_bioz(:, alg)=g_twty_y(alg)*theta_m.^(Tz-20) .* min([(Pz./(P_half_molar(alg)+Pz)),(Nitrogen./(N_half_molar(alg)+Nitrogen))]')'.* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z(:, alg); 0]);
        Loss_bioz(:, alg)=m_twty_y(alg)*theta_m.^(Tz-20);
        R_bioz(:, alg) = Growth_bioz(:, alg)-Loss_bioz(:, alg);
        
    elseif mylake_params.Si_limited(alg) == 1 & mylake_params.N_limited(alg) == 1 % Si Limited, N Limited & P-limited
        % Growth_bioz(:, alg)=g_twty_y(alg)*theta_m.^(Tz-20) .* (Pz./(P_half_molar(alg)+Pz)).* (Nitrogen./(N_half_molar(alg)+Nitrogen)) .* (Si./(Si_half_molar(alg)+Si)).* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z(:, alg); 0]);
        Growth_bioz(:, alg)=g_twty_y(alg)*theta_m.^(Tz-20) .* min([(Pz./(P_half_molar(alg)+Pz)),(Nitrogen./(N_half_molar(alg)+Nitrogen)),(Si./(Si_half_molar(alg)+Si))]')'.* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z(:, alg); 0]);
        Loss_bioz(:, alg)=m_twty_y(alg)*theta_m.^(Tz-20);
        R_bioz(:, alg) = Growth_bioz(:, alg)-Loss_bioz(:, alg);
        
    else % only P limited
        Growth_bioz(:, alg)=g_twty_y(alg)*theta_m.^(Tz-20) .* (Pz./(P_half_molar(alg)+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z(:, alg); 0]);
        Loss_bioz(:, alg)=m_twty_y(alg)*theta_m.^(Tz-20);
        R_bioz(:, alg) = Growth_bioz(:, alg)-Loss_bioz(:, alg);
        
    end
%     if DayFrac > 0.6
%        'DF in run_wc_chem'; 
%     end
end

%flocculation
if (floculation_switch==1) %Fokema
    dfloc_DOC = 0.030 .* DOC;
else
    dfloc_DOC = 0;
end

% TP generalise for N algal species
R_dChlz_growth =  zeros(length(Chl{1}), length(g_twty_y));
for ii = 1: length(g_twty_y)
    R_dChlz_growth(:, ii) =  (Chl{ii}+1e-5) .* Growth_bioz(:, ii) - Chl{ii} .* Loss_bioz(:, ii); % 1e-5 to insure having spores in water-column with 0 input from catchment.
end
%Oxygen production in phytoplankton growth
R_Alg_tot_growth = sum(R_dChlz_growth,2);

% New chemistry

tot_FeOH3 = PP + Fe3;

f_O2    = O2 ./  (Km_O2 + O2) ;
f_NO3   = NO3 ./  (Km_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2) ;
f_FeOH3 = tot_FeOH3 ./  (Km_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2) ;
f_SO4 = SO4 ./ (Km_SO4 + SO4) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
f_CH4 = 1 - f_O2 - f_NO3 - f_FeOH3 - f_SO4;
f_CH4 = f_CH4.*(f_CH4>0);
Sum_H2S = H2S + HS;

part_PO4ads_tot_Fe = PP ./ (tot_FeOH3+1e-16); % Avoid division by zero

R1a =  0; % k_POP_q10  .* Chl .* f_O2 * accel;
R1b =  0; % k_POP_q10  .* Cz .* f_O2 * accel;
R1c =  k_DOP_q10  .* DOP .* f_O2 .* accel;
R1d =  k_DOC_q10 .* DOC .* f_O2 .* accel;
R1e =  k_POP_q10 .* POP .* f_O2 .* accel;
R1f =  k_POC_q10 .* POC .* f_O2 .* accel;

R2a =  0; % k_POP_q10  .* Chl .* f_NO3 .* accel;
R2b =  0; % k_POP_q10  .* Cz .* f_NO3 .* accel;
R2c =  k_DOP_q10  .* DOP .* f_NO3 .* accel;
R2d =  k_DOC_q10 .* DOC .* f_NO3 .* accel;
R2e =  k_POP_q10 .* POP .* f_NO3 .* accel;
R2f =  k_POC_q10 .* POC .* f_NO3 .* accel;

R3a_Fe =  (1 - part_PO4ads_tot_Fe) .* 0; % k_POP_q10  .* Chl .* f_FeOH3;
R3b_Fe =  (1 - part_PO4ads_tot_Fe) .* 0; % k_POP_q10  .* Cz .* f_FeOH3;
R3c_Fe =  (1 - part_PO4ads_tot_Fe) .* k_DOP_q10  .* DOP .* f_FeOH3;
R3d_Fe =  (1 - part_PO4ads_tot_Fe) .* k_DOC_q10 .* DOC .* f_FeOH3;
R3e_Fe =  (1 - part_PO4ads_tot_Fe) .* k_POP_q10 .* POP .* f_FeOH3;
R3f_Fe =  (1 - part_PO4ads_tot_Fe) .* k_POC_q10 .* POC .* f_FeOH3;
R3a_P =  part_PO4ads_tot_Fe .* 0; % k_POP_q10  .* Chl .* f_FeOH3;
R3b_P =  part_PO4ads_tot_Fe .* 0; % k_POP_q10  .* Cz .* f_FeOH3;
R3c_P =  part_PO4ads_tot_Fe .* k_DOP_q10  .* DOP .* f_FeOH3;
R3d_P =  part_PO4ads_tot_Fe .* k_DOC_q10 .* DOC .* f_FeOH3;
R3e_P =  part_PO4ads_tot_Fe .* k_POP_q10 .* POP .* f_FeOH3;
R3f_P =  part_PO4ads_tot_Fe .* k_POC_q10 .* POC .* f_FeOH3;
R3a = R3a_P + R3a_Fe;
R3b = R3b_P + R3b_Fe;
R3c = R3c_P + R3c_Fe;
R3d = R3d_P + R3d_Fe;
R3e = R3e_P + R3e_Fe;
R3f = R3f_P + R3f_Fe;

R5a =  0; % k_POP_q10  .* Chl .* f_SO4;
R5b =  0; % k_POP_q10  .* Cz .* f_SO4;
R5c =  k_DOP_q10  .* DOP .* f_SO4;
R5d =  k_DOC_q10 .* DOC .* f_SO4;
R5e =  k_POP_q10 .* POP .* f_SO4;
R5f =  k_POC_q10 .* POC .* f_SO4;

R6a =  0; % k_POP_q10  .* Chl .* f_CH4;
R6b =  0; % k_POP_q10  .* Cz .* f_CH4;
R6c =  k_DOP_q10  .* DOP .* f_CH4;
R6d =  k_DOC_q10 .* DOC .* f_CH4;
R6e =  k_POP_q10 .* POP .* f_CH4;
R6f =  k_POC_q10 .* POC .* f_CH4;

Ra  = R1a+R2a+R3a+R5a+R6a;
Rb  = R1b+R2b+R3b+R5b+R6b;
Rc  = R1c+R2c+R3c+R5c+R6c;
Rd  = R1d+R2d+R3d+R5d+R6d;
Re  = R1e+R2e+R3e+R5e+R6e;
Rf  = R1f+R2f+R3f+R5f+R6f;

R1  = R1a+R1b+R1c+R1d+R1e+R1f;
R2  = R2a+R2b+R2c+R2d+R2e+R2f;
R3  = R3a+R3b+R3c+R3d+R3e+R3f;
R5  = R5a+R5b+R5c+R5d+R5e+R5f;
R6  = R6a+R6b+R6c+R6d+R6e+R6f;

R11  = k_tsox .* O2 .* Sum_H2S;
R12  = k_tS_Fe .* Fe3 .* Sum_H2S;
R13  = k_Feox .* Fe2 .* O2;
% NOTE: Due to the reaction is too fast and could cause overshooting:
% we need to make this check if R*dt > Conc of source:
% R13 = (R13.*dt < Fe2).*R13 + (R13.*dt > Fe2).* Fe2 ./ (dt) * 0.5;
% R13 = (R13.*dt < O2).*R13 + (R13.*dt > O2).* O2 ./ (dt) * 0.5;

% R14  = k_amox .* O2 ./ (Km_oxao + O2) .* NH4 ./ (Km_amao + NH4);   % NOTE: Doesnt work - Highly unstable.
R14 = k_amox  .* O2 .* NH4;
% R14 = (R14.*dt < NH4).*R14 + (R14.*dt > NH4).* NH4 ./ (dt) * 0.5;
% R14 = (R14.*dt < O2).*R14 + (R14.*dt > O2).* O2 ./ (dt) * 0.5;


R15 = k_ch4_o2 .* CH4aq .* O2;
R16 = k_ch4_so4 .* CH4aq .* SO4;

CH4_solubility = Kh_CH4*(1+z.*0.1);
CH4_over_sat = CH4aq - CH4_solubility;
R17 = k_ch4_dis .* CH4_over_sat .* (CH4_over_sat > 0);

CO2_solubility = Kh_CO2*(1+z.*0.1);
CO2_over_sat = CO2 - CO2_solubility;

R21a = 0; %k_oms .* Sum_H2S .* Chlz;
R21b = 0; %k_oms .* Sum_H2S .* Cz;
R21c = k_oms .* Sum_H2S .* DOP;
R21d = k_oms .* Sum_H2S .* DOC;
R21e = k_oms .* Sum_H2S .* POP;
R21f = k_oms .* Sum_H2S .* POC;

S0 = 0;
S8 = 0;
R22a = k_Spre * S0;
R22b = k_Sdis .* S8;

FeS2 = 0;
R23 = k_fes2ox .* FeS2 .* O2;
R24 = k_fespre .* FeS .* S0;
R25 = k_fesox * O2 .* FeS;
R26 = k_fes2pre .* FeS .* Sum_H2S;

H3O = 10.^mylake_temp_results.pHz*1e3;
Sat_FeS = Fe2*1e-3 .* Sum_H2S*1e-3 ./ (H3O*1e-3+1e-16).^2 ./ K_FeS;
R27a = k_fe_pre .* (Sat_FeS-1) .* (Sat_FeS > 1);
R27b  =  k_fe_dis .* FeS .* (1-Sat_FeS) .* (Sat_FeS < 1); %

Sat_CC = Ca2*1e-3 .* CO3*1e-3 / K_CC;
R28a = k_CCpre .* (Sat_CC-1) .* (Sat_CC > 1);
R28b  =  k_CCdis .* CaCO3 .* (1-Sat_CC) .* (Sat_CC < 1); %

FeCO3 =0;
Sat_FC = Fe2*1e-3 .* CO3*1e-3 / K_FC;
R29a = k_FCpre .* (Sat_FC-1) .* (Sat_FC > 1);
R29b  =  k_FCdis .* FeCO3 .* (1-Sat_FC) .* (Sat_FC < 1); %

% Phosphorus

R31a = k_pdesorb_a * Fe3 .* Pz;
% R31b = 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx2*R3c_P + Cx3*R3d_P + Cx1*R3f_P); % f_pfe .* (4 * R3 + 2 * R12);
% R31b = (R31b.*dt < PO4adsa).*R31b + (R31b.*dt > PO4adsa).* PO4adsa ./ (dt) * 0.5;

FeOOH = 0;
R32a = k_pdesorb_b * FeOOH .* Pz;
% R32b = 4 * (Cx2*R4a_P + Cx3*R4b_P + Cx2*R4c_P + Cx3*R4d_P + Cx1*R4f_P);
% R32b = (R32b.*dt < PO4adsb).*R32b + (R32b.*dt > PO4adsb).* PO4adsb ./ (dt) * 0.5;

Fe3PO42 = 0;
Sat_viv = (Pz*1e-3).^2.*(Fe2*1e-3).^3./(H3O*1e-3+1e-16).^2./K_viv;
R33a = k_viv_pre .* (Sat_viv-1) .* (Sat_viv > 1);
R33b = k_viv_dis .* Fe3PO42 .* (1 - Sat_viv) .* (Sat_viv < 1);

Ca3PO42 = 0;
Sat_apa = (Pz*1e-3).^2.*(Ca2*1e-3).^3./(H3O*1e-3+1e-16).^2./K_apa;
R34a = k_apa_pre .* (Sat_apa-1) .* (Sat_apa > 1);
R34b = k_apa_dis .* Ca3PO42 .* (1 - Sat_apa) .* (Sat_apa < 1);

% R35 disabled now (no solid species Al=PO4)
R35a = 0; % k_pdesorb_c .* PO4 .* AlOH3;
R35b = 0;

% saving rates?
if sediment_params.rate_estimator_switch
    r.R_dDOP  = R_dDOP; r.Growth_bioz = Growth_bioz; r.Loss_bioz = Loss_bioz; r.R_bioz  = R_bioz; r.dfloc_DOC  = dfloc_DOC; r.dfloc_DOC  = dfloc_DOC; r.R_dChlz_growth  = R_dChlz_growth; r.R_Alg_tot_growth  = R_Alg_tot_growth; r.tot_FeOH3  = tot_FeOH3; r.f_O2     = f_O2; r.f_NO3    = f_NO3; r.f_FeOH3  = f_FeOH3; r.f_SO4  = f_SO4; r.f_CH4  = f_CH4; r.f_CH4  = f_CH4; r.Sum_H2S  = Sum_H2S; r.part_PO4ads_tot_Fe  = part_PO4ads_tot_Fe; r.R1a  = R1a; r.R1b  = R1b; r.R1c  = R1c; r.R1d  = R1d; r.R1e  = R1e; r.R1f  = R1f; r.R2a  = R2a; r.R2b  = R2b; r.R2c  = R2c; r.R2d  = R2d; r.R2e  = R2e; r.R2f  = R2f; r.R3a_Fe  = R3a_Fe; r.R3b_Fe  = R3b_Fe; r.R3c_Fe  = R3c_Fe; r.R3d_Fe  = R3d_Fe; r.R3e_Fe  = R3e_Fe; r.R3f_Fe  = R3f_Fe; r.R3a_P  = R3a_P; r.R3b_P  = R3b_P; r.R3c_P  = R3c_P; r.R3d_P  = R3d_P; r.R3e_P  = R3e_P; r.R3f_P  = R3f_P; r.R3a  = R3a; r.R3b  = R3b; r.R3c  = R3c; r.R3d  = R3d; r.R3e  = R3e; r.R3f  = R3f; r.R5a  = R5a; r.R5b  = R5b; r.R5c  = R5c; r.R5d  = R5d; r.R5e  = R5e; r.R5f  = R5f; r.R6a  = R6a; r.R6b  = R6b; r.R6c  = R6c; r.R6d  = R6d; r.R6e  = R6e; r.R6f  = R6f; r.Ra   = Ra; r.Rb   = Rb; r.Rc   = Rc; r.Rd   = Rd; r.Re   = Re; r.Rf   = Rf; r.R1   = R1; r.R2   = R2; r.R3   = R3; r.R5   = R5; r.R6   = R6; r.R11   = R11; r.R12   = R12; r.R13   = R13; r.R14  = R14; r.R15  = R15; r.R16  = R16; r.CH4_solubility  = CH4_solubility; r.CH4_over_sat  = CH4_over_sat; r.R17  = R17; r.CO2_solubility  = CO2_solubility; r.CO2_over_sat  = CO2_over_sat; r.R21a  = R21a; r.R21b  = R21b; r.R21c  = R21c; r.R21d  = R21d; r.R21e  = R21e; r.R21f  = R21f; r.S0  = S0; r.S8  = S8; r.R22a  = R22a; r.R22b  = R22b; r.FeS2  = FeS2; r.R23  = R23; r.R24  = R24; r.R25  = R25; r.R26  = R26; r.H3O  = H3O; r.Sat_FeS  = Sat_FeS; r.R27a  = R27a; r.R27b   = R27b; r.Sat_CC  = Sat_CC; r.R28a  = R28a; r.R28b    = R28b; r.Sat_FC  = Sat_FC; r.R29a  = R29a; r.R29b   = R29b; r.R31a   = R31a; r.R32a  = R32a; r.Sat_viv  = Sat_viv; r.R33a  = R33a; r.R33b  = R33b; r.Ca3PO42  = Ca3PO42; r.Sat_apa  = Sat_apa; r.R34a  = R34a; r.R34b  = R34b; r.R35a  = R35a; r.R35b  = R35b;
else
    r.R1=0;
end

dcdt(:,1)  = -0.25 * R13  - R15 - 2 * R14  - 2* R11 - (Cx1*R1a + Cx1*R1b + Cx2*R1c + Cx3*R1d+ Cx2*R1e+ Cx3*R1f) + Cx1 * R_Alg_tot_growth - 3 * R25  - 5*R23 ; % O2
dcdt(:,2)  = -Rd - Cx3*R21d - dfloc_DOC;% DOC
dcdt(:,3)  = - 0.8*(Cx1*R2a + Cx1*R2b + Cx2*R2c + Cx3*R2d+ Cx2*R2e + Cx3*R2f) + R14 - Ny1 * R_Alg_tot_growth; % NO3z
dcdt(:,4)  = - 4*(Cx1*R3a_Fe + Cx1*R3b_Fe + Cx2*R3c_Fe + Cx3*R3d_Fe+ Cx2*R3e_Fe+ Cx3*R3f_Fe) - 2*R12  + R13 - R31a; % Fe3
dcdt(:,5)  = - 0.5*(Cx1*R5a + Cx1*R5b + Cx2*R5c + Cx3*R5d+ Cx2*R5e+ Cx3*R5f) + R11 - R16 + 2*R23; % SO4
dcdt(:,6)  =  (Ny1 * Ra + Ny1 * Rb + Ny2 * Rc + Ny3 * Rd+ Ny2 * Re+ Ny3 * Rf) - R14;% NH4
dcdt(:,7)  = 4*(Cx1*R3a + Cx1*R3b + Cx2*R3c + Cx3*R3d + Cx2*R3e+ Cx3*R3f) + 2*R12 - R13 + R27b - R27a - R29a + R29b -3*R33a + 3*R33b; % Fe2
dcdt(:,8)  =  - H2S./(Sum_H2S+1e-16).* (R11 + R12 + R27b + R27a + R26 + Cx1*R21a + Cx1*R21b + Cx2*R21c + Cx3*R21d + Cx2*R21e + Cx3*R21f);% H2S
dcdt(:,9) = 0.5*(Cx1*R5a + Cx1*R5b + Cx2*R5c + Cx3*R5d+ Cx2*R5e+ Cx3*R5f) - HS./(Sum_H2S+1e-16).* (R11 + R12 + R27b + R27a + R26 + Cx1*R21a + Cx1*R21b + Cx2*R21c + Cx3*R21d + Cx2*R21e + Cx3*R21f);% HS
dcdt(:,10) = (Pz1 * Ra + Pz1 * Rb + Pz2 * Rc + Pz3 * Rd+ Pz2 * Re+ Pz3 * Rf) + 4*(Cx1*R3a_P + Cx1*R3b_P + Cx2*  Cx3*R3d_P+ Cx2*R3e_P+ Cx3*R3f_P) - R31a - R32a - R35a + R35b - 2 * R34a + 2 * R34b -2*R33a + 2*R33b + R_dDOP - Pz1 * R_Alg_tot_growth;% Pz
dcdt(:,11) = -R35a ;% Al3
dcdt(:,12) = R31a - 4*(Cx1*R3a_P + Cx1*R3b_P + Cx2*  Cx3*R3d_P+ Cx2*R3e_P+ Cx3*R3f_P);% PP or PO4adsa
dcdt(:,13) = -3*R34a + 3*R34b - R28a + R28b;% Ca2
dcdt(:,14) = ((2*R13 + 2*R14 + R15)+ (Cx1 - Ny1 + 2*Pz1)*R1a + (Cx1 - Ny1 + 2*Pz1)*R1b + (Cx2 - Ny2 + 2*Pz2)*R1c + (Cx3 - Ny3 + 2*Pz3)*R1d + (Cx2 - Ny2 + 2*Pz2)*R1e + (Cx2 - Ny2 + 2*Pz2)*R1f + (0.2*Cx1 - Ny1 + 2*Pz1)*R2a + (0.2*Cx1 - Ny1 + 2*Pz1)*R2b + (0.2*Cx2 - Ny2 + 2*Pz2)*R2c + (0.2*Cx3 - Ny3 + 2*Pz3)*R2d + (0.2*Cx2 - Ny2 + 2*Pz2)*R2e + (0.2*Cx3 - Ny3 + 2*Pz3)*R2f + (0.5 * Cx1 - Ny1 - 2*Pz1)*R6a + (0.5 * Cx1 - Ny1 - 2*Pz1)*R6b + (0.5*Cx2 - Ny2 + 2*Pz2)*R6c + (0.5*Cx3 - Ny3 + 2*Pz3)*R6d + (0.5*Cx2 - Ny2 + 2*Pz2)*R6e + (0.5 * Cx3 - Ny3 - 2*Pz3)*R6f).*(CO2_over_sat<0)- R16 -(Ny1-2*Pz1)*R_Alg_tot_growth  - (7*Cx1 + Ny1 + 2*Pz1)*R3a - (7*Cx1 + Ny1 + 2*Pz1)*R3b - (7*Cx2 + Ny2 + 2*Pz2)*R3c - (7*Cx3 + Ny3 + 2*Pz3)*R3d - (7*Cx2 + Ny2 + 2*Pz2)*R3e - (7*Cx3 + Ny3 + 2*Pz3)*R3f - (Ny1 - 2*Pz1)*R5a - (Ny1 - 2*Pz1)*R5b - (Ny2 - 2*Pz2)*R5c - (Ny3 - 2*Pz3)*R5d - (Ny2 - 2*Pz2)*R5e - (Ny3 - 2*Pz3)*R5f ;% CO2
dcdt(:,15) = -Rc - Cx2*R21c - R_dDOP;% DOP
dcdt(:,16) = -Rf - Cx3*R21f + dfloc_DOC; % POC
dcdt(:,17) = - Re - Cx2*R21e;% POP
dcdt(:,18) = 0.5*(Cx1*R6a + Cx1*R6b + Cx2*R6c + Cx3*R6d+ Cx2*R6e + Cx3*R6f).* (CH4_over_sat < 0) -R15 - R16 + R17;% CH4aq
dcdt(:,19) = 0.5*(Cx1*R6a + Cx1*R6b + Cx2*R6c + Cx3*R6d+ Cx2*R6e + Cx3*R6f).* (CH4_over_sat > 0) - R17;% CH4g
dcdt(:,20) = (Ny1+2*Pz1)*R_Alg_tot_growth - (Ny1 - 2*Pz1)*R1a - (Ny1 - 2*Pz1)*R1b - (Ny2 - 2*Pz2)*R1c - (Ny3 - 2*Pz3)*R1d - (Ny2 - 2*Pz2)*R1e - (Ny3 - 2*Pz3)*R1f + (0.8*Cx1 + Ny1 - 2*Pz1)*R2a + (0.8*Cx1 + Ny1 - 2*Pz1)*R2b + (0.8*Cx2 + Ny2 - 2*Pz2)*R2c + (0.8*Cx3 + Ny3 - 2*Pz3)*R2d + (0.8*Cx2 + Ny2 - 2*Pz2)*R2e + (0.8*Cx3 + Ny3 - 2*Pz3)*R2f + (8*Cx1+Ny1-2*Pz1)*R3a + (8*Cx1+Ny1-2*Pz1)*R3b + (8*Cx2+Ny2-2*Pz2)*R3c + (8*Cx3+Ny3-2*Pz3)*R3d + (8*Cx2+Ny2-2*Pz2)*R3e + (8*Cx3+Ny3-2*Pz3)*R3f + (Cx1+Ny1-2*Pz1)*R5a + (Cx1+Ny1-2*Pz1)*R5b + (Cx2+Ny2-2*Pz2)*R5c + (Cx3+Ny3-2*Pz3)*R5d + (Cx2+Ny2-2*Pz2)*R5e + (Cx3+Ny3-2*Pz3)*R5f + (Ny1-2*Pz1)*R6a + (Ny1-2*Pz1)*R6b + (Ny2-2*Pz2)*R6c + (Ny3-2*Pz3)*R6d + (Ny2-2*Pz2)*R6e + (Ny3-2*Pz3)*R6f   - 2*R13 - 2*R14 + 2*R16; % HCO3
dcdt(:,21) = 0; % CO3
dcdt(:,22) = R28a - R28b; % CaCO3
dcdt(:,23) = - R24 - 4*R25 -R26 + R27a - R27b; % FeS
%TP
if sum(mylake_params.Si_limited)>0;
    % P usage from diatoms only
    PUsi = (Pz1 * Ra + Pz1 * Rb + Pz2 * Rc + Pz3 * Rd+ Pz2 * Re+ Pz3 * Rf) + 4*(Cx1*R3a_P + Cx1*R3b_P + Cx2*  Cx3*R3d_P+ Cx2*R3e_P+ Cx3*R3f_P) - R31a - R32a - R35a + R35b - 2 * R34a + 2 * R34b -2*R33a + 2*R33b + R_dDOP - Pz1 * sum(R_dChlz_growth(:, logical(mylake_params.Si_limited)),2);
    % TP multiply P use by molar ratio of Si:P
    dcdt(:,24) =  PUsi * 20; %Si  diatoms Redfield molar C:N:P:Si ratios 106:16:1:20 [Redfield et al. , 1963; Conley et al. , 1989].
else
    dcdt(:,24) = 0;
end

% TP Generalise for multiple species
for ii = 1:length(g_twty_y)
    dcdt(:,24+ii) = -Rb - Cx1*R21b + R_dChlz_growth(:, ii);% Chl_N
end
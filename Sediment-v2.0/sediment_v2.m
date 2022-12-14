function [ sediment_bioirrigation_fluxes, sediment_transport_fluxes, sediment_concentrations, sediment_additional_results] = sediment_v2(sediment_concentrations, sediment_params, sediment_matrix_templates, sediment_bc)
  % SEDIMENTS This function models the chemical process in the sediment

  % model domain:
  years = sediment_params.years; %1 day
  ts    = sediment_params.ts; % time step
  alfax = sediment_params.alfax;
  D_O2  = sediment_params.D_O2;
  D_NO3 = sediment_params.D_NO3;
  D_SO4 = sediment_params.D_SO4;
  D_NH4 = sediment_params.D_NH4;
  D_Fe2 = sediment_params.D_Fe2;
  D_H2S = sediment_params.D_H2S;
  D_S0  = sediment_params.D_S0;
  D_PO4 = sediment_params.D_PO4;
  D_Ca2 = sediment_params.D_Ca2;
  D_HS  = sediment_params.D_HS;
  D_DOP = sediment_params.D_DOP;
  D_DOC = sediment_params.D_DOC;
  D_CO2 = sediment_params.D_CO2;
  D_CO2g = sediment_params.D_CO2g;
  D_HCO3 = sediment_params.D_HCO3;
  D_CO3 = sediment_params.D_CO3;
  D_NH3 = sediment_params.D_NH3;
  D_CH4aq = sediment_params.D_CH4aq;
  D_CH4g = sediment_params.D_CH4g;
  phi    = sediment_params.phi;
  x  = sediment_params.x;
  SO4_bottom_fx = sediment_params.SO4_bottom_fx;


  dx = x(2)-x(1);
  xz = x'/100;
  % time domain:
  t     = 0:ts:years; % years
  m     = size(t,2); %steps in time
  dt    = t(2)-t(1);
  n = sediment_params.n;


  % Preallocation:
  O2 = zeros(n, m);
  POP = zeros(n, m);
  POC = zeros(n, m);
  NO3 = zeros(n, m);
  FeOH3 = zeros(n, m);
  SO4 = zeros(n, m);
  Fe2 = zeros(n, m);
  FeOOH = zeros(n, m);
  FeS = zeros(n, m);
  S0 = zeros(n, m);
  PO4 = zeros(n, m);
  S8 = zeros(n, m);
  FeS2 = zeros(n, m);
  AlOH3 = zeros(n, m);
  PO4adsa = zeros(n, m);
  PO4adsb = zeros(n, m);
  Ca2 = zeros(n, m);
  Ca3PO42 = zeros(n, m);
  OMS = zeros(n, m);
  CaCO3 = zeros(n, m);
  CO2 = zeros(n, m);
  CO3 = zeros(n, m);
  HCO3 = zeros(n, m);
  NH3 = zeros(n, m);
  NH4 = zeros(n, m);
  HS = zeros(n, m);
  H2S = zeros(n, m);
  CO2g = zeros(n, m);
  DOP = zeros(n, m);
  DOC = zeros(n, m);
  Chl = zeros(n, m);
  CH4aq = zeros(n, m);
  CH4g = zeros(n, m);
  FeCO3 = zeros(n, m);
  Fe3PO42 = zeros(n, m);
  PO4adsc = zeros(n, m);

  O2(:,1) = sediment_concentrations.O2;
  POP(:,1) = sediment_concentrations.POP;
  POC(:,1) = sediment_concentrations.POC;
  NO3(:,1) = sediment_concentrations.NO3;
  FeOH3(:,1) = sediment_concentrations.FeOH3;
  SO4(:,1) = sediment_concentrations.SO4;
  Fe2(:,1) = sediment_concentrations.Fe2;
  FeOOH(:,1) = sediment_concentrations.FeOOH;
  FeS(:,1) = sediment_concentrations.FeS;
  S0(:,1) = sediment_concentrations.S0;
  PO4(:,1)  = sediment_concentrations.PO4;
  S8(:,1) = sediment_concentrations.S8;
  FeS2(:,1)  = sediment_concentrations.FeS2;
  AlOH3(:,1) = sediment_concentrations.AlOH3;
  PO4adsa(:,1) = sediment_concentrations.PO4adsa;
  PO4adsb(:,1) = sediment_concentrations.PO4adsb;
  Ca2(:,1) = sediment_concentrations.Ca2;
  Ca3PO42(:,1) = sediment_concentrations.Ca3PO42;
  OMS(:,1) = sediment_concentrations.OMS;
  CaCO3(:,1) = sediment_concentrations.CaCO3;
  CO2(:,1) = sediment_concentrations.CO2;
  CO3(:,1) = sediment_concentrations.CO3;
  HCO3(:,1) = sediment_concentrations.HCO3;
  NH3(:,1) = sediment_concentrations.NH3;
  NH4(:,1) = sediment_concentrations.NH4;
  HS(:,1) = sediment_concentrations.HS;
  H2S(:,1) = sediment_concentrations.H2S;
  CO2g(:,1) = sediment_concentrations.CO2g;
  DOP(:,1) = sediment_concentrations.DOP;
  DOC(:,1) = sediment_concentrations.DOC;
  Chl(:,1) = sediment_concentrations.Chl;
  CH4aq(:,1) = sediment_concentrations.CH4aq;
  CH4g(:,1) = sediment_concentrations.CH4g;
  FeCO3(:,1) = sediment_concentrations.FeCO3;
  Fe3PO42(:,1) = sediment_concentrations.Fe3PO42;
  PO4adsc(:,1) = sediment_concentrations.PO4adsc;


  H3O = ones(n, 1).* sediment_bc.H3O_c;





  % Allocation of the memory and formation of template matrix:
  % =======================================================================================================

  % Solid species: row #1 of the cell "sediment_matrix_templates" is the solid template matrix
  [Solid_AL, Solid_AR]= sediment_matrix_templates{1,1:2};

  % Solute species:
  [O2_AL, O2_AR] = sediment_matrix_templates{2,1:2};
  [NO3_AL, NO3_AR] = sediment_matrix_templates{3,1:2};
  [SO4_AL, SO4_AR] = sediment_matrix_templates{4,1:2};
  [NH4_AL, NH4_AR] = sediment_matrix_templates{5,1:2};
  [Fe2_AL, Fe2_AR] = sediment_matrix_templates{6,1:2};
  [H2S_AL, H2S_AR] = sediment_matrix_templates{7,1:2};
  [S0_AL, S0_AR] = sediment_matrix_templates{8,1:2};
  [PO4_AL, PO4_AR] = sediment_matrix_templates{9,1:2};
  [Ca2_AL, Ca2_AR] = sediment_matrix_templates{10,1:2};
  [HS_AL, HS_AR] = sediment_matrix_templates{11,1:2};
  [CO2_AL, CO2_AR] = sediment_matrix_templates{12,1:2};
  [CO3_AL, CO3_AR] = sediment_matrix_templates{13,1:2};
  [HCO3_AL, HCO3_AR] = sediment_matrix_templates{14,1:2};
  [NH3_AL, NH3_AR] = sediment_matrix_templates{15,1:2};
  [CO2g_AL, CO2g_AR] = sediment_matrix_templates{16,1:2};
  [DOP_AL, DOP_AR] = sediment_matrix_templates{17,1:2};
  [DOC_AL, DOC_AR] = sediment_matrix_templates{18,1:2};
  [CH4aq_AL, CH4aq_AR] = sediment_matrix_templates{19,1:2};
  [CH4g_AL, CH4g_AR] = sediment_matrix_templates{20,1:2};


  % Solving equations!!!
  % =========================================================================================================

  for i=2:m


    % =======================================================================================================
    % pH Module only once a day
    % =======================================================================================================
    if sediment_params.pH_algorithm ~= 0 & i == 2
      [H3O] = pH_module(sediment_params.pH_algorithm, H3O, CO2g(:,i-1), HCO3(:,i-1), CO2(:,i-1), CO3(:,i-1), NH3(:,i-1), NH4(:,i-1), HS(:,i-1), H2S(:,i-1), Fe2(:,i-1), Ca2(:,i-1), NO3(:,i-1), SO4(:,i-1), PO4(:,i-1), FeS(:,i-1), FeS2(:,i-1), FeOH3(:,i-1), FeOOH(:,i-1), Ca3PO42(:,i-1), PO4adsa(:,i-1), PO4adsb(:,i-1), sediment_bc.T, sediment_params);
    end

    Ct = CO2(:,i-1) + HCO3(:,i-1) + CO3(:,i-1);
    a = alpha(-log10(H3O)+3, sediment_params.aq_system.carb_acid.pKs);
    CO2(:,i-1) = Ct.*a(:,1);
    HCO3(:,i-1) = Ct.*a(:,2);
    CO3(:,i-1) = Ct.*a(:,3);

    Nt = NH3(:,i-1) + NH4(:,i-1);
    a = alpha(-log10(H3O)+3, sediment_params.aq_system.amonia.pKs);
    NH4(:,i-1) = Nt.*a(:,1);
    NH3(:,i-1) = Nt.*a(:,2);

    St = HS(:,i-1) + H2S(:,i-1);
    a = alpha(-log10(H3O)+3, sediment_params.aq_system.sulf.pKs);
    H2S(:,i-1) = St.*a(:,1);
    HS(:,i-1) = St.*a(:,2);

    % -log10(H3O)+3

    a = alpha(-log10(H3O)+3, sediment_params.aq_system.p_acid.pKs);
    sediment_params.HPO4 = a(:,3).*PO4(:,i-1);
    sediment_params.H2PO4 = a(:,2).*PO4(:,i-1);

    sediment_params.H3O = H3O;

    % =======================================================================================================
    % Solving Reaction eq-s
    % =======================================================================================================
    C0 = [O2(:,i-1), POP(:,i-1), POC(:,i-1), NO3(:,i-1), FeOH3(:,i-1), SO4(:,i-1), NH4(:,i-1), Fe2(:,i-1), FeOOH(:,i-1), H2S(:,i-1), HS(:,i-1), FeS(:,i-1), S0(:,i-1), PO4(:,i-1), S8(:,i-1), FeS2(:,i-1), AlOH3(:,i-1), PO4adsa(:,i-1), PO4adsb(:,i-1), Ca2(:,i-1), Ca3PO42(:,i-1), OMS(:,i-1), CaCO3(:,i-1), CO2(:,i-1), CO3(:,i-1), HCO3(:,i-1), NH3(:,i-1), CO2g(:,i-1), DOP(:,i-1), DOC(:,i-1), Chl(:,i-1), CH4aq(:,i-1), CH4g(:,i-1), FeCO3(:,i-1), Fe3PO42(:,i-1), PO4adsc(:,i-1)];

      if any(any(isnan(C0)))
          error('NaN')
      end



      int_method = 0;
      [C_new, rates(i-1)] = sediments_chemical_reactions_module(sediment_params, C0,dt, int_method);

      if any(any(isnan(C_new)))
          if sediment_params.pH_algorithm == 0
            error('NaN')
          else
            error('NaN, try to disable pH algorithm')
          end
      end

      O2(:,i)      = C_new(:,1);
      POP(:,i)      = C_new(:,2);
      POC(:,i)     = C_new(:,3);
      NO3(:,i)     = C_new(:,4);
      FeOH3(:,i)   = C_new(:,5);
      SO4(:,i)     = C_new(:,6);
      NH4(:,i)     = C_new(:,7);
      Fe2(:,i)     = C_new(:,8);
      FeOOH(:,i)   = C_new(:,9);
      H2S(:,i)     = C_new(:,10);
      HS(:,i)      = C_new(:,11);
      FeS(:,i)     = C_new(:,12);
      S0(:,i)      = C_new(:,13);
      PO4(:,i)     = C_new(:,14);
      S8(:,i)      = C_new(:,15);
      FeS2(:,i)    = C_new(:,16);
      AlOH3(:,i)   = C_new(:,17);
      PO4adsa(:,i) = C_new(:,18);
      PO4adsb(:,i) = C_new(:,19);
      Ca2(:,i)     = C_new(:,20);
      Ca3PO42(:,i) = C_new(:,21);
      OMS(:,i)     = C_new(:,22);
      CaCO3(:,i)   = C_new(:,23);
      CO2(:,i)     = C_new(:,24);
      CO3(:,i)     = C_new(:,25);
      HCO3(:,i)    = C_new(:,26);
      NH3(:,i)     = C_new(:,27);
      CO2g(:,i)   = C_new(:,28);
      DOP(:,i)   = C_new(:,29);
      DOC(:,i)   = C_new(:,30);
      Chl(:,i)   = C_new(:,31);
      CH4aq(:,i)   = C_new(:,32);
      CH4g(:,i)   = C_new(:,33);
      FeCO3(:,i)   = C_new(:,34);
      Fe3PO42(:,i)   = C_new(:,35);
      PO4adsc(:,i)   = C_new(:,36);


    % =======================================================================================================
    % Solving Transport eq-s
    % =======================================================================================================

      % dissolved species
      O2(:,i) = pde_solver_dirichlet(O2_AL, O2_AR, O2(:,i), sediment_bc.O2_c);
      NO3(:,i) = pde_solver_dirichlet(NO3_AL, NO3_AR, NO3(:,i), sediment_bc.NO3_c);
      NH4(:,i) = pde_solver_dirichlet(NH4_AL, NH4_AR, NH4(:,i), sediment_bc.NH4_c);
      Fe2(:,i) = pde_solver_dirichlet(Fe2_AL, Fe2_AR, Fe2(:,i), sediment_bc.Fe2_c);
      H2S(:,i) = pde_solver_dirichlet(H2S_AL, H2S_AR, H2S(:,i), sediment_bc.H2S_c);
      HS(:,i) = pde_solver_dirichlet(HS_AL, HS_AR, HS(:,i), sediment_bc.HS_c);
      S0(:,i) = pde_solver_dirichlet(S0_AL, S0_AR, S0(:,i), sediment_bc.S0_c);
      PO4(:,i) = pde_solver_dirichlet(PO4_AL, PO4_AR, PO4(:,i), sediment_bc.PO4_c);
      Ca2(:,i) = pde_solver_dirichlet(Ca2_AL, Ca2_AR, Ca2(:,i), sediment_bc.Ca2_c);
      CO2(:,i) = pde_solver_dirichlet(CO2_AL, CO2_AR, CO2(:,i), sediment_bc.CO2_c);
      CO3(:,i) = pde_solver_dirichlet(CO3_AL, CO3_AR, CO3(:,i), sediment_bc.CO3_c);
      HCO3(:,i) = pde_solver_dirichlet(HCO3_AL, HCO3_AR, HCO3(:,i), sediment_bc.HCO3_c);
      NH3(:,i) = pde_solver_dirichlet(NH3_AL, NH3_AR, NH3(:,i), sediment_bc.NH3_c);
      CO2g(:,i) = pde_solver_dirichlet(CO2g_AL, CO2g_AR, CO2g(:,i), sediment_bc.CO2g_c);
      DOP(:,i) = pde_solver_dirichlet(DOP_AL, DOP_AR, DOP(:,i), sediment_bc.DOP_c);
      DOC(:,i) = pde_solver_dirichlet(DOC_AL, DOC_AR, DOC(:,i), sediment_bc.DOC_c);
      CH4aq(:,i) = pde_solver_dirichlet(CH4aq_AL, CH4aq_AR, CH4aq(:,i), sediment_bc.CH4aq_c);
      CH4g(:,i) = pde_solver_dirichlet(CH4g_AL, CH4g_AR, CH4g(:,i), sediment_bc.CH4g_c);


      % solid species
      POP(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, POP(:,i), sediment_bc.POP_fx, sediment_params.solid_flux_coef);
      POC(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, POC(:,i), sediment_bc.POC_fx, sediment_params.solid_flux_coef);
      FeOH3(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, FeOH3(:,i), sediment_bc.FeOH3_fx, sediment_params.solid_flux_coef);
      FeOOH(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, FeOOH(:,i), sediment_bc.FeOOH_fx, sediment_params.solid_flux_coef);
      FeS(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, FeS(:,i), sediment_bc.FeS_fx, sediment_params.solid_flux_coef);
      S8(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, S8(:,i), sediment_bc.S8_fx, sediment_params.solid_flux_coef);
      FeS2(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, FeS2(:,i), sediment_bc.FeS2_fx, sediment_params.solid_flux_coef);
      AlOH3(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, AlOH3(:,i), sediment_bc.AlOH3_fx, sediment_params.solid_flux_coef);
      PO4adsa(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, PO4adsa(:,i), sediment_bc.PO4adsa_fx, sediment_params.solid_flux_coef);
      PO4adsb(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, PO4adsb(:,i), sediment_bc.PO4adsb_fx, sediment_params.solid_flux_coef);
      Ca3PO42(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, Ca3PO42(:,i), sediment_bc.Ca3PO42_fx, sediment_params.solid_flux_coef);
      OMS(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, OMS(:,i), sediment_bc.OMS_fx, sediment_params.solid_flux_coef);
      CaCO3(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, CaCO3(:,i), sediment_bc.CaCO3_fx, sediment_params.solid_flux_coef);
      Chl(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, Chl(:,i), sediment_bc.Chl_fx, sediment_params.solid_flux_coef);
      FeCO3(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, FeCO3(:,i), sediment_bc.FeCO3_fx, sediment_params.solid_flux_coef);
      Fe3PO42(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, Fe3PO42(:,i), sediment_bc.Fe3PO42_fx, sediment_params.solid_flux_coef);
      PO4adsc(:,i) = pde_solver_neumann(Solid_AL, Solid_AR, PO4adsc(:,i), sediment_bc.PO4adsc_fx, sediment_params.solid_flux_coef);

      % custom boundary conditions
      SO4(:,i) = pde_solver_custom(SO4_AL, SO4_AR, SO4(:,i), sediment_bc.SO4_c, SO4_bottom_fx, sediment_params.SO4_flux_coef);

      % Estimate fluxes:
      sediment_bioirrigation_fluxes.O2(i-1)   = integrate_over_depth_2(bioirrigation(O2(:, i), alfax, phi), x);
      sediment_bioirrigation_fluxes.PO4(i-1)  = integrate_over_depth_2(bioirrigation(PO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.Fe2(i-1)  = integrate_over_depth_2(bioirrigation(Fe2(:, i)*sediment_params.Kd_fe2,  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NO3(i-1)  = integrate_over_depth_2(bioirrigation(NO3(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.NH4(i-1)  = integrate_over_depth_2(bioirrigation(NH4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.SO4(i-1)  = integrate_over_depth_2(bioirrigation(SO4(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.HS(i-1)  = integrate_over_depth_2(bioirrigation(HS(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.H2S(i-1)  = integrate_over_depth_2(bioirrigation(H2S(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOP(i-1) = integrate_over_depth_2(bioirrigation(DOP(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.DOC(i-1) = integrate_over_depth_2(bioirrigation(DOC(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CH4aq(i-1) = integrate_over_depth_2(bioirrigation(CH4aq(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CH4g(i-1) = integrate_over_depth_2(bioirrigation(CH4g(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CO2(i-1) = integrate_over_depth_2(bioirrigation(CO2(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.HCO3(i-1) = integrate_over_depth_2(bioirrigation(HCO3(:, i),  alfax,  phi), x);
      sediment_bioirrigation_fluxes.CO3(i-1) = integrate_over_depth_2(bioirrigation(CO3(:, i),  alfax,  phi), x);

      sediment_transport_fluxes.POP(i-1)          = -sediment_bc.POP_fx;
      sediment_transport_fluxes.Chl(i-1)          = -sediment_bc.Chl_fx;
      sediment_transport_fluxes.POC(i-1)          = -sediment_bc.POC_fx;
      sediment_transport_fluxes.FeOH3(i-1)        = -sediment_bc.FeOH3_fx;
      % sediment_transport_fluxes.FeOOH(i-1)        = -sediment_bc.FeOOH_fx;
      sediment_transport_fluxes.AlOH3(i-1)        = -sediment_bc.AlOH3_fx;
      sediment_transport_fluxes.PO4adsa(i-1)      = -sediment_bc.PO4adsa_fx;
      sediment_transport_fluxes.PO4adsb(i-1)      = -sediment_bc.PO4adsb_fx;
      sediment_transport_fluxes.FeCO3(i-1)        = -sediment_bc.FeCO3_fx;
      sediment_transport_fluxes.Fe3PO42(i-1)      = -sediment_bc.Fe3PO42_fx;
      sediment_transport_fluxes.FeS(i-1)          = -sediment_bc.FeS_fx;
      sediment_transport_fluxes.PO4adsc(i-1)      = -sediment_bc.PO4adsc_fx;
      % sediment_transport_fluxes.FeS2(i-1)          = -sediment_bc.FeS2_fx;
      % sediment_transport_fluxes.S8(i-1)          = -sediment_bc.S8_fx;
      % sediment_transport_fluxes.Ca3PO42(i-1)          = -sediment_bc.Ca3PO42_fx;
      % sediment_transport_fluxes.OMS(i-1)          = -sediment_bc.OMS_fx;
      sediment_transport_fluxes.CaCO3(i-1)          = -sediment_bc.CaCO3_fx;
      sediment_transport_fluxes.O2(i-1)           = top_sediment_diffusion_flux(O2(:, i), D_O2, dx, phi) + top_sediment_advection_flux(O2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.PO4(i-1)          = top_sediment_diffusion_flux(PO4(:, i), D_PO4, dx, phi) + top_sediment_advection_flux(PO4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.NO3(i-1)          = top_sediment_diffusion_flux(NO3(:, i), D_NO3, dx, phi) + top_sediment_advection_flux(NO3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.Fe2(i-1)          = top_sediment_diffusion_flux(Fe2(:, i)*sediment_params.Kd_fe2, D_Fe2, dx, phi) + top_sediment_advection_flux(Fe2(:, i)*sediment_params.Kd_fe2, sediment_params.w, phi);
      sediment_transport_fluxes.NH4(i-1)          = top_sediment_diffusion_flux(NH4(:, i), D_NH4, dx, phi) + top_sediment_advection_flux(NH4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.SO4(i-1)          = top_sediment_diffusion_flux(SO4(:, i), D_SO4, dx, phi) + top_sediment_advection_flux(SO4(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.HS(i-1)          = top_sediment_diffusion_flux(HS(:, i), D_HS, dx, phi) + top_sediment_advection_flux(HS(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.H2S(i-1)          = top_sediment_diffusion_flux(H2S(:, i), D_H2S, dx, phi) + top_sediment_advection_flux(H2S(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.DOP(i-1)         = top_sediment_diffusion_flux(DOP(:, i), D_DOP, dx, phi) + top_sediment_advection_flux(DOP(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.DOC(i-1)         = top_sediment_diffusion_flux(DOC(:, i), D_DOC, dx, phi) + top_sediment_advection_flux(DOC(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CH4aq(i-1)         = top_sediment_diffusion_flux(CH4aq(:, i), D_CH4aq, dx, phi) + top_sediment_advection_flux(CH4aq(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CH4g(i-1)         = top_sediment_diffusion_flux(CH4g(:, i), D_CH4g, dx, phi) + top_sediment_advection_flux(CH4g(:, i), sediment_params.w_CH4g, phi);  % NOTE: rising velocity of methane
      sediment_transport_fluxes.CO2(i-1)         = top_sediment_diffusion_flux(CO2(:, i), D_CO2, dx, phi) + top_sediment_advection_flux(CO2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.HCO3(i-1)         = top_sediment_diffusion_flux(HCO3(:, i), D_HCO3, dx, phi) + top_sediment_advection_flux(HCO3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CO3(i-1)         = top_sediment_diffusion_flux(CO3(:, i), D_CO3, dx, phi) + top_sediment_advection_flux(CO3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.Ca2(i-1)         = top_sediment_diffusion_flux(Ca2(:, i), D_Ca2, dx, phi) + top_sediment_advection_flux(Ca2(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.NH3(i-1)         = top_sediment_diffusion_flux(NH3(:, i), D_NH3, dx, phi) + top_sediment_advection_flux(NH3(:, i), sediment_params.w, phi);
      sediment_transport_fluxes.CO2g(i-1)         = top_sediment_diffusion_flux(CO2g(:, i), D_CO2g, dx, phi) + top_sediment_advection_flux(CO2g(:, i), sediment_params.w, phi);
  end

    % convert flux to [mg/m2/d]
    sediment_transport_fluxes.O2      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.O2), 31998);
    sediment_transport_fluxes.POP     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.POP), 30973.762);
    sediment_transport_fluxes.POC     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.POC), 12010.7);
    sediment_transport_fluxes.NO3     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.NO3), 62004);
    sediment_transport_fluxes.FeOH3   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeOH3), 106867);
    sediment_transport_fluxes.SO4     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.SO4), 96062);
    sediment_transport_fluxes.NH4     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.NH4), 18038);
    sediment_transport_fluxes.Fe2     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Fe2), 55845);
    % sediment_transport_fluxes.FeOOH   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeOOH), 88851.7);
    sediment_transport_fluxes.H2S     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.H2S), 34072);
    sediment_transport_fluxes.HS      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.HS), 33072);
    sediment_transport_fluxes.FeS     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeS), 87910);
    % sediment_transport_fluxes.S0      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.S0), 32065.0);
    sediment_transport_fluxes.PO4     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4), 30973.762);
    % sediment_transport_fluxes.S8      = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.S8), 256520);
    % sediment_transport_fluxes.FeS2    = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeS2), 119975.0);
    sediment_transport_fluxes.AlOH3   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.AlOH3), 78003.6);
    sediment_transport_fluxes.PO4adsa = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4adsa), 30973.762);
    sediment_transport_fluxes.PO4adsb = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4adsb), 30973.762);
    sediment_transport_fluxes.PO4adsc = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.PO4adsc), 30973.762);
    sediment_transport_fluxes.Ca2     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Ca2), 40078);
    % sediment_transport_fluxes.Ca3PO42 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Ca3PO42), 310176.7);
    % sediment_transport_fluxes.OMS     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.OMS), 12010.7);
    sediment_transport_fluxes.CaCO3   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CaCO3), 100086.9);
    sediment_transport_fluxes.CO2     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CO2), 44009.5);
    sediment_transport_fluxes.CO3     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CO3), 60008.9);
    sediment_transport_fluxes.HCO3    = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.HCO3), 61016.8);
    % sediment_transport_fluxes.NH3     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.NH3), 17038);
    sediment_transport_fluxes.CO2g    = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CO2g), 44009.5);
    sediment_transport_fluxes.DOP     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.DOP), 30973.762);
    sediment_transport_fluxes.DOC     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.DOC), 12010.7);
    sediment_transport_fluxes.Chl     = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Chl), 30973.762);
    sediment_transport_fluxes.CH4aq   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CH4aq), 16042.5);
    sediment_transport_fluxes.CH4g    = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.CH4g), 16042.5);
    sediment_transport_fluxes.FeCO3   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.FeCO3), 115853.9);
    sediment_transport_fluxes.Fe3PO42 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_transport_fluxes.Fe3PO42), 357477.7);





  sediment_bioirrigation_fluxes.O2   = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.O2), 31998);
  sediment_bioirrigation_fluxes.PO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.PO4), 30973.762);
  sediment_bioirrigation_fluxes.Fe2 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.Fe2), 55845);
  sediment_bioirrigation_fluxes.NO3  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NO3), 62004);
  sediment_bioirrigation_fluxes.NH4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.NH4), 18038);
  sediment_bioirrigation_fluxes.SO4  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.SO4), 33072);
  sediment_bioirrigation_fluxes.HS  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.HS), 33072);
  sediment_bioirrigation_fluxes.H2S  = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.H2S), 34072);
  sediment_bioirrigation_fluxes.DOP = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOP), 30973.762);
  sediment_bioirrigation_fluxes.DOC = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.DOC), 12010.7);
  sediment_bioirrigation_fluxes.CH4aq = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CH4aq), 16042.5);
  sediment_bioirrigation_fluxes.CH4g = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CH4g), 16042.5);
  sediment_bioirrigation_fluxes.CO2 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CO2), 44009.5);
  sediment_bioirrigation_fluxes.HCO3 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.HCO3), 61016.8);
  sediment_bioirrigation_fluxes.CO3 = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(mean(sediment_bioirrigation_fluxes.CO3), 60008.9);


  sediment_concentrations.O2 = O2(:,end);
  sediment_concentrations.POP = POP(:,end);
  sediment_concentrations.POC = POC(:,end);
  sediment_concentrations.NO3 = NO3(:,end);
  sediment_concentrations.FeOH3 = FeOH3(:,end);
  sediment_concentrations.SO4 = SO4(:,end);
  sediment_concentrations.NH4 = NH4(:,end);
  sediment_concentrations.Fe2 = Fe2(:,end);
  sediment_concentrations.Fe2d = Fe2(:,end)*sediment_params.Kd_fe2;
  sediment_concentrations.FeOOH = FeOOH(:,end);
  sediment_concentrations.H2S = H2S(:,end);
  sediment_concentrations.HS = HS(:,end);
  sediment_concentrations.FeS = FeS(:,end);
  sediment_concentrations.S0 = S0(:,end);
  sediment_concentrations.PO4 = PO4(:,end);
  sediment_concentrations.S8 = S8(:,end);
  sediment_concentrations.FeS2 = FeS2(:,end);
  sediment_concentrations.AlOH3 = AlOH3(:,end);
  sediment_concentrations.PO4adsa = PO4adsa(:,end);
  sediment_concentrations.PO4adsb = PO4adsb(:,end);
  sediment_concentrations.Ca2 = Ca2(:,end);
  sediment_concentrations.Ca3PO42 = Ca3PO42(:,end);
  sediment_concentrations.OMS = OMS(:,end);
  sediment_concentrations.H3O = H3O;
  sediment_concentrations.CaCO3 = CaCO3(:,end);
  sediment_concentrations.CO2 = CO2(:,end);
  sediment_concentrations.CO3 = CO3(:,end);
  sediment_concentrations.HCO3 = HCO3(:,end);
  sediment_concentrations.NH3 = NH3(:,end);
  sediment_concentrations.CO2g = CO2g(:,end);
  sediment_concentrations.DOP = DOP(:,end);
  sediment_concentrations.DOC = DOC(:,end);
  sediment_concentrations.Chl = Chl(:,end);
  sediment_concentrations.CH4aq = CH4aq(:,end);
  sediment_concentrations.CH4g = CH4g(:,end);
  sediment_concentrations.FeCO3 = FeCO3(:,end);
  sediment_concentrations.Fe3PO42 = Fe3PO42(:,end);
  sediment_concentrations.PO4adsc = PO4adsc(:,end);

  sediment_concentrations.pH = -log10(H3O)+3;


  % Estimate average rate during the day
  if sediment_params.rate_estimator_switch
    rate_number = fieldnames(rates);
    for i = 1:numel(rate_number)
        r.(rate_number{i}) = 0;
        for j=1:m-1
            r.(rate_number{i}) = r.(rate_number{i}) + rates(j).(rate_number{i});
        end
        r.(rate_number{i}) = r.(rate_number{i})/(m-1);
    end

    sediment_additional_results.rates = r;

  else
    sediment_additional_results.rates = false;
  end

    if any(isnan(sediment_transport_fluxes.O2))| any(isnan(sediment_bc.POP_fx))| any(isnan(sediment_bc.POC_fx))| any(isnan(sediment_bc.FeOH3_fx))| any(isnan(O2)) | any(isnan(POP)) | any(isnan(POC)) | any(isnan(NO3)) | any(isnan(FeOH3)) | any(isnan(SO4)) | any(isnan(NH4)) | any(isnan(Fe2)) | any(isnan(FeOOH)) | any(isnan(H2S)) | any(isnan(HS)) | any(isnan(FeS)) | any(isnan(S0)) | any(isnan(PO4)) | any(isnan(S8)) | any(isnan(FeS2)) | any(isnan(AlOH3)) | any(isnan(PO4adsa)) | any(isnan(PO4adsb)) | any(isnan(H3O)) | any(isnan(Ca2)) | any(isnan(Ca3PO42)) | any(isnan(OMS)) | any(isnan(CaCO3)) | any(isnan(HCO3)) | any(isnan(CO2)) | any(isnan(CO3)) | any(isnan(NH3)) | any(isnan(CO2g)) | any(isnan(Chl)) | any(isnan(CH4aq)) | any(isnan(CH4g))
      error('Breaking out of Sediments function: NaN values');
    end


end






function C_new = pde_solver_dirichlet(AL, AR, C_old, const_bc)
    C_old(1) = const_bc;
    temp = AR*C_old;
    C_new = AL\ temp;
    C_new = (C_new>0).*C_new;
end


function C_new = pde_solver_neumann(AL, AR, C_old, flux_bc, coef)
    temp = AR*C_old;
    temp(1) = temp(1) + flux_bc * coef;
    C_new = AL\ temp;
    C_new = (C_new>0).*C_new;
end

function C_new = pde_solver_custom(AL, AR, C_old, const_bc, flux_bc, coef)
    % solution for custom boundary conditions
    % Dirichlet at the top and Neumann at the bottom
    % in this example it is Sulfate in Vansjo model
    C_old(1) = const_bc;
    temp = AR*C_old;
    temp(end) = temp(end) + flux_bc * coef;
    C_new = AL\ temp;
    C_new = (C_new>0).*C_new;
end

function [C_new, rates] = sediments_chemical_reactions_module(sediment_params, C0,dt, method)
    if method == 0
        [C_new, rates] = rk4(sediment_params, C0,dt);
    elseif method == 1
        [C_new, rates] = butcher5(sediment_params, C0,dt);
    end
    C_new = (C_new>0).*C_new;
end

%% rk4: Runge-Kutta 4th order integration
function [C_new, rates] = rk4(sediment_params, C0, dt)
    % ts - how many time steps during 1 day

    [dcdt_1, r_1] = sediment_rates(sediment_params, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = sediment_rates(sediment_params, C0+0.5.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = sediment_rates(sediment_params, C0+0.5.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = sediment_rates(sediment_params, C0+k_3, dt);
    k_4 = dt.*dcdt_4;
    C_new = C0 + (k_1+2.*k_2+2.*k_3+k_4)/6;
    C0 = C_new;

    if sediment_params.rate_estimator_switch
      % average rate
      rate_number = fieldnames(r_1);
      for fld_idx = 1:numel(rate_number)
        rates.(rate_number{fld_idx}) = (r_1.(rate_number{fld_idx}) + 2*r_2.(rate_number{fld_idx}) + 2*r_3.(rate_number{fld_idx}) + r_4.(rate_number{fld_idx}))/6;
      end
    else
      rates = false;
    end
end

%% butcher5: Butcher's Fifth-Order Runge-Kutta
function [C_new, rates] = butcher5(sediment_params, C0,dt)

    [dcdt_1, r_1] = sediment_rates(sediment_params, C0, dt);
    k_1 = dt.*dcdt_1;
    [dcdt_2, r_2] = sediment_rates(sediment_params, C0 + 1/4.*k_1, dt);
    k_2 = dt.*dcdt_2;
    [dcdt_3, r_3] = sediment_rates(sediment_params, C0 + 1/8.*k_1 + 1/8.*k_2, dt);
    k_3 = dt.*dcdt_3;
    [dcdt_4, r_4] = sediment_rates(sediment_params, C0 - 1/2.*k_2 + k_3, dt);
    k_4 = dt.*dcdt_4;
    [dcdt_5, r_5] = sediment_rates(sediment_params, C0 + 3/16.*k_1 + 9/16.*k_4, dt);
    k_5 = dt.*dcdt_5;
    [dcdt_6, r_6] = sediment_rates(sediment_params, C0 - 3/7.*k_1 + 2/7.*k_2 + 12/7.*k_3 - 12/7.*k_4 + 8/7.*k_5, dt);
    k_6 = dt.*dcdt_6;
    C_new = C0 + (7.*k_1 + 32.*k_3 + 12.*k_4 + 32.*k_5 + 7.*k_6)/90;
    C0 = C_new;

    % average rate
    if sediment_params.rate_estimator_switch
      rate_number = fieldnames(r_1);
      for fld_idx = 1:numel(rate_number)
        rates.(rate_number{fld_idx}) = (7*r_1.(rate_number{fld_idx}) + 32*r_3.(rate_number{fld_idx}) + 12*r_4.(rate_number{fld_idx}) + 32*r_5.(rate_number{fld_idx}) + 7*r_6.(rate_number{fld_idx}))/90;
      end
    else
      rates = false;
    end
end

%% top_sediment_rate_to_flux: returns the flux of species at SWI converted to th units used in WC [ mg m-2 d-1 ].
function [flux] = top_sediment_rate_to_flux(R, dx)
    % TODO: Check units here!!
    % R - rate [umol cm-3 y-1]
    % dx - mesh size [cm];
    flux = integrate_over_depth(R, dx);
end

function [int_rate] = integrate_over_depth(R, dx)
  %% integrate_over_depth: integrates the rates of reaction over the depth and return the average values for the current day
  % int_rate  [umol cm-2 year -1 ]
  % R - rate of interest
  % dx - the mesh size
  % m - number of time step in 1 day
  int_rate = sum(daily_average(R),1)*dx;
end

function [int_rate] = integrate_over_depth_2(R, z)
  %% integrate_over_depth_2: integrates the rates of reaction over the depth and return the average values for the current day using trapezoidal rule
  % int_rate  [umol cm-2 year -1 ]
  % R - rate of interest
  % z - the depth
  if size(R,1) == 1
    int_rate = 0;
  else
    int_rate = trapz(z,R);
  end
end

%% daily_average: returns the average rate during 1 run (usually it is during 1 day if run with MyLake)
function [averaged] = daily_average(R)
  % R - rate of interest
  averaged = sum(R,2)/size(R,2);
end

function [flux] = convert_flux_umol_per_cm2_y_to_mg_per_m2_d(fx, M_C)
  flux = fx * M_C * 10^4 / 365 / 10^6; % [umol/cm^2/y] -> [mg/m2/d] NOTE: devision by 10^6 because molar in mg/mol (eg, 31000 for P)
end

function [flux] = top_sediment_diffusion_flux(C, D, dx, phi)
  % calculates flux of the particular dissolved specie through the top boundary of the sediment
  % in [ mg m-2 d-1 ] units
  % C(1) - BC of dissolved species
  % C - concentration
  % D - diffusion coefficient
  % M_C - molar mass in [ mg mol-1]
  % phi - porosity (no porosity because C is the concentration in pores (not bulk))

  % fourth-order
  flux = D * (-25 * phi(1)*C(1) + 48 * phi(2)*C(2) - 36 * phi(3)*C(3) + 16 * phi(4)*C(4) - 3 * phi(5)*C(5)) / dx / 12;  %  [umol/cm^2/y]

  % third order
  % flux = D * (-11 * C(1) + 18 * C(2) - 9 * C(3) + 2 * C(4)) / dx / 6;  %  [umol/cm^2/y]
  % flux = 0;  %  [umol/cm^2/y]

  % second order
  % flux = D * (-3 * C(1) + 4 * C(2) - C(3)) / dx / 2;  %  [umol/cm^2/y]

  % first order
  % flux = D * (C(3) - C(1)) / 2 / dx;  %  [umol/cm^2/y]
end

function [flux] = top_sediment_advection_flux(C, w, phi)
  % calculates flux of the particular dissolved specie through the top boundary of the sediment
  % in [ mg m-2 d-1 ] units
  % C - concentration
  % phi - porosity (no porosity because C is the concentration in pores (not bulk))
  % w - advection velocity
  % minus sign because of velocity signs

  % fourth-order
  flux = - phi(1) * w * C(1) ;  %  [umol/cm^2/y]
end


function bioR = bioirrigation(C, alfax, phi)
  % bioirrigation rate is the "artificial" function represents the bioirrigation by organisms in the sediment (worms etc) implemented according to Boudreau, B.P., 1999.
  % Co - bc wc-sediment value of current species
  % C - concentration profile of current species
  % phi - porosity
  Co = C(1);
  bioR = alfax .* (C - Co) .* phi;
  % NOTE:Disabled?
  % bioR = 0;
end

function [H3O] = pH_module(algorithm, H3O, CO2g, HCO3, CO2, CO3, NH3, NH4, HS, H2S, Fe2, Ca2, NO3, SO4, PO4, FeS, FeS2, FeOH3, FeOOH, Ca3PO42, PO4adsa, PO4adsb, Temperature, sediment_params)
  %% pH_module: pH equilibrium function
  % 0. No pH module
  % 1. Phreeqc
  % 2. New algorithm by Markelov (under test)
  % 3. Phreeqc Py
  phi = sediment_params.phi;

    if algorithm == 1 %
        Kw=10^(-14);
        OH = Kw./H3O/1e-3;
        T = Temperature*ones(size(H3O));
        P = sediment_params.pressure*ones(size(H3O));
        in =[H3O HCO3 CO2 CO3 NH3 NH4 HS H2S OH CO2g Fe2 Ca2 NO3 SO4 PO4 FeS FeS2 FeOH3 FeOOH Ca3PO42 PO4adsa.*(1-phi)./phi PO4adsb.*(1-phi)./phi];
        [pH_est] = pH_phreeqc(size(H3O,1),in);
        H3O = 10.^(-pH_est')*1e3;

    elseif algorithm == 2
      for i=size(H3O,1):-1:1
        sediment_params.aq_system.carb_acid.conc = 1e-3*(CO2(i)+HCO3(i)+CO3(i));
        sediment_params.aq_system.amonia.conc = 1e-3*(NH4(i)+NH3(i));
        sediment_params.aq_system.sulf.conc = 1e-3*(H2S(i)+HS(i));
        sediment_params.aq_system.ca.conc = 1e-3*(Ca2(i));
        sediment_params.aq_system.fe2.conc = 1e-3*(Fe2(i)); % 0.5 accounts for sorptio
        sediment_params.aq_system.no3.conc = 1e-3*(NO3(i));
        sediment_params.aq_system.so4.conc = 1e-3*(SO4(i));
        sediment_params.aq_system.p_acid.conc = 1e-3*(PO4(i));

        if i == size(H3O,1)
          pHs = linspace(0,14,1400)';
        else
          pHs = linspace(pHz(i+1)-0.1,pHz(i+1)+0.1,20)';
        end
        pHz(i) = new_pH_module(sediment_params.aq_system, pHs);
        H3O(i) = 10^(-pHz(i))*1e3;
      end
    end
end


function [dcdt, r] = sediment_rates(sediment_params, C, dt)
% parameters for water-column chemistry
% NOTE: the rates are the same as in sediments. Units are per "year" due to time step is in year units too;


    dcdt=zeros(size(C));

    if any(isnan(C))
      error('Breaking out of Sediments function: NaN values');
    end

    O2      = C(:,1) .* (C(:,1)>0) ;
    POP      = C(:,2) .* (C(:,2)>0) ;
    POC     = C(:,3) .* (C(:,3)>0) ;
    NO3     = C(:,4) .* (C(:,4)>0) ;
    FeOH3   = C(:,5) .* (C(:,5)>0) ;
    SO4     = C(:,6) .* (C(:,6)>0) ;
    NH4     = C(:,7) .* (C(:,7)>0) ;
    Fe2     = C(:,8) .* (C(:,8)>0) ;
    FeOOH   = C(:,9) .* (C(:,9)>0) ;
    H2S     = C(:,10) .* (C(:,10)>0) ;
    HS      = C(:,11) .* (C(:,11)>0) ;
    FeS     = C(:,12) .* (C(:,12)>0) ;
    S0      = C(:,13) .* (C(:,13)>0) ;
    PO4     = C(:,14) .* (C(:,14)>0) ;
    S8      = C(:,15) .* (C(:,15)>0) ;
    FeS2    = C(:,16) .* (C(:,16)>0) ;
    AlOH3   = C(:,17) .* (C(:,17)>0) ;
    PO4adsa = C(:,18) .* (C(:,18)>0) ;
    PO4adsb = C(:,19) .* (C(:,19)>0) ;
    Ca2     = C(:,20) .* (C(:,20)>0) ;
    Ca3PO42 = C(:,21) .* (C(:,21)>0) ;
    OMS     = C(:,22) .* (C(:,22)>0) ;
    CaCO3   = C(:,23) .* (C(:,23)>0) ;
    CO2     = C(:,24) .* (C(:,24)>0) ;
    CO3     = C(:,25) .* (C(:,25)>0) ;
    HCO3    = C(:,26) .* (C(:,26)>0) ;
    NH3     = C(:,27) .* (C(:,27)>0) ;
    CO2g   = C(:,28) .* (C(:,28)>0) ;
    DOP    = C(:,29) .* (C(:,29)>0) ;
    DOC    = C(:,30) .* (C(:,30)>0) ;
    Chl    = C(:,31) .* (C(:,31)>0) ;
    CH4aq    = C(:,32) .* (C(:,32)>0) ;
    CH4g    = C(:,33) .* (C(:,33)>0) ;
    FeCO3    = C(:,34) .* (C(:,34)>0) ;
    Fe3PO42    = C(:,35) .* (C(:,35)>0) ;
    PO4adsc    = C(:,36) .* (C(:,36)>0) ;

    H3O = sediment_params.H3O;
    HPO4 = sediment_params.HPO4;
    H2PO4 = sediment_params.H2PO4;

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
    CH4_solubility = sediment_params.CH4_solubility;
    CO2_solubility = sediment_params.CO2_solubility;
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
    alfax = sediment_params.alfax;
    phi    = sediment_params.phi;
    Kd_fe2 = sediment_params.Kd_fe2;

    Fe_dis = Fe2*Kd_fe2;


    tot_FeOH3 = PO4adsa + FeOH3;
    tot_FeOOH = PO4adsb + FeOOH;

    f_O2    = O2 ./  (Km_O2 + O2);
    f_NO3   = NO3 ./  (Km_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_FeOH3 = tot_FeOH3 ./  (Km_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_FeOOH = tot_FeOOH ./  (Km_FeOOH + tot_FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_SO4   = SO4 ./ (Km_SO4 + SO4 ) .* Kin_FeOOH ./ (Kin_FeOOH + FeOOH) .* Kin_FeOH3 ./ (Kin_FeOH3 + tot_FeOH3) .* Kin_NO3 ./ (Kin_NO3 + NO3) .* Kin_O2 ./ (Kin_O2 + O2);
    f_CH4 = 1 - f_O2 - f_NO3 - f_FeOH3 - f_FeOOH - f_SO4;
    f_CH4 = f_CH4.*(f_CH4>0);



    R1a = k_POP.*POP .* f_O2 * accel;
    R1b = k_POC.*POC .* f_O2 * accel;
    R1c = k_DOP.* DOP .* f_O2 * accel;
    R1d = k_DOC.*DOC .* f_O2 * accel;
    R1f = k_Chl.*Chl .* f_O2 * accel;

    R2a = k_POP.*POP .* f_NO3 * accel;
    R2b = k_POC.*POC .* f_NO3 * accel;
    R2c = k_DOP.*DOP .* f_NO3 * accel;
    R2d = k_DOC.*DOC .* f_NO3 * accel;
    R2f = k_Chl.*Chl .* f_NO3 * accel;


    part_PO4ads_tot_FeOH3 = PO4adsa ./ (tot_FeOH3+1e-16); % ratio of ads P to total Fe(III)
    R3a = k_POP.*POP .* f_FeOH3;
    R3b = k_POC .*POC .* f_FeOH3;
    R3c = k_DOP .*DOP .* f_FeOH3;
    R3d = k_DOC .*DOC .* f_FeOH3;
    R3f = k_Chl .*Chl .* f_FeOH3;
    R3a_Fe = (1 - part_PO4ads_tot_FeOH3) .* R3a;
    R3b_Fe = (1 - part_PO4ads_tot_FeOH3) .* R3b;
    R3c_Fe = (1 - part_PO4ads_tot_FeOH3) .* R3c;
    R3d_Fe = (1 - part_PO4ads_tot_FeOH3) .* R3d;
    R3f_Fe = (1 - part_PO4ads_tot_FeOH3) .* R3f;
    R3a_P = part_PO4ads_tot_FeOH3 .* R3a;
    R3b_P = part_PO4ads_tot_FeOH3 .* R3b;
    R3c_P = part_PO4ads_tot_FeOH3 .* R3c;
    R3d_P = part_PO4ads_tot_FeOH3 .* R3d;
    R3f_P = part_PO4ads_tot_FeOH3 .* R3f;



    part_PO4ads_tot_FeOOH = PO4adsb ./ (tot_FeOOH+1e-16); % ratio of ads P to total Fe(III)
    R4a = k_POP .*POP .* f_FeOOH;
    R4b = k_POC .*POC .* f_FeOOH;
    R4c = k_DOP .*DOP .* f_FeOOH;
    R4d = k_DOC .*DOC .* f_FeOOH;
    R4f = k_Chl .*Chl .* f_FeOOH;
    R4a_Fe = (1 - part_PO4ads_tot_FeOOH) .* R4a ;
    R4b_Fe = (1 - part_PO4ads_tot_FeOOH) .* R4b ;
    R4c_Fe = (1 - part_PO4ads_tot_FeOOH) .* R4c ;
    R4d_Fe = (1 - part_PO4ads_tot_FeOOH) .* R4d ;
    R4f_Fe = (1 - part_PO4ads_tot_FeOOH) .* R4f ;
    R4a_P = part_PO4ads_tot_FeOOH .* R4a;
    R4b_P = part_PO4ads_tot_FeOOH .* R4b;
    R4c_P = part_PO4ads_tot_FeOOH .* R4c;
    R4d_P = part_PO4ads_tot_FeOOH .* R4d;
    R4f_P = part_PO4ads_tot_FeOOH .* R4f;

    R5a = k_POP.*POP .* f_SO4  ;
    R5b = k_POC.*POC .* f_SO4 ;
    R5c = k_DOP.*DOP .* f_SO4 ;
    R5d = k_DOC.*DOC .* f_SO4 ;
    R5f = k_Chl.*Chl .* f_SO4 ;

    % NOTE: Disabled for the moment
    R6a = 0; %k_POP.*POP .* f_CH4;
    R6b = 0; %k_POC.*POC .* f_CH4;
    R6c = 0; %k_DOP.*DOP .* f_CH4;
    R6d = 0; %k_DOC.*DOC .* f_CH4;
    R6f = 0; %k_Chl.*Chl .* f_CH4;



    Ra = R1a+R2a+R3a+R4a+R5a+R6a;
    Rb = R1b+R2b+R3b+R4b+R5b+R6b;
    Rc = R1c+R2c+R3c+R4c+R5c+R6c;
    Rd = R1d+R2d+R3d+R4d+R5d+R6d;
    Rf = R1f+R2f+R3f+R4f+R5f+R6f;

    R1 = R1a+R1b+R1c+R1d+R1f;
    R2 = R2a+R2b+R2c+R2d+R2f;
    R3 = R3a+R3b+R3c+R3d+R3f;
    R4 = R4a+R4b+R4c+R4d+R4f;
    R5 = R5a+R5b+R5c+R5d+R5f;
    R6 = R6a+R6b+R6c+R6d+R6f;


    Sum_H2S = H2S + HS;;
    R11 = k_tsox * O2 .* Sum_H2S;
    R12 = k_tS_Fe * FeOH3 .*  Sum_H2S;

    R13 = k_Feox .* Fe_dis .* O2;
    % NOTE: Due to the reaction is too fast and could cause overshooting:
    % we need to make this check if R*dt > Conc of source:
    % if R*dt > Conc then R13 = C/dt
    % if R*dt < Conc then R13 = R13
    % R13 = (R13.*dt < Fe_dis/50).*R13 + (R13.*dt > Fe_dis/50).* R13 ./ 1000;
    % R13 = (R13.*dt < Fe_dis).*R13 + (R13.*dt > Fe_dis).* Fe_dis ./ (dt) * 0.5;
    % R13 = (R13.*dt < O2).*R13 + (R13.*dt > O2).* O2 ./ (dt) * 0.5;

    % R14 = k_amox * O2 ./ (Km_oxao + O2) .* (NH4 ./ (Km_amao + NH4)); % NOTE: Doesnt work - Highly unstable.
    R14 = k_amox  .* NH4 .* O2;
    % R14 = (R14.*dt < NH4).*R14 + (R14.*dt > NH4).* NH4 ./ (dt) * 0.5;
    % R14 = (R14.*dt < O2).*R14 + (R14.*dt > O2).* O2 ./ (dt) * 0.5;

    R15 = k_ch4_o2 .* CH4aq .* O2;
    R16 = k_ch4_so4 .* CH4aq .* SO4;


    CH4_over_sat = CH4aq - CH4_solubility;
    R17 = k_ch4_dis .* CH4_over_sat .* (CH4_over_sat > 0);

    CO2_over_sat = CO2 - CO2_solubility;

    % Minerals and solids:

    R21a = k_oms * Sum_H2S .* POP;
    R21b = k_oms * Sum_H2S .* POC;
    R21c = k_oms * Sum_H2S .* DOP;
    R21d = k_oms * Sum_H2S .* DOC;
    R21f = k_oms * Sum_H2S .* Chl;

    % NOTE: Could cause instability. These rates are too high when pH > 7
    R22a = k_Spre * S0;
    R22b = k_Sdis .* S8;

    R23 = k_fes2ox .* FeS2 .* O2;
    R24 = k_fespre .* FeS .* S0;
    R25 = k_fesox * O2 .* FeS;
    R26 = k_fes2pre .* FeS .* Sum_H2S;

    Sat_FeS = Fe_dis*1e-3 .* Sum_H2S*1e-3 ./ (H3O*1e-3+1e-16).^2 ./ K_FeS;
    R27a = k_fe_pre .* (Sat_FeS-1) .* (Sat_FeS > 1);
    R27b  =  k_fe_dis .* FeS .* (1-Sat_FeS) .* (Sat_FeS < 1); %


    Sat_CC = Ca2*1e-3 .* CO3*1e-3 / K_CC;
    R28a = k_CCpre .* (Sat_CC-1) .* (Sat_CC > 1);
    R28b  =  k_CCdis .* CaCO3 .* (1-Sat_CC) .* (Sat_CC < 1); %

    Sat_FC = Fe_dis*1e-3 .* CO3*1e-3 / K_FC;
    R29a = k_FCpre .* (Sat_FC-1) .* (Sat_FC > 1);
    R29b  =  k_FCdis .* FeCO3 .* (1-Sat_FC) .* (Sat_FC < 1); %



    % Phosphorus

    R31a = k_pdesorb_a * FeOH3 .* PO4;
    R31b = 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx2*R3c_P + Cx3*R3d_P + Cx1*R3f_P); % f_pfe .* (4 * R3 + 2 * R12);
    % R31b = (R31b.*dt < PO4adsa).*R31b + (R31b.*dt > PO4adsa).* PO4adsa ./ (dt) * 0.5;

    R32a = k_pdesorb_b * FeOOH .* PO4;
    R32b = 4 * (Cx2*R4a_P + Cx3*R4b_P + Cx2*R4c_P + Cx3*R4d_P + Cx1*R4f_P);
    % R32b = (R32b.*dt < PO4adsb).*R32b + (R32b.*dt > PO4adsb).* PO4adsb ./ (dt) * 0.5;


    Sat_viv = (HPO4*1e-3).^2.*(Fe_dis*1e-3).^3./(H3O*1e-3+1e-16).^2./K_viv;
    R33a = k_viv_pre .* (Sat_viv-1) .* (Sat_viv > 1);
    R33b = k_viv_dis .* Fe3PO42 .* (1 - Sat_viv) .* (Sat_viv < 1);


    % Sat_apa = (HPO4*1e-3).^2.*(Ca2*1e-3).^3./(H3O*1e-3+1e-16).^2./K_apa;
    Sat_apa = (H2PO4*1e-3).^2.*(Ca2*1e-3).^3./(H3O*1e-3+1e-16).^2./K_apa;
    R34a = k_apa_pre .* (Sat_apa-1) .* (Sat_apa > 1);
    R34b = k_apa_dis .* Ca3PO42 .* (1 - Sat_apa) .* (Sat_apa < 1);

    % R35 disabled now
    R35a = k_pdesorb_c .* PO4 .* AlOH3; %
    R35b = 0;



    % saving all rates
    if sediment_params.rate_estimator_switch
      r.f_O2 = f_O2; r.f_NO3 = f_NO3; r.f_FeOH3 = f_FeOH3; r.f_FeOOH = f_FeOOH; r.f_SO4 = f_SO4; r.f_CH4 = f_CH4; r.f_CH4 = f_CH4; r.R1a = R1a; r.R1b = R1b; r.R1c = R1c; r.R1d = R1d; r.R1f = R1f; r.R2a = R2a; r.R2b = R2b; r.R2c = R2c; r.R2d = R2d; r.R2f = R2f; r.R3a = R3a; r.R3b = R3b; r.R3c = R3c; r.R3d = R3d; r.R3f = R3f; r.R3a_Fe = R3a_Fe; r.R3b_Fe = R3b_Fe; r.R3c_Fe = R3c_Fe; r.R3d_Fe = R3d_Fe; r.R3f_Fe = R3f_Fe; r.R3a_P = R3a_P; r.R3b_P = R3b_P; r.R3c_P = R3c_P; r.R3d_P = R3d_P; r.R3f_P = R3f_P; r.R4a = R4a; r.R4b = R4b; r.R4c = R4c; r.R4d = R4d; r.R4f = R4f; r.R4a_Fe = R4a_Fe; r.R4b_Fe = R4b_Fe; r.R4c_Fe = R4c_Fe; r.R4d_Fe = R4d_Fe; r.R4f_Fe = R4f_Fe; r.R4a_P = R4a_P; r.R4b_P = R4b_P; r.R4c_P = R4c_P; r.R4d_P = R4d_P; r.R4f_P = R4f_P; r.R5a = R5a; r.R5b = R5b; r.R5c = R5c; r.R5d = R5d; r.R5f = R5f; r.R6a = R6a; r.R6b = R6b; r.R6c = R6c; r.R6d = R6d; r.R6f = R6f; r.Ra = Ra; r.Rb = Rb; r.Rc = Rc; r.Rd = Rd; r.Rf = Rf; r.R1 = R1; r.R2 = R2; r.R3 = R3; r.R4 = R4; r.R5 = R5; r.R6 = R6; r.Sum_H2S = Sum_H2S; r.R11 = R11; r.R12 = R12; r.R13 = R13; r.R14 = R14; r.R15 = R15; r.R16 = R16; r.R17 = R17; r.R21a = R21a; r.R21b = R21b; r.R21c = R21c; r.R21d = R21d; r.R21f = R21f; r.R22a = R22a; r.R22b = R22b; r.R23 = R23; r.R24 = R24; r.R25 = R25; r.R26 = R26; r.R27a = R27a; r.R27b = R27b; r.R28a = R28a; r.R28b = R28b; r.R29a = R29a; r.R29b = R29b; r.R31a = R31a; r.R31b = R31b; r.R32a = R32a; r.R32b = R32b; r.R33a = R33a; r.R33b = R33b; r.R34a = R34a; r.R34b = R34b; r.R35a = R35a; r.R35b = R35b;
    else
      r.R1a = 0;
    end

      % save only specific rates:
    % r.R28a = R28a;
    % r.R28b = R28b;
    % r.R31a = R31a;
    % r.R31b = R31b;
    % r.R32a = R32a;
    % r.R32b = R32b;
    % r.Sat_viv = Sat_viv;
    % r.R33a = R33a;
    % r.R33b = R33b;
    % r.Sat_apa = Sat_apa;
    % r.R34a = R34a;
    % r.R34b = R34b;
    % r.R35a = R35a;
    % r.R35b = R35b;


    % for stoichiometry check:
    % Canavan, R. W., Slomp, C. P., Jourabchi, P., Van Cappellen, P., Laverman, A. M., & van den Berg, G. A. (2006). Organic matter mineralization in sediment of a coastal freshwater lake and response to salinization. Geochimica Et Cosmochimica Acta, 70(11), 2836???2855. http://doi.org/10.1016/j.gca.2006.03.012


    % NOTE: R.*F the rate is written per solid;
    % R./F the rate is written per aqueous;
    F = (1-phi) ./ phi;

    dcdt(:,1)  = - bioirrigation(O2, alfax, phi) +  -0.25 * R13  - R15 - 2 * R14  - 2* R11 - (Cx2*R1a + Cx3*R1b+Cx1*R1f) .* F - (Cx2*R1c + Cx3*R1d) - 3 * R25.*F  - 5*R23; % O2 (aq)
    dcdt(:,2)  = -Ra - Cx1*R21a./F; % POP (solid)
    dcdt(:,3)  = -Rb - Cx1*R21b./F; % POC (solid)
    dcdt(:,4)  = - bioirrigation(NO3, alfax, phi) +  - 0.8*(Cx2*R2a+Cx2*R2b+Cx1*R2f) .* F - 0.8*(Cx2*R2c+Cx2*R2d)+ R14; % NO3(aq)
    dcdt(:,5)  = -4 * (Cx2*R3a_Fe+Cx3*R3b_Fe+Cx1*R3f_Fe) - 4*(Cx2*R3c_Fe + Cx3*R3d_Fe)./F - 2*R12 + R13./ F - R31a; % FeOH3(solid)
    dcdt(:,6)  = - bioirrigation(SO4, alfax, phi) +  - 0.5*(Cx2*R5a + Cx3*R5b+ Cx1*R5f) .* F -0.5*(Cx2*R5c + Cx3*R5d)+ R11 - R16 + 2*R23; % SO4(aq)
    dcdt(:,7)  = - bioirrigation(NH4, alfax, phi) +  (Ny2 * Ra + Ny3 * Rb+ Ny1 * Rf) .* F + (Ny2 * Rc + Ny3 * Rd) - R14; % NH4(aq)
    dcdt(:,8)  = - bioirrigation(Fe_dis, alfax, phi) +  4*(Cx2*R3a+Cx3*R3b+Cx1*R3f) .* F + 4* (Cx2*R3c + Cx3*R3d) + 4*(Cx2*R4a + Cx3*R4b+ Cx1*R4f) .* F + 4 * (Cx2*R4c + Cx3*R4d) + 2*R12.*F - R13 + R27b.*F - R27a - R29a + R29b.*F -3*R33a + 3*R33b.*F; % Fe2(aq)
    dcdt(:,9)  = -4 * (Cx2*R4a_Fe + Cx3*R4b_Fe + Cx1*R4f_Fe) - 4*(Cx2*R4c_Fe + Cx3*R4d_Fe)./F + R25 - R32a + R23./F; % FeOOH(solid)
    dcdt(:,10) = - bioirrigation(H2S, alfax, phi) - H2S./(Sum_H2S+1e-16).* (R11 + R12.*F + R27b.*F + R27a + R26 + Cx2*R21a + Cx3*R21b + Cx2*R21c + Cx3*R21d + Cx1*R21f); % H2S(aq)
    dcdt(:,11) = - bioirrigation(HS, alfax, phi) +  0.5*(Cx2*R5a + Cx3*R5b + Cx1*R5f) .* F + 0.5 * (Cx2*R5c + Cx3*R5d)  + R16 - HS./(Sum_H2S+1e-16).* (R11 + R12.*F + R27b.*F + R27a + R26 + Cx2*R21a + Cx3*R21b + Cx2*R21c + Cx3*R21d + Cx1*R21f) ; % HS(aq)
    dcdt(:,12) =  - R24 - 4*R25 -R26./F + R27a./F - R27b ; % FeS(solid)
    dcdt(:,13) = - R24.*F - R22a + R12.*F + R22b.*F; % S0(aq)
    dcdt(:,14) = - bioirrigation(PO4, alfax, phi) +  (Pz2 * Ra + Pz3 * Rb + Pz1 * Rf) .* F + (Pz2 * Rc + Pz3 * Rd) + 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx1*R3f_P) .*F + 4 * (Cx2*R3c_P + Cx3*R3d_P)  + 4 * (Cx2*R4a_P + Cx3*R4b_P + Cx1*R4f_P).*F  + 4 * (Cx2*R4c_P + Cx3*R4d_P) - R31a.*F - R32a.*F - R35a.*F + R35b.*F - 2 * R34a + 2 * R34b.*F -2*R33a + 2*R33b.*F; % PO4(aq)
    dcdt(:,15) = 4*R25 - R22b + R22a./F; % S8(solid)
    dcdt(:,16) = + R24 + R26./F - R23./F; % FeS2(solid)
    dcdt(:,17) = -R35a; % AlOH3(s)
    dcdt(:,18) = R31a - 4 * (Cx2*R3a_P + Cx3*R3b_P + Cx1*R3f_P) - 4* (Cx2*R3c_P + Cx3*R3d_P)./F; % PO4adsa(s)
    dcdt(:,19) = R32a - 4 * (Cx2*R4a_P + Cx3*R4b_P + Cx1*R4f_P) - 4 * (Cx2*R4c_P + Cx3*R4d_P)./F; % PO4adsb(s)
    dcdt(:,20) = - bioirrigation(Ca2, alfax, phi) -3*R34a + 3*R34b.*F - R28a + R28b.*F; % Ca2(aq)
    dcdt(:,21) = R34a./F - R34b; % Ca3PO42(s)
    dcdt(:,22) = (Cx2*R21a + Cx3*R21b + Cx2*R21c + Cx3*R21d + Cx1*R21f)./F; % OMS(s)
    dcdt(:,23) = R28a.*F - R28b; % CaCO3(s)
    dcdt(:,24) = - bioirrigation(CO2, alfax, phi) + (R15 +  ((Cx2 - Ny2 + 2*Pz2)*R1a + (Cx3 - Ny3 + 2*Pz3)*R1b + (Cx1 - Ny1 + 2*Pz1)*R1f + (0.2*Cx2 - Ny2 + 2*Pz2)*R2a +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2b +  (0.2*Cx1 - Ny1 + 2*Pz1)*R2f  + (0.5 * Cx2 - Ny2 - 2*Pz2)*R6a + (0.5 * Cx3 - Ny3 - 2*Pz3)*R6b + (0.5 * Cx1 - Ny1 - 2*Pz1)*R6f) .* F  +  (Cx2 - Ny2 + 2*Pz2)*R1c + (Cx3 - Ny3 + 2*Pz3)*R1d + (0.2*Cx2 - Ny2 + 2*Pz2)*R2c +  (0.2*Cx3 - Ny3 + 2*Pz3)*R2d - (7*Cx2 + Ny2 + 2*Pz2)*(R3c+R4c) + (0.5*Cx2 - Ny2 + 2*Pz2)*R6c + (0.5*Cx3 - Ny3 + 2*Pz3)*R6d).*(CO2_over_sat<0) + (- (7*Cx2 + Ny2 + 2*Pz2)*(R3a+R4a) - (7*Cx3 + Ny3 + 2*Pz3)*(R3b+R4b) - (7*Cx1 + Ny1 + 2*Pz1)*(R3f+R4f)  - (Ny2 - 2*Pz2)*R5a - (Ny3 - 2*Pz3)*R5b - (Ny1 - 2*Pz1)*R5f) .* F  - (7*Cx3 + Ny3 + 2*Pz3)*(R3d+R4d)  - (Ny2 - 2*Pz2)*R5c - (Ny3 - 2*Pz3)*R5d - R16 + 2*R13 + 2*R14;  % CO2 (aq) NOTE: we need to add gas (Rename CO2g to gas)
    dcdt(:,25) = - bioirrigation(CO3, alfax, phi) - R28a + R28b.*F - R29a + R29b.*F; % CO3(aq)
    dcdt(:,26) = - bioirrigation(HCO3, alfax, phi)+  ( - (Ny2 - 2*Pz2)*R1a - (Ny3 - 2*Pz3)*R1b - (Ny1 - 2*Pz1)*R1f +  (0.8*Cx2 + Ny2 - 2*Pz2)*R2a + (0.8*Cx3 + Ny3 - 2*Pz3)*R2b + (0.8*Cx1 + Ny1 - 2*Pz1)*R2f  + (8*Cx2+Ny2-2*Pz2)*(R3a + R4a) +(8*Cx3+Ny3-2*Pz3)*(R3b + R4b) +(8*Cx1+Ny1-2*Pz1)*(R3f + R4f)  + (Cx2+Ny2-2*Pz2)*R5a + (1*Cx3+Ny3-2*Pz3)*R5b + (1*Cx1+Ny1-2*Pz1)*R5f + (Ny2-2*Pz2)*R6a + (Ny3-2*Pz3)*R6b + (Ny1-2*Pz1)*R6f) .* F - (Ny2 - 2*Pz2)*R1c - (Ny3 - 2*Pz3)*R1d + (0.8*Cx2 + Ny2 - 2*Pz2)*R2c + (0.8*Cx3 + Ny3 - 2*Pz3)*R2d + (8*Cx2+Ny2-2*Pz2)*(R3c + R4c) +(8*Cx3+Ny3-2*Pz3)*(R3d + R4d) + (Cx2+Ny2-2*Pz2)*R5c + (1*Cx3+Ny3-2*Pz3)*R5d + (Ny2-2*Pz2)*R6c + (Ny3-2*Pz3)*R6d  -  2*R13 - 2*R14 + 2*R16; % HCO3(aq)
    dcdt(:,27) = - bioirrigation(NH3, alfax, phi) ; % NH3(aq)
    dcdt(:,28) = - bioirrigation(CO2g, alfax, phi) ; % CO2g
    dcdt(:,29) = - bioirrigation(DOP, alfax, phi)  - Rc - Cx1*R21c; % DOP (aq)
    dcdt(:,30) = - bioirrigation(DOC, alfax, phi)  - Rd - Cx1*R21d; % DOC (aq)
    dcdt(:,31) = -Rf - Cx1*R21f./F; % Chl (s)
    dcdt(:,32) = - bioirrigation(CH4aq, alfax, phi) + 0.5*(Cx2*R6a+Cx2*R6b+Cx1*R6f).* F .* (CH4_over_sat < 0) + 0.5*(Cx2*R6c + Cx3*R6d).* (CH4_over_sat < 0) - R15 - R16 - R17;  % CH4aq
    dcdt(:,33) = - bioirrigation(CH4g, alfax, phi) + 0.5*(Cx2*R6a+Cx2*R6b+Cx1*R6f).* F .* (CH4_over_sat > 0) + 0.5*(Cx2*R6c + Cx3*R6d).* (CH4_over_sat > 0) + R17;  % CH4g
    dcdt(:,34) = R29a./F - R29b;  % FeCO3
    dcdt(:,35) = R33a./F - R33b;  % Fe3PO42
    dcdt(:,36) = R35a - R35b;  % P=Al(OH)3
end


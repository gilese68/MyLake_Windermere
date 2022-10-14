classdef output_matrices
    %OUTPUT_MATRICES Initialise and updates output matrices
    %   Output matrices are initialise at the beginning of the simulation
    %   and updated at each time step
    
    properties
        tanks
    end
    
    methods
        function obj = output_matrices(tt, Nz_tanks, Nt, N_algae)
            %OUTPUT_MATRICES Construct an instance of this class
            % Allocate and initialise output data matrices

            obj.tanks = struct();
            for t=1:Nt
                Nz = Nz_tanks(t);
                obj.tanks(t).Qst = zeros(3,length(tt));
                obj.tanks(t).Kzt = zeros(Nz,length(tt));
                obj.tanks(t).rhozt = zeros(Nz,length(tt));
                obj.tanks(t).Tzt = zeros(Nz,length(tt));
                for k=1:N_algae
                    obj.tanks(t).Chlzt{k} = zeros(Nz,length(tt));
                end
                obj.tanks(t).POCzt = zeros(Nz,length(tt));
                obj.tanks(t).Pzt = zeros(Nz,length(tt));
                obj.tanks(t).PPzt = zeros(Nz,length(tt));
                obj.tanks(t).DOPzt = zeros(Nz,length(tt));
                obj.tanks(t).DOCzt = zeros(Nz,length(tt));
                obj.tanks(t).DICzt = zeros(Nz,length(tt));
                obj.tanks(t).CO2zt = zeros(Nz,length(tt));
                obj.tanks(t).HCO3zt = zeros(Nz,length(tt));
                obj.tanks(t).CO3zt = zeros(Nz,length(tt));
                obj.tanks(t).O2zt = zeros(Nz,length(tt));

                obj.tanks(t).NO3zt  = zeros(Nz,length(tt));
                obj.tanks(t).NH4zt  = zeros(Nz,length(tt));
                obj.tanks(t).SO4zt  = zeros(Nz,length(tt));
                obj.tanks(t).HSzt  = zeros(Nz,length(tt));
                obj.tanks(t).H2Szt  = zeros(Nz,length(tt));
                obj.tanks(t).Fe2zt  = zeros(Nz,length(tt));
                obj.tanks(t).Ca2zt  = zeros(Nz,length(tt));
                obj.tanks(t).pHzt  = zeros(Nz,length(tt));
                obj.tanks(t).CH4aqzt  = zeros(Nz,length(tt));
                obj.tanks(t).Fe3zt  = zeros(Nz,length(tt));
                obj.tanks(t).Al3zt  = zeros(Nz,length(tt));
                obj.tanks(t).FeSzt  = zeros(Nz,length(tt));
                obj.tanks(t).CaCO3zt  = zeros(Nz,length(tt));
                obj.tanks(t).CH4gzt  = zeros(Nz,length(tt));
                obj.tanks(t).POPzt  = zeros(Nz,length(tt));
                obj.tanks(t).Sizt = zeros(Nz, length(tt));

                obj.tanks(t).H_sw_zt  = zeros(Nz,length(tt));
                obj.tanks(t).H_sw_zt_2  = zeros(Nz,length(tt));
                obj.tanks(t).PAR_zt  = zeros(Nz,length(tt));

                obj.tanks(t).O2_diffzt = zeros(Nz,length(tt));
                obj.tanks(t).O2_sat_relt = zeros(Nz,length(tt));
                obj.tanks(t).O2_sat_abst = zeros(Nz,length(tt));
                obj.tanks(t).Qzt_sed = zeros(Nz,length(tt));
                obj.tanks(t).lambdazt = zeros(Nz,length(tt));
                obj.tanks(t).P3zt_sed = zeros(Nz,length(tt),4); %3-D
                obj.tanks(t).P3zt_sed_sc = zeros(Nz,length(tt),3); %3-D
                obj.tanks(t).His = zeros(8,length(tt)); %NEW!!!
                obj.tanks(t).MixStat = zeros(23,length(tt));
                % Fokema
                obj.tanks(t).CDOMzt=zeros(Nz,length(tt));
                obj.tanks(t).DOCzt1=zeros(Nz,length(tt)); %Fokema-model subpool 1
                obj.tanks(t).DOCzt2=zeros(Nz,length(tt)); %Fokema-model subpool 2
                obj.tanks(t).DOCzt3=zeros(Nz,length(tt)); %Fokema-model subpool 3
                obj.tanks(t).DOC1tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 1
                obj.tanks(t).DOC2tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 2 fraction
                obj.tanks(t).DOC3tfrac=zeros(Nz,length(tt)); %Fokema-model subpool 3 fraction
                obj.tanks(t).Daily_BB1t=zeros(Nz,length(tt)); %Fokema-model subpool 1 daily bacterial decomposition
                obj.tanks(t).Daily_BB2t=zeros(Nz,length(tt)); %Fokema-model subpool 2 daily bacterial decomposition
                obj.tanks(t).Daily_BB3t=zeros(Nz,length(tt)); %Fokema-model subpool 3 daily bacterial decomposition
                obj.tanks(t).Daily_PBt=zeros(Nz,length(tt)); %Fokema-model daily photobleaching

                obj.tanks(t).surfaceflux = zeros(1,length(tt)); %CO2 surface flux
                obj.tanks(t).CO2_eqt = zeros(1,length(tt));     %CO2 equilibrium concentration
                obj.tanks(t).CO2_ppmt = zeros(1,length(tt));    %CO2 fraction in air
                obj.tanks(t).K0t = zeros(1,length(tt));         %CO2 solubility coefficient

                obj.tanks(t).O2fluxt = zeros(1,length(tt));     %oxygen surface flux
                obj.tanks(t).O2_eqt = zeros(1,length(tt));      %O2 equilibrium concentration
                obj.tanks(t).K0_O2t = zeros(1,length(tt));      %O2 solubility coefficient
                obj.tanks(t).dO2Chlt = zeros(Nz,length(tt));    %Oxygen change due to phytoplankton (mg m-3))
                obj.tanks(t).dO2BODt = zeros(Nz,length(tt));    %Oxygen consumption due to BOD (mg m-3))
                obj.tanks(t).dfloc_DOC =  zeros(Nz,length(tt));  % floculation rates
                obj.tanks(t).testi1t = zeros(Nz,length(tt));
                obj.tanks(t).testi2t = zeros(Nz,length(tt));
                obj.tanks(t).testi3t = zeros(Nz,length(tt));
                obj.tanks(t).lvlDzt = zeros(1,length(tt));
            end
        end
        
        function obj = update_output_matrices(obj, i, photobleaching, tanks, Nt, dt, N_algae)
            %UPDATE_OUTPUT_MATRICES at the end of each time step
            for t = 1:Nt
                obj.tanks(t).Qst(:,i) = tanks(t).results.Qs;
                obj.tanks(t).rhozt(:,i) = tanks(t).results.rho;
                obj.tanks(t).Kzt(:,i) = tanks(t).results.Kz;
                obj.tanks(t).Tzt(:,i) = tanks(t).profiles.Tz;
                for k = 1:N_algae
                    obj.tanks(t).Chlzt{k}(:,i) = tanks(t).profiles.Chlz{k};
                end
                obj.tanks(t).POCzt(:,i) = tanks(t).profiles.POCz;
                obj.tanks(t).Pzt(:,i) = tanks(t).profiles.Pz;
                obj.tanks(t).PPzt(:,i) = tanks(t).profiles.PPz;
                obj.tanks(t).DOPzt(:,i) = tanks(t).profiles.DOPz;
                obj.tanks(t).DOCzt(:,i) = tanks(t).profiles.DOCz;
                obj.tanks(t).DICzt(:,i) = tanks(t).profiles.DICz;
                obj.tanks(t).CO2zt(:,i) = tanks(t).profiles.CO2z;
                obj.tanks(t).HCO3zt(:,i) = tanks(t).profiles.HCO3z;
                obj.tanks(t).CO3zt(:,i) = tanks(t).profiles.CO3z;
                obj.tanks(t).O2zt(:,i) = tanks(t).profiles.O2z;

                obj.tanks(t).NO3zt(:,i) = tanks(t).profiles.NO3z;
                obj.tanks(t).NH4zt(:,i) = tanks(t).profiles.NH4z;
                obj.tanks(t).SO4zt(:,i) = tanks(t).profiles.SO4z;
                obj.tanks(t).HSzt(:,i) = tanks(t).profiles.HSz;
                obj.tanks(t).H2Szt(:,i) = tanks(t).profiles.H2Sz;
                obj.tanks(t).Fe2zt(:,i) = tanks(t).profiles.Fe2z;
                obj.tanks(t).Ca2zt(:,i) = tanks(t).profiles.Ca2z;
                obj.tanks(t).pHzt(:,i) = tanks(t).profiles.pHz;
                obj.tanks(t).CH4aqzt(:,i) = tanks(t).profiles.CH4aqz;
                obj.tanks(t).Fe3zt(:,i) = tanks(t).profiles.Fe3z;
                obj.tanks(t).Al3zt(:,i) = tanks(t).profiles.Al3z;
                obj.tanks(t).FeSzt(:,i) = tanks(t).profiles.FeSz;
                obj.tanks(t).CaCO3zt(:,i) = tanks(t).profiles.CaCO3z;
                obj.tanks(t).CH4gzt(:,i) = tanks(t).profiles.CH4gz;
                obj.tanks(t).POPzt(:,i) = tanks(t).profiles.POPz;
                obj.tanks(t).Sizt(:,i) = tanks(t).profiles.Siz;

                % TODO: [LC] generalise for N species
                obj.tanks(t).H_sw_zt(:,i) = tanks(t).results.H_sw_z(:, 1);
                obj.tanks(t).H_sw_zt_2(:,i) = tanks(t).results.H_sw_z(:, 2);
                obj.tanks(t).PAR_zt(:,i) = tanks(t).profiles.PAR_z;

                obj.tanks(t).lvlDzt(:,i) = tanks(t).results.lvlDz;

                % O2diffzt(:,i) = O2_diff;

                obj.tanks(t).O2_sat_relt(:,i) = tanks(t).results.O2_sat_rel;
                obj.tanks(t).O2_sat_abst(:,i) = tanks(t).results.O2_sat_abs;
                obj.tanks(t).BODzt = 0; %for compatibility with the other code

                %Fokema
                %CDOMzt(:,i)=CDOMz;
                if (photobleaching==1)
                    obj.tanks(t).DOCzt1(:,i) = tanks(t).results.DOCz1; %Fokema-model DOC subpool 1
                    obj.tanks(t).DOCzt2(:,i) = tanks(t).results.DOCz2; %Fokema-model DOC subpool 2
                    obj.tanks(t).DOCzt3(:,i) = tanks(t).results.DOCz3; %Fokema-model DOC subpool 3
                    obj.tanks(t).DOC1tfrac(:,i) = tanks(t).results.DOC1frac; %Fokema-model subpool 1 fraction
                    obj.tanks(t).DOC2tfrac(:,i) = tanks(t).results.DOC2frac; %Fokema-model subpool 2 fraction
                    obj.tanks(t).DOC3tfrac(:,i) = tanks(t).results.DOC3frac; %Fokema-model subpool 3 fraction
                    obj.tanks(t).Daily_BB1t(:,i) = tanks(t).results.Daily_BB1; %Fokema-model subpool 1 daily bacterial decomposition
                    obj.tanks(t).Daily_BB2t(:,i) = tanks(t).results.Daily_BB2; %Fokema-model subpool 2 daily bacterial decomposition
                    obj.tanks(t).Daily_BB3t(:,i) = tanks(t).results.Daily_BB3; %Fokema-model subpool 3 daily bacterial decomposition
                    obj.tanks(t).Daily_PBt(:,i) = tanks(t).results.Daily_PB; %Fokema-model daily photobleaching
                end

                obj.tanks(t).Qzt_sed(:,i) = tanks(t).results.Qz_sed./(60*60*24*dt); % (W m-2)
                obj.tanks(t).lambdazt(:,i) = tanks(t).results.lambdaz_wtot_avg;

                obj.tanks(t).surfaceflux(1,i) = tanks(t).results.surfflux; %Carbon dioxide surface flux

                obj.tanks(t).His(:,i) = tanks(t).results.His;

                obj.tanks(t).MixStat(:,i) = tanks(t).results.MixStat;
            end
        end
        
        function results = save_results(obj, Nt)
            %SAVE_RESULTS rearrange output matrices into results structure
            for t=1:Nt
                results.tanks(t).Qst = obj.tanks(t).Qst;
                results.tanks(t).K = obj.tanks(t).Kzt;
                results.tanks(t).rho = obj.tanks(t).rhozt;
                results.tanks(t).T = obj.tanks(t).Tzt;
                results.concentrations.tanks(t).P = obj.tanks(t).Pzt;
                results.concentrations.tanks(t).PP = obj.tanks(t).PPzt;
                results.concentrations.tanks(t).Chl = obj.tanks(t).Chlzt;
                results.concentrations.tanks(t).POP = obj.tanks(t).POPzt;
                results.concentrations.tanks(t).DOP = obj.tanks(t).DOPzt;
                results.concentrations.tanks(t).DOC = obj.tanks(t).DOCzt;
                results.concentrations.tanks(t).DIC = obj.tanks(t).DICzt;
                results.concentrations.tanks(t).CO2aq = obj.tanks(t).CO2zt;
                results.concentrations.tanks(t).HCO3 = obj.tanks(t).HCO3zt;
                results.concentrations.tanks(t).CO3 = obj.tanks(t).CO3zt;
                results.concentrations.tanks(t).O2 = obj.tanks(t).O2zt;
                results.concentrations.tanks(t).NO3 = obj.tanks(t).NO3zt;
                results.concentrations.tanks(t).NH4 = obj.tanks(t).NH4zt;
                results.concentrations.tanks(t).Fe3 = obj.tanks(t).Fe3zt;
                results.concentrations.tanks(t).Fe2 = obj.tanks(t).Fe2zt;
                results.concentrations.tanks(t).SO4 = obj.tanks(t).SO4zt;
                results.concentrations.tanks(t).HS = obj.tanks(t).HSzt;
                results.concentrations.tanks(t).H2S = obj.tanks(t).H2Szt;
                results.concentrations.tanks(t).Ca2 = obj.tanks(t).Ca2zt;
                results.concentrations.tanks(t).pH = obj.tanks(t).pHzt;
                results.concentrations.tanks(t).POC = obj.tanks(t).POCzt;
                results.concentrations.tanks(t).Al3 = obj.tanks(t).Al3zt;
                results.concentrations.tanks(t).FeS = obj.tanks(t).FeSzt;
                results.concentrations.tanks(t).CaCO3 = obj.tanks(t).CaCO3zt;
                results.concentrations.tanks(t).CH4aq = obj.tanks(t).CH4aqzt;
                results.concentrations.tanks(t).CH4g = obj.tanks(t).CH4gzt;
                results.concentrations.tanks(t).Si = obj.tanks(t).Sizt;
                results.tanks(t).O2_sat_relt = obj.tanks(t).O2_sat_relt;
                results.tanks(t).O2_sat_abst = obj.tanks(t).O2_sat_abst;
                results.tanks(t).Qzt_sed = obj.tanks(t).Qzt_sed;
                results.tanks(t).lambdazt = obj.tanks(t).lambdazt;
                results.tanks(t).P3zt_sed = obj.tanks(t).P3zt_sed;
                results.tanks(t).P3zt_sed_sc = obj.tanks(t).P3zt_sed_sc;
                results.tanks(t).His = obj.tanks(t).His;
                results.tanks(t).MixStat = obj.tanks(t).MixStat;
                results.tanks(t).H_sw = obj.tanks(t).H_sw_zt;
                results.tanks(t).H_sw_2 = obj.tanks(t).H_sw_zt_2;
                results.tanks(t).PAR_zt = obj.tanks(t).PAR_zt;
                results.tanks(t).lvlDzt = obj.tanks(t).lvlDzt;
                results.tanks(t).surfaceflux = obj.tanks(t).surfaceflux;
            end
        end
    end
end


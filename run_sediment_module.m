function [sediment_bioirrigation_fluxes, sediment_transport_fluxes,sediment_concentrations, sediment_additional_results, profiles] = ...
            run_sediment_module(profiles, mylake_temp_results, sediment_params, mylake_params, sediment_concentrations, sediment_matrix_templates, N_algae)
    for k=1:N_algae
        mylake_temp_results.Chlz{k} = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Chlz{k}, 30973.762);
    end
    mylake_temp_results.POPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.POPz, 30973.762);
    mylake_temp_results.DOCz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.DOCz, 12010.7);
    mylake_temp_results.DOPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.DOPz, 30973.762);
    mylake_temp_results.O2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.O2z, 31998.8);
    mylake_temp_results.Pz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Pz, 30973.762);
    mylake_temp_results.Fe2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Fe2z, 55845);
    mylake_temp_results.NO3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.NO3z, 62004);
    mylake_temp_results.NH4z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.NH4z, 18038);
    mylake_temp_results.SO4z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.SO4z, 96062);
    mylake_temp_results.HS = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.HSz, 33072);
    mylake_temp_results.Fe3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Fe3z, 106867.0);
    mylake_temp_results.Ca2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Ca2z, 40078.0);
    mylake_temp_results.Al3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.Al3z, 78003.6);
    mylake_temp_results.PPz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.PPz, 30973.762);
    mylake_temp_results.POCz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.POCz,  12010.7);
    mylake_temp_results.CH4aqz = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CH4aqz,  16042.5);
    mylake_temp_results.CO2z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CO2z,  44009.5);
    mylake_temp_results.HCO3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.HCO3z,  61016.8);
    mylake_temp_results.CO3z = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CO3z,  60008.9);
    mylake_temp_results.HSz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.HSz, 33072.9);
    mylake_temp_results.FeSz     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.FeSz, 87910.0);
    mylake_temp_results.CaCO3z     = convert_mg_per_qubic_m_to_umol_per_qubic_cm(profiles.CaCO3z, 100086.9);

    if sediment_params.effective_depth == -1
        effective_depth = TCz;
    else
        effective_depth = sediment_params.effective_depth;
    end
    sediment_bc = update_sediment(mylake_temp_results, mylake_params, sediment_params, effective_depth);

    % Running sediment module
    [sediment_bioirrigation_fluxes, sediment_transport_fluxes, sediment_concentrations, sediment_additional_results] = sediment_v2(...
        sediment_concentrations, sediment_params, sediment_matrix_templates, sediment_bc);

    mylake_prev_results.Tz = profiles.Tz;
    mylake_prev_results.O2z = profiles.O2z;
    mylake_prev_results.Pz = profiles.Pz;
    mylake_prev_results.Fe2z = profiles.Fe2z;
    mylake_prev_results.NO3z = profiles.NO3z;
    mylake_prev_results.NH4z = profiles.NH4z;
    mylake_prev_results.SO4z = profiles.SO4z;
    mylake_prev_results.HSz = profiles.HSz;
    mylake_prev_results.DOPz = profiles.DOPz;
    mylake_prev_results.DOCz = profiles.DOCz;
    mylake_prev_results.CH4aqz = profiles.CH4aqz;
    mylake_prev_results.CH4gz = profiles.CH4gz;
    mylake_prev_results.CO2z = profiles.CO2z;
    mylake_prev_results.HCO3z = profiles.HCO3z;
    mylake_prev_results.CO3z = profiles.CO3z;
    mylake_prev_results.Chlz = profiles.Chlz;
    mylake_prev_results.POPz = profiles.POPz;
    mylake_prev_results.POCz = profiles.POCz;
    mylake_prev_results.Fe3z = profiles.Fe3z;
    mylake_prev_results.FeSz = profiles.FeSz;
    mylake_prev_results.Al3z = profiles.Al3z;
    mylake_prev_results.PPz = profiles.PPz;
    mylake_prev_results.CaCO3z = profiles.CaCO3z;

    if any(isnan(profiles.O2z)) | any(isnan(profiles.Pz)) | any(isnan(profiles.Fe2z)) | any(isnan(profiles.NO3z)) | any(isnan(profiles.NH4z))
        error('NaN')
    end

    % Update WC:  [sediment] ----> [WC]
    [mylake_new_resutls] = update_wc(mylake_prev_results, sediment_concentrations, sediment_transport_fluxes, sediment_bioirrigation_fluxes, mylake_params, sediment_params, effective_depth);

    profiles.O2z = mylake_new_resutls.O2z;
    profiles.Pz = mylake_new_resutls.Pz;
    profiles.Fe2z = mylake_new_resutls.Fe2z;
    profiles.NO3z = mylake_new_resutls.NO3z;
    profiles.NH4z = mylake_new_resutls.NH4z;
    profiles.SO4z = mylake_new_resutls.SO4z;
    profiles.HSz = mylake_new_resutls.HSz;
    profiles.DOPz = mylake_new_resutls.DOPz;
    profiles.DOCz = mylake_new_resutls.DOCz;
    profiles.CH4aqz = mylake_new_resutls.CH4aqz;
    profiles.CH4gz = mylake_new_resutls.CH4gz;
    profiles.CO2z = mylake_new_resutls.CO2z;
    profiles.HCO3z = mylake_new_resutls.HCO3z;
    profiles.CO3z = mylake_new_resutls.CO3z;
    profiles.POPz = mylake_new_resutls.POPz;
    profiles.POCz = mylake_new_resutls.POCz;
    profiles.Fe3z = mylake_new_resutls.Fe3z;
    profiles.FeSz = mylake_new_resutls.FeSz;
    profiles.Al3z = mylake_new_resutls.Al3z;
    profiles.PPz = mylake_new_resutls.PPz;
    profiles.CaCO3z = mylake_new_resutls.CaCO3z;
    for k=1:N_algae
        profiles.Chlz{k} = mylake_new_resutls.Chlz{k};
    end

    profiles.DICz = profiles.CO2z + profiles.HCO3z + profiles.CO3z;

function C = convert_mg_per_qubic_m_to_umol_per_qubic_cm(C,M_C)
    C = C./M_C;
%end of function
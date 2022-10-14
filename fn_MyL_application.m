function [MyLake_results, Sediment_results] = fn_MyL_application(m_start,m_stop, K_sediments, K_lake, tanks_params, algal_params, name_of_scenario, run_ID, clim_ID, is_save_results, init_file)
% This is the main MyLake application configuration file. 

global sed_par_file lake_par_file Eevapor
Eevapor=0;

% writing sediments parameters file
%calibration_k_values = [(1:length(K_sediments))',cell2mat(K_sediments(:,1)) ];
calibration_k_values = zeros(length(K_sediments), 7);
calibration_k_values(:,1) = (1:length(K_sediments))';
for ii=1:size(K_sediments,1)
    calibration_k_values(ii, 2:end) = cell2mat(K_sediments(ii,1));
end

%% generates unique files for temporary parameter files
sed_par_file = tempname;
lake_par_file = tempname;
dlmwrite(sed_par_file, calibration_k_values,'delimiter','\t');

% Sets the Lake parameters read from load_params.m in a cell array
lake_params = zeros(length(K_lake), tanks_params{1}+1);
lake_params(:,1) = (1:length(K_lake))';
for ii=1:size(K_lake,1)
    for jj = 1:size(K_lake{1, 1},2) % jj = 1:size(K_lake{1, 1},2)
        lake_params(ii, jj+1) = K_lake{ii}{jj};
        
    end
end

% writing TEMP lake parameter file FOR use in HPC or similar parallel-type 
% system to reduces read-write traffic on network
fid=fopen(lake_par_file,'wt');
fprintf(fid,'\n\n');
% dlmwrite(lake_par_file, [[1:length(K_lake)]',data_lake{2},data_lake{3},data_lake{4},(1:length(K_lake))'],'delimiter','\t','-append'); % 1:length(K_lake) is the length of the parameter file.
dlmwrite(lake_par_file, lake_params,'delimiter','\t','-append'); 
fclose(fid);

% INPUT TIMESERIES FILE
inputfile=name_of_scenario;

%%%%%%%%%%%%%% Solve the model given the specified files %%%%%%%%%%%%%%%%%%
parafile=lake_par_file;
initfile={['IO/',init_file]};

[MyLake_results_basin1, sediment_results_basin1] ...
    = solvemodel_v3(m_start, m_stop, initfile, inputfile, parafile, tanks_params, algal_params);
MyLake_results.basin1 = MyLake_results_basin1; Sediment_results.basin1 = sediment_results_basin1;

if is_save_results
    %disp('Saving sediment and water-column profiles for basin 1');
    sediment_save_init_conc(Sediment_results.basin1, 1)
    MyLake_save_result_for_init_conc(MyLake_results.basin1, 1)
else
    %disp('Skipping saving the results and initial concentrations');
end

fclose('all'); % close all open files

end

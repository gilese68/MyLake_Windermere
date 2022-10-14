function [Chlorophyll, Temperature, Mixed_Depth] = Results_Processing(Sn,Sn2,N,Tk,ScenarioNum,Type)
%Results_Processing

if Type == 1
    for hh = 1:Sn+1 % Number of reductions
        for ii = 1:N % All the acceptable simulations
            
            % Load simulated
            formatSpec = '%d\\Wind_%d_%d';
            str = sprintf(formatSpec,ScenarioNum,hh,ii);
            load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\Windermere_Test\Wind_Ex_',str],'MyLake_results');
            
            for tnk_plot = 1:Tk % which tanks to plot
                
                % Chlorophyll
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla = MyLake_results.basin1.concentrations.tanks(tnk_plot).Chl; % All Chlorophyll (by species)
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum = sum(cat(3,Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{:}),3); % Total chlorophyll
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Top_Mean = nanmean(Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(1:6,:),1); % Surface to 5m mean total chlorophyll
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Diatoms = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{1} + Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{2};
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Green = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{3} + Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{4};
                Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Cyano = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{5} + Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla{6};
                
                % Temperature
                Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp = MyLake_results.basin1.tanks(tnk_plot).T; % All temperatures
                Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp_Top_Mean = nanmean(Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(1:6,:),1); % Surface to 5m mean temperature
                
                % Mixed depth
                Mixed_Depth.tanks(tnk_plot).Sim(hh).Mxd(ii,:) = MyLake_results.basin1.tanks(tnk_plot).MixStat(12,:);
                Mxd_lyr(ii,:) = Mixed_Depth.tanks(tnk_plot).Sim(hh).Mxd(ii,:);
                Mxd_lyr(isnan(Mxd_lyr)) = 17;
                Mxd_lyr = 17-Mxd_lyr;
                Mixed_Depth.tanks(tnk_plot).Sim(hh).metric(ii,:)= trapz(Mxd_lyr(ii,:));
                
%                 Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_1m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(2,:);
%                 Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_3m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(4,:);
%                 Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_5m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(6,:);
%                 Temperature.tanks(tnk_plot).Reduction(hh).Temp_1m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(2,:);
%                 Temperature.tanks(tnk_plot).Reduction(hh).Temp_3m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(4,:);
%                 Temperature.tanks(tnk_plot).Reduction(hh).Temp_5m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(6,:);
%                 
            end
        end
        
    end
end

if Type == 2
    for gg = 1:length(Sn2) % Number of variable scenarios (e.g. inflow)
        for hh = 1:Sn+1 % Number of reductions
            for ii = 1:N % All the acceptable simulations
                
                % Load simulated
                formatSpec = '%d\\QE2_%d_%d_%d';
                str = sprintf(formatSpec,ScenarioNum,gg,hh,ii);
                load(['\\lancs\luna\FST\LEC\Users\exleyg\Parameter_Distributions\QE2_Scenario_',str],'MyLake_results');
                
                for tnk_plot = 1:Tk % which tanks to plot
                    
                    % Chlorophyll
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla = MyLake_results.basin1.concentrations.tanks(tnk_plot).Chl; % All Chlorophyll (by species)
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla_Sum = sum(cat(3,Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{:}),3); % Total chlorophyll
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla_Top_Mean = nanmean(Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla_Sum(1:6,:),1); % Surface to 5m mean total chlorophyll
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Diatoms = Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{1} + Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{2};
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Green = Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{3} + Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{4};
                    Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Cyano = Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{5} + Chlorophyll.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Chla{6};
                    
                    % Temperature
                    Temperature.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Temp = MyLake_results.basin1.tanks(tnk_plot).T; % All temperatures
                    Temperature.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Temp_Top_Mean = nanmean(Temperature.tanks(tnk_plot).Scenario(gg).Reduction(hh).RunNum(ii).Temp(1:6,:),1); % Surface to 5m mean temperature
                    
                    % Mixed depth
                    Mixed_Depth.tanks(tnk_plot).Scenario(gg).Sim(hh).Mxd(ii,:) = MyLake_results.basin1.tanks(tnk_plot).MixStat(12,:);
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(tnk_plot).Scenario(gg).Sim(hh).Mxd(ii,:);
                    Mxd_lyr(isnan(Mxd_lyr)) = 17;
                    Mxd_lyr = 17-Mxd_lyr;
                    Mixed_Depth.tanks(tnk_plot).Scenario(gg).Sim(hh).metric(ii,:) = trapz(Mxd_lyr(ii,:));
                    
                    %             Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_1m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(2,:);
                    %             Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_3m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(4,:);
                    %             Chlorophyll.tanks(tnk_plot).Reduction(hh).Chl_5m(ii,:) = Chlorophyll.tanks(tnk_plot).Reduction(hh).RunNum(ii).Chla_Sum(6,:);
                    %             Temperature.tanks(tnk_plot).Reduction(hh).Temp_1m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(2,:);
                    %             Temperature.tanks(tnk_plot).Reduction(hh).Temp_3m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(4,:);
                    %             Temperature.tanks(tnk_plot).Reduction(hh).Temp_5m(ii,:) = Temperature.tanks(tnk_plot).Reduction(hh).RunNum(ii).Temp(6,:);
                    %
                end
            end
            
        end
        
    end
end

end
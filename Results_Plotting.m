function [PLP_rec, Mxd_tot, P1_rec] = Results_Plotting(Chlorophyll, Temperature, Mixed_Depth, m_start, m_stop, N, ScenarioNum, Sn, Sn2, Tk)
%Results_Plotting Plot output from MyLake runs
%   Plotting routine for MyLake runs

filepath = '\\lancs\\luna\\FST\\LEC\\Users\\exleyg\\Parameter_Distributions'; % Set the file path to save figures
GLF = ones(N,1); % Needs weightings but set to one here

%% Read in LoA
N_Par = 9; % Number of parameters

% Load limits of acceptability
load 'IO\'Chl_Obs_LoA_mod
Chl_Obs_LoA = Chl_Obs_LoA_mod;
load 'IO\'Temp_Obs_LoA
load 'IO\'Mxd_Obs_LoA
SimDate = datenum(m_start):datenum(m_stop);

% Time series for observed chlorophyll saved into a cell array
for kk = 2:2:(N_Par*2)
    TimeS_Chl = (Chl_Obs_LoA(~isnan(Chl_Obs_LoA(:,kk))));
    if kk <=2
        TS_Chl = {TimeS_Chl};
    else
        TS_Chl = [TS_Chl TimeS_Chl];
    end
end

% Time series for observed temperature saved into a cell array
for kk = 2:2:(N_Par*2)
    TimeS_Temp = (Temp_Obs_LoA(~isnan(Temp_Obs_LoA(:,kk))));
    if kk <=2
        TS_Temp = {TimeS_Temp};
    else
        TS_Temp = [TS_Temp TimeS_Temp];
    end
end

% Time series for estimated mixed depth saved into a cell array
TimeS_Mxd = (Mxd_Obs_LoA(~isnan(Mxd_Obs_LoA(:,2))));
        TS_Mxd = {TimeS_Mxd};


%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 28);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:AB366";

% Specify column names and types
opts.VariableNames = ["Date", "CHLOROPHYLL_ONSITE1", "CHLOROPHYLL_ONSITE3", "CHLOROPHYLL_ONSITE5", "CHLOROPHYLL_ONSITE7", "CHLOROPHYLL_ONSITE9", "CHLOROPHYLL_ONSITE11", "CHLOROPHYLL_ONSITE13", "CHLOROPHYLL_ONSITE15", "CHLOROPHYLL_ONSITE17", "DISSOLVEDOXYGEN_ONSITE1", "DISSOLVEDOXYGEN_ONSITE3", "DISSOLVEDOXYGEN_ONSITE5", "DISSOLVEDOXYGEN_ONSITE7", "DISSOLVEDOXYGEN_ONSITE9", "DISSOLVEDOXYGEN_ONSITE11", "DISSOLVEDOXYGEN_ONSITE13", "DISSOLVEDOXYGEN_ONSITE15", "DISSOLVEDOXYGEN_ONSITE17", "TEMP_ONSITE1", "TEMP_ONSITE3", "TEMP_ONSITE5", "TEMP_ONSITE7", "TEMP_ONSITE9", "TEMP_ONSITE11", "TEMP_ONSITE13", "TEMP_ONSITE15", "TEMP_ONSITE17"];
opts.SelectedVariableNames = ["Date", "CHLOROPHYLL_ONSITE1", "CHLOROPHYLL_ONSITE3", "CHLOROPHYLL_ONSITE5", "CHLOROPHYLL_ONSITE7", "CHLOROPHYLL_ONSITE9", "CHLOROPHYLL_ONSITE11", "CHLOROPHYLL_ONSITE13", "CHLOROPHYLL_ONSITE15", "CHLOROPHYLL_ONSITE17", "DISSOLVEDOXYGEN_ONSITE1", "DISSOLVEDOXYGEN_ONSITE3", "DISSOLVEDOXYGEN_ONSITE5", "DISSOLVEDOXYGEN_ONSITE7", "DISSOLVEDOXYGEN_ONSITE9", "DISSOLVEDOXYGEN_ONSITE11", "DISSOLVEDOXYGEN_ONSITE13", "DISSOLVEDOXYGEN_ONSITE15", "DISSOLVEDOXYGEN_ONSITE17", "TEMP_ONSITE1", "TEMP_ONSITE3", "TEMP_ONSITE5", "TEMP_ONSITE7", "TEMP_ONSITE9", "TEMP_ONSITE11", "TEMP_ONSITE13", "TEMP_ONSITE15", "TEMP_ONSITE17"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");

% Import the data
ObsQEIIDaily = readtable("C:\Users\exleyg\OneDrive - Lancaster University\MyLake Project\MyLake_MultiTank_V1_29_06_2020 - QEII_2\MyLake_MultiTank_V1_TwoTanks\ObsQEIIDaily.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

% Table to matrix
Chl_Obs = table2array(ObsQEIIDaily(:,2:10));
Temp_Obs = table2array(ObsQEIIDaily(:,20:28));
load 'IO\'Mxd_Obs;
Mxd_Obs = Mxd_Obs;

% %% Means
% 
% % Plot mean top chlorophyll
% for ii = 1:Sn+1
% Chl_Top_Mean_Mean(ii,:) = nanmean([Chlorophyll.tanks(1).Reduction(1).RunNum(ii).Chla_Top_Mean],3);
% end
% 
% figure(1)
% plot(SimDate, Chl_Top_Mean_Mean);
% datetick('x','mmm');
% xlabel('Date');
% ylabel('Mean Surface Chlorophyll-a (µg L^{-1})');
% lgd = legend('0','10','20','30','40','50','60','70','80','90','100');
% lgd;
% title(lgd,'FPV Coverage (%)');
% 
% % Plot mean top temperature
% for ii = 1:Sn
% Temp_Top_Mean_Mean(ii,:) = nanmean([Temperature_1.tanks(1).Reduction(1).RunNum(ii).Temp_Top_Mean],3);  
% end
% 
% figure(2)
% plot(SimDate, Temp_Top_Mean_Mean);
% datetick('x','mmm')
% xlabel('Date')
% ylabel('Water Temperature (^{o}C)')
% lgd = legend('0','10','20','30','40','50','60','70','80','90','100')
% title(lgd,'FPV Coverage (%)')

%% Box plot of all scenarios
% Chlorophyll
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        ppp = Sn2 == 1;
        for ll = [2,4,6] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    % dat(jj,:,:) = Chlorophyll.tanks(1).Reduction(jj).Chl_5m(:,:);
                    % Chlorophyll.tanks(1).Reduction(1).RunNum(1).Chla_Sum
                    if length(Sn2)== 1
                    Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    dat2(ii,:) = Chlorophyll.tanks(mm).Reduction(1).RunNum(ii).Chla_Sum(2,:);    
                    else
                    Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    dat2(ii,:) = Chlorophyll.tanks(mm).Scenario(ppp).Reduction(1).RunNum(ii).Chla_Sum(2,:);
                    end
                end
                dat(jj,:,:) = Chl_lyr;
                TS = find(dat2(2,:) == max(dat2(2,:)));
            end
            figure()
            for kk = 1:Sn+1
                P1 = DeterminePercents(GLF,[0.000001,0.025,0.5,0.975,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
                PLP = mean(P1(TS-5:TS+5,:)); % Get mean of percentiles in the 11 day range
                if length(Sn2) == 1
            PLP_rec.tanks(mm).Chl(ll-1).PLP(kk,1:5) = PLP;
            P1_rec.tanks(mm).Chl(ll-1).P1(kk).P1 = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(kk,1:5) = PLP;
            end
                yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
                lentp = 1;
                %          Median ,pecentile 1,percentile 3,max,min
                SPIRboxplot(PLP(3),PLP(2),PLP(4),PLP(5),PLP(1),lentp,yy,0,'+',1,kk*10);
                % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
                hold on
                % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
            end
           
            if length(Sn2) == 1
                 % Percentage change
            PLP_rec.tanks(mm).Chl(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Chl(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Chl(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Chl(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Chl(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Chl(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Depth: %d m, Tank %d';
                A1 = ll-1;
                A2 = mm;
                str = sprintf(formatSpec,A1,A2);
                title({'Chlorophyll-a distribution on day of maximum Chlorophyll-a', str})
                            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
            else
                % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Chl(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
                A1 = ll-1;
                A2 = mm;
                A3 = Sn2(pp);
                str = sprintf(formatSpec,A1,A2,A3);
                title({'Chlorophyll-a distribution on day of maximum Chlorophyll-a', str})
                xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
            end

            xlabel('Floating solar array coverage (%)')
            ylabel('Chlorophyll-a distribution (µg L^{-1})')
            %ylim([0 35])
            set(gca,'FontSize',12)
            % Save
            if length(Sn2) == 1
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Chl_Box_%dm_Tank%d';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                str = sprintf(formatSpec,B1,B2,B3);
                saveas(gcf,[filepath,str]);
            else
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Chl_Box_%dm_Tank%d_InfwF%g';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                B4 = pp;
                str = sprintf(formatSpec,B1,B2,B3,B4);
                saveas(gcf,[filepath,str]);
            end
        end
    end
end

%% Maximum Temperature
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        ppp = Sn2 == 1;
        for ll = [2,4,6] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                    Temp_lyr(ii,:) = Temperature.tanks(mm).Reduction(jj).RunNum(ii).Temp(ll,:);
                    dat2(ii,:) = Temperature.tanks(mm).Reduction(1).RunNum(ii).Temp(2,:);    
                    else
                    Temp_lyr(ii,:) = Temperature.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Temp(ll,:);
                    dat2(ii,:) = Temperature.tanks(mm).Scenario(ppp).Reduction(1).RunNum(ii).Temp(2,:);
                    end
                end
                dat(jj,:,:) = Temp_lyr;
                TS = find(dat2(2,:) == max(dat2(2,:)));
            end
            figure()
            for kk = 1:Sn+1
                P1 = DeterminePercents(GLF,[0.000001,0.025,0.5,0.975,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
                PLP = mean(P1(TS-5:TS+5,:)); % Get mean of percentiles in the 11 day range
                            if length(Sn2) == 1
            PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(kk,1:5) = PLP;
            P1_rec.tanks(mm).WaterTemp(ll-1).P1(kk).P1 = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(kk,1:5) = PLP;
            end
                yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
                lentp = 1;
                %          Median ,pecentile 1,percentile 3,max,min
                SPIRboxplot(PLP(3),PLP(2),PLP(4),PLP(5),PLP(1),lentp,yy,0,'+',1,kk*10);
                % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
                hold on
                % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
            end
            
            if length(Sn2) == 1
                 % Percentage change
            PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(:,2))./PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(:,4))./PLP_rec.tanks(mm).WaterTemp_Max(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Depth: %d m, Tank %d';
                A1 = ll-1;
                A2 = mm;
                str = sprintf(formatSpec,A1,A2);
                title({'Temperature distribution on day of maximum temperature', str})
                            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
            else
                % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Max(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
                A1 = ll-1;
                A2 = mm;
                A3 = Sn2(pp);
                str = sprintf(formatSpec,A1,A2,A3);
                title({'Temperature distribution on day of maximum temperature', str})
                      xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
            end
            xlabel('Floating solar array coverage (%)')
            ylabel('Water Temperature (^{o}C)')
            %ylim([18 24])
            set(gca,'FontSize',12)
            % Save
            if length(Sn2) == 1
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Box_%dm_Tank%d';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                str = sprintf(formatSpec,B1,B2,B3);
                saveas(gcf,[filepath,str]);
            else
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Box_%dm_Tank%d_InfwF%g';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                B4 = pp;
                str = sprintf(formatSpec,B1,B2,B3,B4);
                saveas(gcf,[filepath,str]);
            end
        end
    end
end

%% Minimum water Temperature
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        ppp = Sn2 == 1;
        for ll = [2,4,6] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                    Temp_lyr(ii,:) = Temperature.tanks(mm).Reduction(jj).RunNum(ii).Temp(ll,:);
                    dat2(ii,:) = Temperature.tanks(mm).Reduction(1).RunNum(ii).Temp(2,:);    
                    else
                    Temp_lyr(ii,:) = Temperature.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Temp(ll,:);
                    dat2(ii,:) = Temperature.tanks(mm).Scenario(ppp).Reduction(1).RunNum(ii).Temp(2,:);
                    end
                end
                dat(jj,:,:) = Temp_lyr;
                TS = find(dat2(2,:) == min(dat2(2,:)));
            end
            figure()
            for kk = 1:Sn+1
                P1 = DeterminePercents(GLF,[0.000001,0.025,0.5,0.975,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
                PLP = mean(P1(TS-5:TS+5,:)); % Get mean of percentiles in the 11 day range
                if length(Sn2) == 1
            PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(kk,1:5) = PLP;
            else
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(kk,1:5) = PLP;
                end
                yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
                lentp = 1;
                %          Median ,pecentile 1,percentile 3,max,min
                SPIRboxplot(PLP(3),PLP(2),PLP(4),PLP(5),PLP(1),lentp,yy,0,'+',1,kk*10);
                % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
                hold on
                % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
            end
            
            if length(Sn2) == 1
                % Percentage change
            PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(:,2))./PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(:,4))./PLP_rec.tanks(mm).WaterTemp_Min(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Depth: %d m, Tank %d';
                A1 = ll-1;
                A2 = mm;
                str = sprintf(formatSpec,A1,A2);
                title({'Temperature distribution on day of minimum temperature', str})
                            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
            else
            % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).WaterTemp_Min(ll-1).PLP(1:end-1,4);
               
                formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
                A1 = ll-1;
                A2 = mm;
                A3 = Sn2(pp);
                str = sprintf(formatSpec,A1,A2,A3);
                title({'Temperature distribution on day of minimum temperature', str})
                      xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
            end
            xlabel('Floating solar array coverage (%)')
            ylabel('Water Temperature (^{o}C)')
            %ylim([18 24])
            set(gca,'FontSize',12)
            % Save
            if length(Sn2) == 1
                formatSpec = '\\QE2_Scenario_%d\\Figures\\TempMin_Box_%dm_Tank%d';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                str = sprintf(formatSpec,B1,B2,B3);
                saveas(gcf,[filepath,str]);
            else
                formatSpec = '\\QE2_Scenario_%d\\Figures\\TempMin_Box_%dm_Tank%d_InfwF%g';
                B1 = ScenarioNum;
                B2 = ll-1;
                B3 = mm;
                B4 = pp;
                str = sprintf(formatSpec,B1,B2,B3,B4);
                saveas(gcf,[filepath,str]);
            end
        end
    end
end

%% Mixed depth
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for jj = 1:Sn+1
            for ii = 1:N
                if length(Sn2)== 1
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(mm).Sim(jj).metric(ii,:);
                    Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Sim(jj).Mxd),2); % Number of stratified days
                else
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).metric(ii,:);
                    Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd),2); % Number of stratified days
                end
            end
            Mxd_lyr_norm = Mxd_lyr/trapz(17*ones(365,1)); % Normalise
            dat(jj,:,:) = Mxd_lyr_norm;
            
            end
            figure()
            for kk = 1:Sn+1
                P1 = DeterminePercents(GLF,[0.000001,0.025,0.5,0.975,0.9999],reshape(dat(kk,:,:),N,1)); % calculate percentiles
                % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
                            if length(Sn2) == 1
            PLP_rec.tanks(mm).Mxd(ll-1).PLP(kk,1:5) = P1;
            P1_rec.tanks(mm).Mxd(ll-1).P1(kk).P1 = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(kk,1:5) = P1;
            end
                yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
                lentp = 1;
                %          Median ,pecentile 1,percentile 3,max,min
                SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
                % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
                hold on
                % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
            end
            
            if length(Sn2) == 1
                % Percentage change
            PLP_rec.tanks(mm).Mxd(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Mxd(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Mxd(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Mxd(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Mxd(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Mxd(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Tank %d';
                A2 = mm;
                str = sprintf(formatSpec,A2);
                title({'Mixed depth metric', str})
                            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
            else
                % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Mxd(ll-1).PLP(1:end-1,4);
           
                formatSpec = 'Tank: %d, Inflow Factor: %g';
                A2 = mm;
                A3 = Sn2(pp);
                str = sprintf(formatSpec,A2,A3);
                title({'Mixed depth metric', str})
                      xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
            end
            xlabel('Floating solar array coverage (%)')
            ylabel('Mixed depth metric (normalised)')
            %ylim([18 24])
            
            % Save
            if length(Sn2) == 1
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Mxd_Box_Tank%d';
                B1 = ScenarioNum;
                B3 = mm;
                str = sprintf(formatSpec,B1,B3);
                saveas(gcf,[filepath,str]);
            else
                formatSpec = '\\QE2_Scenario_%d\\Figures\\Mxd_Box_Tank%d_InfwF%g';
                B1 = ScenarioNum;
                B3 = mm;
                B4 = pp;
                str = sprintf(formatSpec,B1,B3,B4);
                saveas(gcf,[filepath,str]);
            end
        end
end

%% Stratification duration
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for jj = 1:Sn+1
            for ii = 1:N
                if length(Sn2)== 1
                    Mixed_Depth.tanks(mm).Sim(jj).Strat.Mxd(ii,:) = Mixed_Depth.tanks(mm).Sim(jj).Mxd(ii,:);
                    Mixed_Depth.tanks(mm).Sim(jj).Strat.Mxd(isnan(Mixed_Depth.tanks(mm).Sim(jj).Strat.Mxd)) = 17;
                    RLE = rle(Mixed_Depth.tanks(mm).Sim(jj).Strat.Mxd(ii,:)<17);
                    if ~isempty(max(RLE(RLE(:,1)==1,2)))
                        [val,~] = max(RLE(RLE(:,1)==1,2));
                        idx = find(RLE(:,2)== val, 1 );
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.Max_Dur(ii,:) = val;
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.Start(ii,:) = RLE(idx,3);
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.End(ii,:) = RLE(idx,4);
                    else
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.Max_Dur(ii,:) = 0;
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.Start(ii,:) = NaN;
                        Mixed_Depth.tanks(mm).Sim(jj).Strat.End(ii,:) = NaN;
                    end
                    %Strat_dur(ii,:) = Mixed_Depth.tanks(mm).Sim(jj).Strat(ii).Max_Dur;
                else
                    Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Mxd(ii,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Mxd(ii,:);
                    Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Mxd(isnan(Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Mxd)) = 17;
                    RLE = rle(Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Mxd(ii,:)<17);
                    if ~isempty(max(RLE(RLE(:,1)==1,2)))
                        [val,~] = max(RLE(RLE(:,1)==1,2));
                        idx = find(RLE(:,2)== val, 1 );
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Max_Dur(ii,:) = val;
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Start(ii,:) = RLE(idx,3);
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.End(ii,:) = RLE(idx,4);
                    else
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Max_Dur(ii,:) = 0;
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Start(ii,:) = NaN;
                        Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.End(ii,:) = NaN;
                    end
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).metric(ii,:);
                end
            end
            if length(Sn2)== 1
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Sim(jj).Strat.Max_Dur;
            else
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Max_Dur;
            end
        end
        figure()
        for kk = 1:Sn+1
            P1 = DeterminePercents(GLF,[0.000001,0.025,0.5,0.975,0.9999],reshape(dat(kk,:,:),N,1)); % calculate percentiles
            % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
            if length(Sn2) == 1
            PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(kk,1:5) = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(kk,1:5) = P1;
            end
            yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
            lentp = 1;
            %          Median ,pecentile 1,percentile 3,max,min
            SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
            % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
            hold on
            % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
        end
        
        if length(Sn2) == 1
            % Percentage change
            PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Strat_Max_Dur(ll-1).PLP(1:end-1,4);
            
            formatSpec = 'Tank %d';
            A2 = mm;
            str = sprintf(formatSpec,A2);
            title({'Stratification duration', str})
            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
        else
            % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Strat_Max_Dur(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank: %d, Inflow Factor: %g';
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A2,A3);
            title({'Stratification duration', str})
            xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
        end
        xlabel('Floating solar array coverage (%)')
        ylabel('Duration (days)')
        %ylim([18 24])
        set(gca,'FontSize',16)
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratDur_Box_Tank%d';
            B1 = ScenarioNum;
            B3 = mm;
            str = sprintf(formatSpec,B1,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratDur_Box_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Polymixis
% Longest stratified period divide total number of stratified days
clear dat Mxd_tot
for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor 3
    for jj = 1:Sn+1
        if length(Sn2)== 1
            dat(jj,:,:) = Mixed_Depth.tanks(1).Sim(jj).Strat.Max_Dur;
            Mxd_temp = Mixed_Depth.tanks(1).Sim(jj).Mxd;
            Mxd_temp(Mxd_temp==17)= NaN;
            Mxd_tot(:,jj) = sum(~isnan(Mxd_temp),2); % Number of stratified days
        else
            dat(jj,:,:) = Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Strat.Max_Dur;
            Mxd_temp = Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd;
            Mxd_temp(Mxd_temp==17)= NaN;
            Mxd_tot(:,jj) = sum(~isnan(Mxd_temp),2); % Number of stratified days
        end
    end
    Polymixis = dat./Mxd_tot';
end


clear dat
for mm = 1
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        figure()
        for kk = 1:Sn+1
            daystrattemp = Polymixis(kk,(~isnan(Polymixis(kk,:,:))));
            if ~isempty(daystrattemp)
            P1 = DeterminePercents(GLF(~isnan(Polymixis(kk,:))),[0.000001,0.025,0.5,0.975,0.9999],daystrattemp'); % calculate percentiles
            % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
            PLP_rec.tanks(mm).Polymixis(ll-1).PLP(kk,1:5) = P1;
            yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
            lentp = 1;
            %          Median ,pecentile 1,percentile 3,max,min
            SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
            text(kk*10,P1(5)+0.03,num2str(length(daystrattemp)));
            else
                plot(kk*10,1,'r*','MarkerSize',20)
            end
            % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
            hold on
            % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
        end
                % Percentage change
%             PLP_rec.tanks(mm).Polymixis(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Polymixis(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Polymixis(ll-1).PLP(1:end-1,2);
%             PLP_rec.tanks(mm).Polymixis(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Polymixis(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Polymixis(ll-1).PLP(1:end-1,4);
%            
        if length(Sn2) == 1
            formatSpec = 'Tank %d';
            A2 = mm;
            str = sprintf(formatSpec,A2);
            %title({'Polymixis', str})
            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
        else
            formatSpec = 'Tank: %d, Inflow Factor: %g';
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A2,A3);
            title({'Polymixis', str})
            xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
        end
        xlabel('Floating solar array coverage (%)')
        ylabel('Polymixis metric')
        %ylim([18 24])
        set(gca,'FontSize',12)
                yticks([0:0.1:1])
        yticklabels({'0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'})
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Polymixis_Box_Tank%d';
            B1 = ScenarioNum;
            B3 = mm;
            str = sprintf(formatSpec,B1,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Polymixis_Box_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Stratification start day
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for jj = 1:Sn+1
            if length(Sn2)== 1
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Sim(jj).Strat.Start;
            else
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.Start;    
            end
        end
        figure()
        for kk = 1:Sn+1
            daystrattemp = dat(kk,(~isnan(dat(kk,:,:))),:);
            if ~isempty(daystrattemp)
            P1 = DeterminePercents(GLF(~isnan(dat(kk,:,:))),[0.000001,0.025,0.5,0.975,0.9999],daystrattemp'); % calculate percentiles
            % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
            if length(Sn2) == 1
            PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(kk,1:5) = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(kk,1:5) = P1;
            end
            yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
            lentp = 1;
            %          Median ,pecentile 1,percentile 3,max,min
            SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
            text(kk*10,P1(5)+10,num2str(length(daystrattemp)),'FontSize',14);
            else
                plot(kk*10,1,'r*','MarkerSize',20)
            end
            % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
            hold on
            % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
        end
                
        if length(Sn2) == 1
            % Percentage change
            PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Strat_Start(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank %d';
            A2 = mm;
            str = sprintf(formatSpec,A2);
            title({'Stratification onset', str})
            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
        else
            % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Strat_Start(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank: %d, Inflow Factor: %g';
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A2,A3);
            title({'Stratification onset', str})
            xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
        end
        xlabel('Floating solar array coverage (%)')
        ylabel('Stratification Onset (day of year)')
        %ylim([18 24])
        set(gca,'FontSize',16)
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratOnset_Box_Tank%d';
            B1 = ScenarioNum;
            B3 = mm;
            str = sprintf(formatSpec,B1,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratOnset_Box_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Stratification end day
clear dat daystrattemp
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
           for jj = 1:Sn+1
            if length(Sn2)== 1
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Sim(jj).Strat.End;
            else
            Mxd_tot(:,jj) = sum(~isnan(Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd),2); % Number of stratified days
            dat(jj,:,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).Strat.End;    
            end
        end
        figure()
        for kk = 1:Sn+1
            daystrattemp = dat(kk,(~isnan(dat(kk,:,:))),:);
            if ~isempty(daystrattemp)
            P1 = DeterminePercents(GLF(~isnan(dat(kk,:,:))),[0.000001,0.025,0.5,0.975,0.9999],daystrattemp'); % calculate percentiles
            % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
                        if length(Sn2) == 1
            PLP_rec.tanks(mm).Strat_End(ll-1).PLP(kk,1:5) = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(kk,1:5) = P1;
                        end
            yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
            lentp = 1;
            %          Median ,pecentile 1,percentile 3,max,min
            SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
            text(kk*10,P1(5)+10,num2str(length(daystrattemp)),'FontSize',14);
            else
                plot(kk*10,1,'r*','MarkerSize',20)
            end
            % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
            hold on
            % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
        end
               
        if length(Sn2) == 1
             % Percentage change
            PLP_rec.tanks(mm).Strat_End(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Strat_End(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Strat_End(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Strat_End(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Strat_End(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Strat_End(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank %d';
            A2 = mm;
            str = sprintf(formatSpec,A2);
            title({'Stratification overturn', str})
            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
        else
            % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank: %d, Inflow Factor: %g';
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A2,A3);
            title({'Stratification overturn', str})
            xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
        end
        xlabel('Floating solar array coverage (%)')
        ylabel('Stratification Overturn (day of year)')
        %ylim([18 24])
        set(gca,'FontSize',16)
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratOver_Box_Tank%d';
            B1 = ScenarioNum;
            B3 = mm;
            str = sprintf(formatSpec,B1,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratOver_Box_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Total number of stratified days (cumulative stratification)
clear dat daystrattemp dat_temp
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for jj = 1:Sn+1
            if length(Sn2)== 1
                dat_temp = Mixed_Depth.tanks(1).Sim(jj).Mxd;
                dat_temp(dat_temp==17)= NaN;
                dat(jj,:) = sum(~isnan(dat_temp),2)'; % Number of stratified days
                %dat(jj,:) = sum(~isnan(Mixed_Depth.tanks(mm).Sim(jj).Mxd),2)'; % Number of stratified days
                %dat(jj,:,:) = Mixed_Depth.tanks(mm).Sim(jj).Strat.End;
            else
                dat_temp = Mixed_Depth.tanks(1).Scenario(pp).Sim(jj).Mxd;
                dat_temp(dat_temp==17)= NaN;
                dat(jj,:) = sum(~isnan(dat_temp),2)'; % Number of stratified days
            end
        end
        figure()
        for kk = 1:Sn+1
            daystrattemp = dat(kk,(~isnan(dat(kk,:,:))),:);
            if ~isempty(daystrattemp)
            P1 = DeterminePercents(GLF(~isnan(dat(kk,:,:))),[0.000001,0.025,0.5,0.975,0.9999],daystrattemp'); % calculate percentiles
            % P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
            if length(Sn2) == 1
            PLP_rec.tanks(mm).Strat_Cumul(ll-1).PLP(kk,1:5) = P1;
            else
            PLP_rec.tanks(mm).Scenario(pp).Strat_Cumul(ll-1).PLP(kk,1:5) = P1;
            end
            yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
            lentp = 1;
            %          Median ,pecentile 1,percentile 3,max,min
            SPIRboxplot(P1(3),P1(2),P1(4),P1(5),P1(1),lentp,yy,0,'+',1,kk*10);
            % text(kk*10,P1(5)+10,num2str(length(daystrattemp)));
            else
                plot(kk*10,1,'r*','MarkerSize',20)
            end
            % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
            hold on
            % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
        end
               
        if length(Sn2) == 1
             % Percentage change
            PLP_rec.tanks(mm).Strat_End(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Strat_End(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Strat_End(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Strat_End(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Strat_End(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Strat_End(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank %d';
            A2 = mm;
            str = sprintf(formatSpec,A2);
            title({'Cumulative stratification', str})
            xticks([10:10:110])
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
        else
            % Percentage change
            PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(2:Sn+1,6) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(:,2))./PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(1:end-1,2);
            PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(2:Sn+1,7) = 100*diff(PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(:,4))./PLP_rec.tanks(mm).Scenario(pp).Strat_End(ll-1).PLP(1:end-1,4);
           
            formatSpec = 'Tank: %d, Inflow Factor: %g';
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A2,A3);
            title({'Cumulative stratification', str})
            xticks([10:10:50])
            xticklabels({'0','25','50','75','100'})
        end
        xlabel('Floating solar array coverage (%)')
        ylabel('Cumulative stratification (days)')
        %ylim([18 24])
        set(gca,'FontSize',16)
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratTotal_Box_Tank%d';
            B1 = ScenarioNum;
            B3 = mm;
            str = sprintf(formatSpec,B1,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\StratTotal_Box_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end
%%
% for mm = 1:Tk
%     for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
%         ppp = Sn2 == 1;
%             for jj = 1:Sn+1
%                 for ii = 1:N
%                     if length(Sn2)== 1
%                     Temp_lyr(ii,:) = Temperature.tanks(mm).Reduction(jj).RunNum(ii).Mxd;
%                     dat2(ii,:) = Temperature.tanks(mm).Reduction(1).RunNum(ii).Temp(2,:);    
%                     else
%                     Temp_lyr(ii,:) = Temperature.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Temp(ll,:);
%                     dat2(ii,:) = Temperature.tanks(mm).Scenario(ppp).Reduction(1).RunNum(ii).Temp(2,:);
%                     end
%                 end
%                 dat(jj,:,:) = Temp_lyr;
%                 TS = find(dat2(2,:) == max(dat2(2,:)));
%             end
%             figure()
%             for kk = 1:Sn+1
%                 P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
%                 PLP = mean(P1(TS-5:TS+5,:)); % Get mean of percentiles in the 11 day range
%                 yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
%                 lentp = 1;
%                 %          Median ,pecentile 1,percentile 3,max,min
%                 SPIRboxplot(PLP(3),PLP(2),PLP(4),PLP(5),PLP(1),lentp,yy,0,'+',1,kk*10);
%                 % text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
%                 hold on
%                 % plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
%             end
%             if length(Sn2) == 1
%                 formatSpec = 'Depth: %d m, Tank %d';
%                 A1 = ll-1;
%                 A2 = mm;
%                 str = sprintf(formatSpec,A1,A2);
%                 title({'Temperature distribution on day of maximum temperature', str})
%                             xticks([10:10:110])
%             xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%             else
%                 formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
%                 A1 = ll-1;
%                 A2 = mm;
%                 A3 = Sn2(pp);
%                 str = sprintf(formatSpec,A1,A2,A3);
%                 title({'Temperature distribution on day of maximum temperature', str})
%                       xticks([10:10:50])
%             xticklabels({'0','25','50','75','100'})
%             end
%             xlabel('Floating solar array coverage (%)')
%             ylabel('Water Temperature (^{o}C)')
%             %ylim([18 24])
%             
%             % Save
%             if length(Sn2) == 1
%                 formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Box_%dm_Tank%d';
%                 B1 = ScenarioNum;
%                 B2 = ll-1;
%                 B3 = mm;
%                 str = sprintf(formatSpec,B1,B2,B3);
%                 saveas(gcf,[filepath,str]);
%             else
%                 formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Box_%dm_Tank%d_InfwF%g';
%                 B1 = ScenarioNum;
%                 B2 = ll-1;
%                 B3 = mm;
%                 B4 = pp;
%                 str = sprintf(formatSpec,B1,B2,B3,B4);
%                 saveas(gcf,[filepath,str]);
%             end
%         end
%     end
% end

%% Envelope plots
% Chlorophyll
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for ll = [2] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                        Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    else
                        Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    end
                end
                dat(jj,:,:) = Chl_lyr;
            end
            x = 1:365; % set x axis
            clr = distinguishable_colors(11); % get some useful colours
            % plot uncertainty envelopes
            figure()
            for kk = 1:2:Sn+1 %:2:11 % Define which layers to plot here
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),clr(kk,:)); % calculate percentiles)
                hold on
            end
        end
        datetick('x','mmm')
        xlabel('Date')
        if length(Sn2) == 1
            formatSpec = 'Depth: %d m, Tank: %d';
            A1 = ll-1;
            A2 = mm;
            str = sprintf(formatSpec,A1,A2);
            title({'Total Chlorophyll-a', str})
        else
            formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
            A1 = ll-1;
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A1,A2,A3);
            title({'Total Chlorophyll-a', str})
        end
        ylabel('Chlorophyll-a (µg L^{-1})')
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Chl_Envel_%dm_Tank%d';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            str = sprintf(formatSpec,B1,B2,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Chl_Envel_%dm_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B2,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

% Temperature
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for ll = [2] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                        Temp_lyr(ii,:) = Temperature.tanks(mm).Reduction(jj).RunNum(ii).Temp(ll,:);
                    else
                        Temp_lyr(ii,:) = Temperature.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Temp(ll,:);
                    end
                end
                dat(jj,:,:) = Temp_lyr;
            end
            x = 1:365; % set x axis
            clr = distinguishable_colors(11); % get some useful colours
            % plot uncertainty envelopes
            figure()
            for kk = 1:2:Sn+1 %:2:11 % Define which layers to plot here
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),clr(kk,:)); % calculate percentiles)
                hold on
            end
        end
        datetick('x','mmm')
        xlabel('Date')
        if length(Sn2) == 1
            formatSpec = 'Depth: %d m, Tank: %d';
            A1 = ll-1;
            A2 = mm;
            str = sprintf(formatSpec,A1,A2);
            title({'Water temperature', str})
        else
            formatSpec = 'Depth: %d m, Tank: %d, Inflow Factor: %g';
            A1 = ll-1;
            A2 = mm;
            A3 = Sn2(pp);
            str = sprintf(formatSpec,A1,A2,A3);
            title({'Water temperature', str})
        end
        ylabel('Water Temperature (^{o}C)')
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Envel_%dm_Tank%d';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            str = sprintf(formatSpec,B1,B2,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\Temp_Envel_%dm_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B2,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Algal species 
% Subplots
clear dat dat2 dat3
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for ll = [2] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                        Dia_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Diatoms(ll,:);
                        Gre_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Green(ll,:);
                        Cyn_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Cyano(ll,:);
                    else
                        Dia_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Diatoms(ll,:);
                        Gre_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Green(ll,:);
                        Cyn_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Cyano(ll,:);
                    end
                end
                dat(jj,:,:) = Dia_lyr;
                dat2(jj,:,:) = Gre_lyr;
                dat3(jj,:,:) = Cyn_lyr;
            end
            x = 1:365; % set x axis
            clr = distinguishable_colors(3); % get some useful colours
            % plot uncertainty envelopes
            figure()
            for kk = 1:1:Sn+1 %:2:11 % Define which layers to plot here
                subplot(6,2,kk)
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),clr(1,:)); % calculate percentiles)
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat2(kk,:,:),N,365)),clr(2,:));
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat3(kk,:,:),N,365)),clr(3,:));
                hold on
                datetick('x','mmm')
                %xlabel('Date')
                if length(Sn2) == 1
                    formatSpec = 'Floating Solar Coverage: %d %%'; %, Depth: %d m, Tank: %d';
                    A1 = ll-1;
                    A2 = mm;
                    A3_1 = [0:10:100];
                    A3 = A3_1(kk);
                    str = sprintf(formatSpec,A3); %,A1,A2);
                    title(str)
                    ylim([0 ceil(max(dat,[],'all'))]);
                else
                    formatSpec = 'Floating Solar Coverage: %d %%'; %'Depth: %d m, Tank: %d, Inflow Factor: %g';
                    A1 = ll-1;
                    A2 = mm;
                    A3 = Sn2(pp);
                    A4_1 = [0:25:100];
                    A4 = A4_1(kk);
                    str = sprintf(formatSpec,A4); %,A1,A2);
                    title(str)
                    ylim([0 ceil(max(dat,[],'all'))]);
                end
                % ylabel('Chlorophyll-a (µg L^{-1})')
            end
        end
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpecies_Envel_%dm_Tank%d';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            str = sprintf(formatSpec,B1,B2,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpecies_Envel_%dm_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B2,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Algal species as proportions
% Subplots
clear dat dat2 dat3
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for ll = [2] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                        Cyn_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Cyano(ll,:)./Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                        Dia_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Diatoms(ll,:)./Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                        Gre_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Green(ll,:)./Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                    else
                        Cyn_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Cyano(ll,:)./Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                        Dia_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Diatoms(ll,:)./Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                        Gre_lyr_prop(ii,:) = (Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Green(ll,:)./Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:))*100;
                    end
                end
                dat(jj,:,:) = Dia_lyr_prop;
                dat2(jj,:,:) = Cyn_lyr_prop;
                dat3(jj,:,:) = Gre_lyr_prop;
            end
            x = 1:365; % set x axis
            %clr = [0,1,1,[1,0,0],[0,1,0]%distinguishable_colors(3); % get some useful colours
            % plot uncertainty envelopes
            figure()
            for kk = 1:1:Sn+1 %:2:11 % Define which layers to plot here
                subplot(6,2,kk)
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),'r'); % calculate percentiles)
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat2(kk,:,:),N,365)),'c');
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat3(kk,:,:),N,365)),'g');
                hold on
                datetick('x','mmm')
                %xlabel('Date')
                if length(Sn2) == 1
                    formatSpec = 'Floating Solar Coverage: %d %%'; %, Depth: %d m, Tank: %d';
                    A1 = ll-1;
                    A2 = mm;
                    A3_1 = [0:10:100];
                    A3 = A3_1(kk);
                    str = sprintf(formatSpec,A3); %,A1,A2);
                    title(str)
                    % ylim([0 ceil(max(dat,[],'all'))]);
                    ylim([0 100]);
                    yticks([0 50 100]);
                else
                    formatSpec = 'Floating Solar Coverage: %d %%'; %'Depth: %d m, Tank: %d, Inflow Factor: %g';
                    A1 = ll-1;
                    A2 = mm;
                    A3 = Sn2(pp);
                    A4_1 = [0:25:100];
                    A4 = A4_1(kk);
                    str = sprintf(formatSpec,A4); %,A1,A2);
                    title(str);
%                     ylim([0 ceil(max(dat,[],'all'))]);
                    ylim([0 100]);
                    yticks([0 50 100]);
                end
               set(gca,'FontSize',12)
                % ylabel('Chlorophyll-a (µg L^{-1})')
            end
        end
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpeciesProp_Envel_%dm_Tank%d';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            str = sprintf(formatSpec,B1,B2,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpeciesProp_Envel_%dm_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B2,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end

%% Individual total chlorophyll plots
% Subplots
clear dat dat2 dat3
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for ll = [2] % Layers to plot
            for jj = 1:Sn+1
                for ii = 1:N
                    if length(Sn2)== 1
                        Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    else
                        Chl_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Chla_Sum(ll,:);
                    end
                end
                dat(jj,:,:) = Chl_lyr;
            end
            x = 1:365; % set x axis
            clr = distinguishable_colors(11); % get some useful colours
            % plot uncertainty envelopes
            figure()
            for kk = 1:1:Sn+1 %:2:11 % Define which layers to plot here
                subplot(6,2,kk)
                Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),clr(kk,:)); % calculate percentiles)
                hold on
                datetick('x','mmm')
                %xlabel('Date')
                if length(Sn2) == 1
                    formatSpec = 'Floating Solar Coverage: %d %%'; %%, Depth: %d m, Tank: %d';
                    A1 = ll-1;
                    A2 = mm;
                    A3_1 = [0:10:100];
                    A3 = A3_1(kk);
                    str = sprintf(formatSpec,A3);%,A1,A2);
                    title({str}) %'Total Chlorophyll-a',
                    ylim([0 ceil(max(dat,[],'all'))]);
                else
                    formatSpec = 'Floating Solar Coverage: %d %%'; %'Depth: %d m, Tank: %d, Inflow Factor: %g';
                    A1 = ll-1;
                    A2 = mm;
                    A3 = Sn2(pp);
                    A4_1 = [0:25:100];
                    A4 = A4_1(kk);
                    str = sprintf(formatSpec,A4); %,A1,A2);
                    title({str})
                    ylim([0 ceil(max(dat,[],'all'))]);
                end
                %ylabel('Chlorophyll-a (µg L^{-1})')
            end
        end
        % Save
        if length(Sn2) == 1
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlTotSub_Envel_%dm_Tank%d';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            str = sprintf(formatSpec,B1,B2,B3);
            saveas(gcf,[filepath,str]);
        else
            formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlTotSub_Envel_%dm_Tank%d_InfwF%g';
            B1 = ScenarioNum;
            B2 = ll-1;
            B3 = mm;
            B4 = pp;
            str = sprintf(formatSpec,B1,B2,B3,B4);
            saveas(gcf,[filepath,str]);
        end
    end
end


% %% Individual algal species plots
% 
% clear dat dat2 dat3
% for mm = 1:Tk
%     for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
%         for ll = [2] % Layers to plot
%             for jj = 1:Sn+1
%                 for ii = 1:N
%                     if length(Sn2)== 1
%                         Dia_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Diatoms(ll,:);
%                         Gre_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Green(ll,:);
%                         Cyn_lyr(ii,:) = Chlorophyll.tanks(mm).Reduction(jj).RunNum(ii).Cyano(ll,:);
%                     else
%                         Dia_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Diatoms(ll,:);
%                         Gre_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Green(ll,:);
%                         Cyn_lyr(ii,:) = Chlorophyll.tanks(mm).Scenario(pp).Reduction(jj).RunNum(ii).Cyano(ll,:);
%                     end
%                 end
%                 dat(jj,:,:) = Dia_lyr;
%                 dat2(jj,:,:) = Gre_lyr;
%                 dat3(jj,:,:) = Cyn_lyr;
%             end
%             x = 1:365; % set x axis
%             clr = distinguishable_colors(3); % get some useful colours
%             % plot uncertainty envelopes
%             for kk = 1:1:Sn+1 %:2:11 % Define which layers to plot here
%                 figure()
%                 Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat(kk,:,:),N,365)),clr(1,:)); % calculate percentiles)
%                 Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat2(kk,:,:),N,365)),clr(2,:));
%                 Envel_fill(SimDate,DeterminePercents(GLF,[0.000001,0.999999],reshape(dat3(kk,:,:),N,365)),clr(3,:));
%                 hold on
%                 datetick('x','mmm')
%                 xlabel('Date')
%                 if length(Sn2) == 1
%                     formatSpec = 'Floating Solar Coverage: %d %%, \n Depth: %d m, Tank: %d';
%                     A1 = ll-1;
%                     A2 = mm;
%                     A3_1 = [0:10:100];
%                     A3 = A3_1(kk);
%                     str = sprintf(formatSpec,A3,A1,A2);
%                     title({'Functional Group Chlorophyll-a', str})
%                     %ylim([0 ceil(max(dat,[],'all'))]);
%                 else
%                     formatSpec = 'Floating Solar Coverage: %d %%, Depth: %d m, Tank: %d, Inflow Factor: %g';
%                     A1 = ll-1;
%                     A2 = mm;
%                     A3 = Sn2(pp);
%                     A4_1 = [0:25:100];
%                     A4 = A4_1(kk);
%                     str = sprintf(formatSpec,A1,A2,A3);
%                     title({'Functional Group Chlorophyll-a', str})
%                     ylabel('Chlorophyll-a (µg L^{-1})')
%                 end
%                 % Save
%                 if length(Sn2) == 1
%                     formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpecInd_Envel_%dm_FPV%d_Tank%d';
%                     B1 = ScenarioNum;
%                     B2 = ll-1;
%                     B3 = mm;
%                     str = sprintf(formatSpec,B1,B2,A3,B3);
%                     saveas(gcf,[filepath,str]);
%                 else
%                     formatSpec = '\\QE2_Scenario_%d\\Figures\\ChlSpecies_Envel_%dm_Tank%d_InfwF%g';
%                     B1 = ScenarioNum;
%                     B2 = ll-1;
%                     B3 = mm;
%                     B4 = pp;
%                     str = sprintf(formatSpec,B1,B2,B3,B4);
%                     saveas(gcf,[filepath,str]);
%                 end
%             end
%         end
%     end
% end

%% Mixed depth
clear dat
for mm = 1:Tk
    for pp = 1:length(Sn2) % Number of variable scenarios e.g. inflow factor
        for jj = 1:Sn+1
            for ii = 1:N
                if length(Sn2)== 1
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(mm).Sim(jj).metric(ii,:);
                else
                    Mxd_lyr(ii,:) = Mixed_Depth.tanks(mm).Scenario(pp).Sim(jj).metric(ii,:);
                end
            end
            Mxd_lyr_norm = Mxd_lyr/trapz(17*ones(365,1)); % Normalise
    end
    end
end

end

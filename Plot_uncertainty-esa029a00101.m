close all
clear all
addpath('R:\page\PROGRAM_FILES\MATLAB\useful_m_files')

ChlTest = open('Chlorophyll.mat')

N_sim = 77


% get data
dat = ChlTest.Chlorophyll_1.tanks(1).Reduction(1).Chl_1m(:,:);

% get data
dat2 = ChlTest.Chlorophyll_1.tanks(1).Reduction(3).Chl_1m(:,:);

% get data
dat3 = ChlTest.Chlorophyll_1.tanks(1).Reduction(5).Chl_1m(:,:);

figure
plot(dat')

pcts = [0.025,0.975]% Set percentiles to look at
GLF = ones(N,1); % Needs weightings but set to one here
Uncert_percs = DeterminePercents(GLF,pcts,dat)
Uncert_percs2 = DeterminePercents(GLF,pcts,dat2)
Uncert_percs3 = DeterminePercents(GLF,pcts,dat3)

figure
plot(Uncert_percs,'b','linewidth',2)
hold on
plot(Uncert_percs2,'r','linewidth',2)
% Now plot envelope of results

figure

x=1:365;            %#initialize x array
y1=Uncert_percs(:,2)' ;         %#create first curve
y2=Uncert_percs(:,1)';         %#create second curve
plot(x, y1, 'LineWidth', 0.05,'Color',[0.7,0.7,0.7]);
hold on;
plot(x, y2, 'LineWidth', 0.05,'Color',[0.7,0.7,0.7]);
x2 = [x, fliplr(x)];
inBetween = [y1, fliplr(y2)];
fill(x2, inBetween, [0.7,0.7,0.7],'LineStyle','none');
hold on
y1=Uncert_percs2(:,2)' ;         %#create first curve
y2=Uncert_percs2(:,1)';         %#create second curve
plot(x, y1, 'LineWidth', 0.05,'Color',[0.1,0.7,0.7]);
hold on;
plot(x, y2, 'LineWidth', 0.05,'Color',[0.1,0.7,0.7]);
x2 = [x, fliplr(x)];
inBetween = [y1, fliplr(y2)];
fill(x2, inBetween, [0.1,0.7,0.7],'LineStyle','none');
hold on
y1=Uncert_percs3(:,2)' ;         %#create first curve
y2=Uncert_percs3(:,1)';         %#create second curve
plot(x, y1, 'LineWidth', 0.05,'Color',[0.1,0.5,0.9]);
hold on;
plot(x, y2, 'LineWidth', 0.05,'Color',[0.1,0.5,0.9]);
x2 = [x, fliplr(x)];
inBetween = [y1, fliplr(y2)];
fill(x2, inBetween, [0.1,0.5,0.9],'LineStyle','none');

% Look at distribution of esimates for max CHL timestep

TS = 190 % timestep to look at
figure
hist(dat(:,TS),15)
hold on
h2 = histogram(dat2(:,TS),15)
h2.FaceColor=[1,0,0];

P1 = DeterminePercents(GLF,[0.0,0.05,0.5,0.95,1],dat);
P2 = DeterminePercents(GLF,[0.0,0.05,0.5,0.95,1],dat2);

% Boxplots

figure
PLP(1,:) = P1(TS,:);
PLP(2,:) = P2(TS,:);
for kk = 1:2
    yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
    lentp = 1;
    % med(i),q1(i),q3(i),upadj(i),loadj(i)
    SPIRboxplot(PLP(kk,3),PLP(kk,2),PLP(kk,4),PLP(kk,5),PLP(kk,1),lentp,yy,0,'+',1,kk*10);
    t = axis;
    text(kk*10-1,t(3)+1,['Scenario ',num2str(kk)])
    hold on
    
%     plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
end
title('Chl distribution on day of max Chl')

%% Get data for scenarios of interest
for jj = 1:Sn+1
  dat(jj,:,:) = Chlorophyll.tanks(1).Reduction(jj).Chl_1m(:,:);  
end 
figure
for kk = 1:Sn+1
    P1 = DeterminePercents(GLF,[0.000001,0.1,0.5,0.9,0.9999],reshape(dat(kk,:,:),N,365)); % calculate percentiles
    PLP = mean(P1(TS-5:TS+5,:)); % Get mean of percentiles in the 11 day range
    yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
    lentp = 1;
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot(PLP(3),PLP(2),PLP(4),PLP(5),PLP(1),lentp,yy,0,'+',1,kk*10);
    text(kk*10-1,PLP(1)-1,['Scenario ',num2str(kk)])
    hold on
    % modify x-tick labels
%     plot(obs_dat(count,1),obs_dat(count,2),'k*','MarkerSize',20)% observed
end
title('Chl distribution on day of max Chl')



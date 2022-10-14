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
GLF = rand(N_sim,1); % Needs weightings but set to one here
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



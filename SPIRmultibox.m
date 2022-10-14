%% Plot multiple temperature scenarios on a box plot
Sn = 10;
yy = {[]}; % Need blank double brackets to be same as number of plotted boxes
lentp = 1;
G = 10:30:310;
x = [5 35 35 5];
y = [18.05 18.05 23.95 23.95];
patch(x,y,[1 1 0.96],'EdgeColor','none')
hold on
for jj = 3:2:Sn+1
    x = x+60;
y = [18.05 18.05 23.95 23.95];
patch(x,y,[1 1 0.96],'EdgeColor','none')
hold on
end
PLP = PLP_rec_8.tanks(1).WaterTemp_Max(1).PLP; % Scenario 1
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj),'b-',0.8);
    hold on
end
PLP = PLP_rec_5.tanks(1).WaterTemp_Max(1).PLP; % Scenario 2
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+10,'r-',0.8);
    hold on
end
PLP = PLP_rec_6.tanks(1).WaterTemp_Max(1).PLP; % Scenario 3
for jj = 1:Sn+1
    %          Median ,pecentile 1,percentile 3,max,min
    SPIRboxplot_clr(PLP(jj,3),PLP(jj,2),PLP(jj,4),PLP(jj,5),PLP(jj,1),lentp,yy,0,'+',1,G(jj)+20,'g-',0.8);
    hold on
end
title('Temperature distribution on day of maximum temperature')
xticks([5:30:335]) %([20:30:320])
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
xlabel('Floating solar array coverage (%)')
ylabel('Water Temperature (^{o}C)')
%ylim([18 24])
set(gca,'FontSize',12)
N = 10;

xmin = [ 
    [0.25, 1], 
    [0.25, 1], % Fraction of PAR in incoming solar radiation
    [0.7, 1],
    [-2, 1],
    [0.5, 1],
    [0.5, 1],
    [0.25, 1],
    [0.5, 1],
    [0.025, 1], % Non PAR light attenuation coefficient
    [0.025, 1]
    ]; % minimum values
xmax = [
    [1, 2], 
    [0.6, 2], % Fraction of PAR in incoming solar radiation
    [1.1, 2],
    [2, 2],
    [1.5, 2],
    [1.5, 2],
    [4, 2],
    [1.5, 2],
    [0.5, 2], % Non PAR light attenuation coefficient
    [0.5, 2]
    ]; % maximum values  

for ii = 1:length(xmin(1,:));
WindShelter(:,ii) = ((xmax(1,ii)-xmin(1,ii))*rand(N,1)+xmin(1,ii));
FracPar(:,ii) = ((xmax(2,ii)-xmin(2,ii))*rand(N,1)+xmin(2,ii)); % Fraction of PAR in incoming solar radiation
InflowFactor(:,ii) = ((xmax(3,ii)-xmin(3,ii))*rand(N,1)+xmin(3,ii));
InflowTemp(:,ii) = ((xmax(4,ii)-xmin(4,ii))*rand(N,1)+xmin(4,ii));
InflowTP(:,ii) = ((xmax(5,ii)-xmin(5,ii))*rand(N,1)+xmin(5,ii));
InflowNO3(:,ii) = ((xmax(6,ii)-xmin(6,ii))*rand(N,1)+xmin(6,ii));
%InflowNH4(:,ii) = ((xmax(7,ii)-xmin(7,ii))*rand(N,1)+xmin(7,ii)); % Do not vary
InflowSi(:,ii) = ((xmax(8,ii)-xmin(8,ii))*rand(N,1)+xmin(8,ii));
NonPAR_LA(:,ii) = ((xmax(9,ii)-xmin(9,ii))*rand(N,1)+xmin(9,ii)); % Non PAR light attenuation coefficient
PAR_LA(:,ii) = NonPAR_LA(:,ii); %((xmax(10,ii)-xmin(10,ii))*rand(N,1)+xmin(10,ii)); % PAR light attenuation coefficient

Xrec1(ii).Tank = [WindShelter(:,1), FracPar(:,1), InflowFactor(:,1), InflowTemp(:,1), InflowTP(:,1), InflowNO3(:,1),...
    InflowSi(:,1), NonPAR_LA(:,1), PAR_LA(:,1), Loss1, Loss2, Loss3, Growth1, Growth2, Growth3,...
    PARSat1, PARSat2, PARSat3, SetVel1, SetVel2, SetVel3];
end

for ii = 1:N
for jj = 1:length(xmin(1,:));
%********* Overwrite with sampled parameter values here.
lake_params{5,1}{1,jj} = WindShelter(ii,jj); % Wind "shelter/enhancement! factor 0.5
lake_params{10,1}{1,jj} = FracPar(ii,jj);
lake_params{14,1}{1,jj} = InflowFactor(ii,jj);
lake_params{15,1}{1,jj} = InflowTemp(ii,jj);
lake_params{17,1}{1,jj} = InflowTP(ii,jj);
lake_params{23,1}{1,jj} = InflowNO3(ii,jj);
%lake_params{24,1}{1,jj} = InflowNH4(ii,jj);
lake_params{35,1}{1,jj} = InflowSi(ii,jj);
lake_params{37,1}{1,jj} = NonPAR_LA(ii,jj);
lake_params{38,1}{1,jj} = PAR_LA(ii,jj);
end
end


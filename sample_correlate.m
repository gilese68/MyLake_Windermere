function [S1, S2, S3] = sample_correlate(Smin, Smax, N, CC)

R_Sample1 = rand(N,1);
R_Sample2 = rand(N,1);
R_Sample3 = rand(N,1);

S1 = (R_Sample1 * (Smax(:,1)-Smin(:,1))+Smin(:,1));

R_Sample4 = R_Sample1*CC+R_Sample2*(1-CC);
R_Sample5 = R_Sample1*CC+R_Sample3*(1-CC);

S2 = (R_Sample4 * (Smax(:,2)-Smin(:,2))+Smin(:,2)); % Algae 2
S3 = (R_Sample5 * (Smax(:,3)-Smin(:,3))+Smin(:,3)); % Algae 3
close all
clear all

x = 1:365

dat = load('WNBCHL.csv')

intchl = interp1(dat(:,1),dat(:,2),x)
intchl(1:7) = 1
intchl(315:end) = 4
plot(x,intchl)
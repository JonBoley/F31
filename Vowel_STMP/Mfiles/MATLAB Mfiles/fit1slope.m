function [Cfit,MSE,fit]=fit1slope(x,y)
% File: fit1slope.m
%
% File to fit 1 slope to a range of data
% M. Heinz 5/15/2002
%

%% Fit 1st line
Cfit=polyfit(x,y,1);
fit=polyval(Cfit,x);
MSE=mean((fit-y).^2);

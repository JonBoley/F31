function [Cfit,MSE,fit]=fit2slopes(x,y,elbowIND)
% File: fit2slopes.m
%
% File to fit 2 slopes to a range of data
%


%% Fit 1st line
Cfit{1}=polyfit(x(1:elbowIND),y(1:elbowIND),1);
fit{1}=polyval(Cfit{1},x(1:elbowIND));

if(elbowIND==length(x))
	Cfit{2}=[]; fit{2}=[];
	MSE=mean((fit{1}-y).^2);
else
	Cfit{2}=polyfit(x(elbowIND+1:end),y(elbowIND+1:end),1);
	fit{2}=polyval(Cfit{2},x(elbowIND+1:end));
	MSE=mean(([fit{1} fit{2}]-y).^2);
end
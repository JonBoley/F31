% File: setMonitor.m
%
% Date: 29Dec2005 (M. Heinz) 
% 
% sets a global variable depending on whether I am using the laptop screen
% or the docking station screen.  This will be used for placing figures 
% in appropriate pllaces

global MONITOR MONITORtext

MONITORpositions=get(0,'MonitorPositions');

LAPTOP_MONITORpositions=[1 1 1400 1050];
DOCKING_MONITORpositions=[1 1 1280 1024];

MONITORtext={'laptop','docking','unknown'};
if sum(MONITORpositions-LAPTOP_MONITORpositions)==0
	MONITOR='laptop';
elseif sum(MONITORpositions-DOCKING_MONITORpositions)==0
	MONITOR='docking';
else
	MONITOR='unknown';
end

disp(sprintf('Monitor set to: ''%s''',MONITOR))

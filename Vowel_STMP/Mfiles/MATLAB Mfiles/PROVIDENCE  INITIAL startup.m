% File: startup.m
% Startup file for M. Heinz
% Created: 2005 Dec 23
%
% Initial file, until I can get back to Cambridge to get original
% Setup for R03 now
%



global R03_dir
R03_dir='C:\Documents and Settings\Mike\My Documents\Work\R03';
MATLABfiles_dir='C:\Documents and Settings\Mike\My Documents\Work\MATLAB Mfiles';

eval(['cd ''' R03_dir filesep 'Data Analysis'''])
path(strcat(MATLABfiles_dir,filesep,'EDYfiles'),path)
path(strcat(MATLABfiles_dir,filesep,'AFfunctions'),path)
path(MATLABfiles_dir,path)
path(strcat(R03_dir,filesep,'Data Analysis',filesep,'NOHR Mfiles'),path)




global NOHR_dir NOHR_ExpList
NOHR_dir=R03_dir;
NOHR_ExpList={'MH-2004_09_02-ANnorm-NOHR','MH-2005_04_18-ANnorm','MH-2004_11_18-NOHRnorm','MH-2005_07_01-VCNnorm', ...
	'MH-2004_12_14-NOHRdeafCat','MH-2005_07_08-VCNnorm','MH-2004_07_08-ANnorm-NOHR','MH-2005_02_10-ANmodel-ARO', ...
	'MH-2005_07_13-ANdeafcat','MH-2004_08_02-AN-NOHR','MH-2005_03_28-ANnormal'};	

global SavedPICS SavedPICnums SavedPICSuse
global FeaturesText FormsAtHarmonicsText InvertPolarityText
FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
	
FormsAtHarmonicsText={'yes','no'};
InvertPolarityText={'no','yes'};

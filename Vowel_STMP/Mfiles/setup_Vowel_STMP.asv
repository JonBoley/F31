% FILE: setup_Vowel_STMP.m
%
% May 13 2009 - setup for Reiri&Jon
%
% sets all directories, paths etc for whichever machine analysis is being
% run on

global ExpList CALIBpics ROOT_dir STMP_dir InvertPolarityText
global NOHR_dir NOHR_ExpList
ROOT_dir='C:\Research\MATLAB\Vowel_STMP';
STMP_dir=ROOT_dir;   % For compatibility for old code
NOHR_dir=ROOT_dir;   % For compatibility for old code
eval(['cd ''' ROOT_dir filesep 'Mfiles'''])
path(strcat(ROOT_dir,filesep,'Mfiles',filesep,'NOHR Mfiles'),path)
path(strcat(ROOT_dir,filesep,'Mfiles',filesep,'MATLAB Mfiles'),path)
path(strcat(ROOT_dir,filesep,'Mfiles',filesep,'newMfiles'),path) %2009/06/19 added by Reiri Sono
path(strcat(ROOT_dir,filesep,'Mfiles',filesep,'STMPfunctions'),path)
path(strcat(ROOT_dir,filesep,'Mfiles',filesep,'SCCfunctions'),path)
path(strcat(ROOT_dir,filesep,'Mfiles'),path)

% global MATLABfiles_dirs
% 	path(strcat(MATLABfiles_dir,filesep,'NEL_TDT_Mfiles'),path)

ExpList={ ...
	'MH-2004_11_18-NOHRnorm', ...
	'MH-2004_12_14-NOHRdeafCat', ...
	'MH-2005_04_18-ANnorm', ...
	'MH-2006_11_16-ANnorm', ...	
    'MH-2007_04_13-ANnorm_wABR',...
    'MH-2005_03_28-ANnormal',...
    'MH-2005_07_13-ANdeafcat',...
    'JB-2011_01_18-AN-norm-exposed'};
STMP_ExpList=ExpList;  % For compatibility for old code
NOHR_ExpList=ExpList;  % For compatibility for old code
CALIBpics=[1 2 2 2 32];

global FeaturesText FormsAtHarmonicsText 
FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
FormsAtHarmonicsText={'yes','no'};
InvertPolarityText={'no','yes'};
set(0,'DefaultTextInterpreter','tex')

% File: startup.m
% Startup file for M. Heinz
% Created: 2005 Dec 23
%
% Modified from startup file from Cambridge 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INCLUDE={'dummy'};  % Needed for book -keeping
% INCLUDE{end+1}='Recruitment';
% INCLUDE{end+1}='CNexps';
% INCLUDE{end+1}='LoudnessJNDsAnal';
% INCLUDE{end+1}='R03';
% INCLUDE{end+1}='JNDsPaper';
% INCLUDE{end+1}='RecruitPaper';
% INCLUDE{end+1}='ARO2004';
% INCLUDE{end+1}='IHCON2004';
% INCLUDE{end+1}='NoiseCorr';
% INCLUDE{end+1}='NOHRstim';


%  INCLUDE{end+1}='R03anal';  % use 2.04
INCLUDE{end+1}='STMPanal';  % use 2.04
%INCLUDE{end+1}='HY2004data';  % for Heinz and Young 2004 DATA

% INCLUDE{end+1}='SLHS658';



% INCLUDE{end+1}='Minnesota';
% INCLUDE{end+1}='Phase';
% INCLUDE{end+1}='RLFpaper';
% INCLUDE{end+1}='oldLoudness';
% INCLUDE{end+1}='ARO2003';
% INCLUDE{end+1}='ISH2003talk';
% INCLUDE{end+1}='ISH2003paper';

if sum(strcmp('HY2004data',INCLUDE))
	eval(['cd ''C:\Documents and Settings\Mike\My Documents\Work\Research\Heinz and Young 2004 Data'''])
	
	disp('Not setting anything up for Heinz and Young 2004 data')
	break
	%% Don't do anything!!!
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Add General MATLAB Mfile Directories
global MATLABfiles_dir 
MATLABfiles_dir='C:\Documents and Settings\Mike\My Documents\Work\MATLAB Mfiles';

% Took out: June 03, 2007 - put whatever needed in MATLABfiles
path(strcat(MATLABfiles_dir,filesep,'EDYfiles'),path)
path(strcat(MATLABfiles_dir,filesep,'AFfunctions'),path)
path(MATLABfiles_dir,path)

% setMonitor  % used to determine which monitor is being used, so I can place figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% AN Recruitment work
%%% DO FIRST, if any recruitment work is set to be done
if sum(strcmp('Recruitment',INCLUDE))
   global Recruit_dir 
	
%    Recruit_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment';
	Recruit_dir='C:\Documents and Settings\Mike\My Documents\Work\Cambridge PC\mgheinz\Recruitment';
	eval(['cd ''' Recruit_dir ''''])

   global SEVimpExps MILDimpExps NORMExps ALLimpExps
   SEVimpExps={'122001','032002','032702','040302'};
   MILDimpExps={'071002','091102','100902','101602','103002','110602'};
   NORMExps={'012501','013001','040501','041701','102501','110801','112901','061202'};
   ALLimpExps={'122001','032002','032702','040302','062002','071002','072302','091102','100902','101602', ...
         '103002','110602'};
   
   path(path,strcat(Recruit_dir,filesep,'POPData2'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'MGHfunctions'))
   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'POPanal'))

   %   path(path,strcat(Recruit_dir,filesep,'Mfiles',filesep,'MGHfunctions',filesep,'Phase'))
end

if sum(strcmp('R03',INCLUDE))
   global R03_dir
   R03_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Grants\R03 - Apr 22, 2004\Mfiles';
   eval(['cd ''' R03_dir filesep ''''])
   path(R03_dir,path)
end
if sum(strcmp('RLFpaper',INCLUDE))
   global RLFpaper_dir
   RLFpaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RLFslopes';
   eval(['cd ''' RLFpaper_dir filesep 'Mfiles'''])
   path(strcat(RLFpaper_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('LoudnessJNDsAnal',INCLUDE))
   global LoudnessJNDsAnal_dir
   LoudnessJNDsAnal_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment\Mfiles\LoudnessJND_MFiles';
   path(LoudnessJNDsAnal_dir,path)
   eval(['cd ''' LoudnessJNDsAnal_dir ''''])
end

if sum(strcmp('Phase',INCLUDE))
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'Phase'),path)
   disp('IN Phase')
end
if sum(strcmp('JNDsPaper',INCLUDE))
   global JNDsPaper_dir
   JNDsPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\JNDs Paper';
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'LoudnessJND_Mfiles'),path)
end
if sum(strcmp('RecruitPaper',INCLUDE))
   %% NEED TO INCLUDE: 'Recruitment','LoudnessJNDsAnal','RecruitPaper';
   if ~(sum(strcmp(INCLUDE,'Recruitment'))&sum(strcmp(INCLUDE,'LoudnessJNDsAnal'))&sum(strcmp(INCLUDE,'RecruitPaper')))
      error('Fix INCLUDES for RecruitPaper!!!!')
   end
   
   global RecruitPaper_dir
   RecruitPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RecruitmentAN Paper\Sub2';
   %    RecruitPaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\RecruitmentAN Paper\Sub2\ExtraAnalysis\checkFlatHFloss';
   
   %    path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'LoudnessJND_Mfiles'),path)
   path(strcat(RecruitPaper_dir,filesep,'Mfiles'),path)
   eval(['cd ''' RecruitPaper_dir filesep 'Mfiles'''])
end
if sum(strcmp('ARO2003',INCLUDE))
   global ARO2003_dir 
   ARO2003_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\ARO2003';
   eval(['cd ''' ARO2003_dir filesep 'Mfiles'''])
   path(strcat(ARO2003_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('ARO2004',INCLUDE))
   global ARO2004_dir 
   ARO2004_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\ARO2004';
   eval(['cd ''' ARO2004_dir filesep 'Mfiles'''])
   path(strcat(ARO2004_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('IHCON2004',INCLUDE))
   global IHCON2004_dir 
   IHCON2004_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\IHCON2004';
   eval(['cd ''' IHCON2004_dir filesep 'Talk' filesep 'Mfiles'''])
   path(strcat(IHCON2004_dir,filesep,'Talk',filesep,'Mfiles'),path)
   global LoudnessJNDsAnal_dir
   LoudnessJNDsAnal_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Recruitment\Mfiles\LoudnessJND_MFiles';
   path(LoudnessJNDsAnal_dir,path)
end
if sum(strcmp('ISH2003talk',INCLUDE))
   global ISHtalk_dir
   ISHtalk_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\ISH2003\Talk';
   eval(['cd ''' ISHtalk_dir filesep 'Mfiles'''])
   path(strcat(ISHtalk_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('ISH2003paper',INCLUDE))
   global ISHpaper_dir 
   ISHpaper_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Papers\ISH2003\Paper';
   path(strcat(ISHpaper_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('oldLoudness',INCLUDE))
   oldLoudness_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\Presentations\LabMtg_012403';
   path(strcat(oldLoudness_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('NoiseCorr',INCLUDE))
   path(strcat(Recruit_dir,filesep,'Mfiles',filesep,'NoiseCorr'),path)
end

if sum(strcmp('Minnesota',INCLUDE))
   global Minn_dir
   Minn_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\MinnesotaBME\Interview 051203\Seminar';
   path(strcat(Minn_dir,filesep,'Mfiles'),path)
end
if sum(strcmp('CNexps',INCLUDE))
   global CNexps_dir
   CNexps_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps';
   eval(['cd ''' CNexps_dir filesep 'Data Analysis'''])
   path(strcat(CNexps_dir,filesep,'Data Analysis'),path)
   
   GE_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps\NEL folders\Users\GE';
   GEMH_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\CNexps\NEL folders\Users\GE_MH';
   path(GE_dir,path)
   path(GEMH_dir,path)
   
end
if sum(strcmp('NOHRstim',INCLUDE))
   global NOHR_dir
   NOHR_dir='C:\Documents and Settings\Michael G. Heinz\My Documents\mgheinz\NOHR';
   eval(['cd ''' NOHR_dir filesep 'Stimuli'''])
   path(strcat(NOHR_dir,filesep,'Stimuli'),path)
   path(strcat(NOHR_dir,filesep,'Stimuli',filesep,'Vowel synthesis'),path)

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
end
if sum(strcmp('R03anal',INCLUDE))
	%%% clean up NOHR/R03 names eventually, but left as is for now!!!!	
	global R03_dir
	R03_dir='C:\Documents and Settings\Mike\My Documents\Work\Research\R03 Experiments';	
	global NOHR_dir NOHR_ExpList NOHR_IMPExpList
	eval(['cd ''' R03_dir filesep 'Data Analysis'''])
% 	path(strcat(R03_dir,filesep,'Data Analysis',filesep,'NOHR Mfiles'),path)
	path(strcat(R03_dir,filesep,'Data Analysis',filesep,'STMP Mfiles'),path)
	
	NOHR_dir=R03_dir;  % FIX THIS EVENTUALLY!!!
	NOHR_ExpList={'MH-2004_07_08-ANnorm-NOHR','MH-2004_08_02-AN-NOHR','MH-2004_09_02-ANnorm-NOHR','MH-2004_11_18-NOHRnorm', ...
		'MH-2004_12_14-NOHRdeafCat','MH-2005_02_10-ANmodel-ARO','MH-2005_03_28-ANnormal','MH-2005_04_18-ANnorm', ...
		'MH-2005_07_01-VCNnorm','MH-2005_07_08-VCNnorm','MH-2005_07_13-ANdeafcat'};
	NOHR_IMPExpList={'MH-2004_12_14-NOHRdeafCat','MH-2005_07_13-ANdeafcat'};

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
	set(0,'DefaultTextInterpreter','tex')
end

if sum(strcmp('SLHS658',INCLUDE))
	% 	SLHS658_dir='C:\Documents and Settings\Mike\My Documents\Work\Courses\SLHS 658 - Spring 2006\In-class Demos\Week1';
	SLHS658_dir='C:\Documents and Settings\Mike\My Documents\Work\Courses\SLHS 658 - Spring 2006\SLHS 658 Toolbox';
	cd(SLHS658_dir)
	% 	ANmodel_dir='C:\Documents and Settings\Mike\My Documents\Work\Courses\SLHS 658 - Spring 2006\In-class Demos\Week1\ARLO';
	% 	path(ANmodel_dir,path)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp('STMPanal',INCLUDE))
	%%% clean up NOHR/R03 names eventually, but left as is for now!!!!	
	global STMP_dir
	STMP_dir='C:\Documents and Settings\Mike\My Documents\Work\Research\R03 Experiments';	
	global STMP_ExpList STMP_IMPExpList
	eval(['cd ''' STMP_dir filesep 'Data Analysis'''])
% 	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'NOHR Mfiles'),path)
	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'STMP Mfiles'),path)
	path(strcat(STMP_dir,filesep,'Data Analysis',filesep,'STMP Mfiles',filesep,'STMPfunctions'),path)
	
	STMP_ExpList={'MH-2004_07_08-ANnorm-NOHR','MH-2004_08_02-AN-NOHR','MH-2004_09_02-ANnorm-NOHR','MH-2004_11_18-NOHRnorm', ...
		'MH-2004_12_14-NOHRdeafCat','MH-2005_02_10-ANmodel-ARO','MH-2005_03_28-ANnormal','MH-2005_04_18-ANnorm', ...
		'MH-2005_07_01-VCNnorm','MH-2005_07_08-VCNnorm','MH-2005_07_13-ANdeafcat','MH-2006_10_19-ANnorm', ...
      'MH-2006_11_03-ANnorm','MH-2006_11_16-ANnorm','MH-2007_02_23-ANnorm','MH-2007_03_09-ANnorm','MH-2007_04_13-ANnorm_wABR', ...
      'MH-2007_04_20-ANnorm_wABRs','AN-2007_06_20-modelSTMP'};
   STMP_IMPExpList={'MH-2004_12_14-NOHRdeafCat','MH-2005_07_13-ANdeafcat'};

   global FeaturesText FormsAtHarmonicsText InvertPolarityText
   FeaturesText={'F0','T0','F1','T1','F2','T2','F3','T3','TN','NO'};
   FormsAtHarmonicsText={'yes','no'};
   InvertPolarityText={'no','yes'};
	set(0,'DefaultTextInterpreter','tex')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Write out what was INCLUDED
disp(sprintf('\n**********\nINCLUDING:\n**********'))
disp(strvcat(INCLUDE{2:end}))

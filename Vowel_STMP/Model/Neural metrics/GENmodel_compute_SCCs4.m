% File: testGEN_compute_SCCs.m
% M. Heinz
% May 21, 2008
% updated: July 29, 2008
%
% Runs CCCanal.m (general function) for any two conditions (CF, SR, stim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all

LOADexistingDATA=0;
DATAfilename='test_BOOT_BBNAA_Jon';  % test any two stim
% DATAfilename='GENmod_CHIM_TFS16_20reps';  % test any two stim

global ROOT_dir
global MATLABfiles_dir

% path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'ARLO'),path)
path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'Jackson Spike Generator FAST'),path)
path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'ZB0607_v2'),path)
STIM_dir=fullfile(ROOT_dir,'Data Analysis','STMP Mfiles','BBNstim');
path(STIM_dir,path);

DATAdir='C:\Documents and Settings\Mike\My Documents\Work\Research\R03 Experiments\Data Analysis\STMP Mfiles\testBBN_data';
eval(['cd ''' DATAdir ''''])

disp(sprintf('\n *** STARTING GENmodel_compute_SCCs ["%s"]',DATAfilename))

if ~LOADexistingDATA

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% SET PARAMETERS
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%% MODEL params - general (for both conditions)
	ANmodel_Fs_Hz=100000;
	Nreps=20;

	%% STIM params - general (for both conditions)
	OALevel_dBSPL=60;

	%% Condition parameters
	%%%%%%%%%%%%%
	% Condition A
	stim_FileName_A='BBN_A.wav';
% 	stim_FileName_A='Speech.wav';
	CF_A_kHz=.550;
	SR_A_sps=50;
	Cohc_A=1.0;
	Cihc_A=1.0;

	%%%%%%%%%%%%%
	% Condition B
	stim_FileName_B='BBN_A.wav';
% 	stim_FileName_B='CHSpeechFS-16.wav';
	CF_B_kHz=CF_A_kHz;
	SR_B_sps=SR_A_sps;
	Cohc_B=Cohc_A;
	Cihc_B=Cihc_A;

% 	%%%%%%%%%%%%%
% 	% Condition A
% 	stim_FileName_A='BBN_A.wav';
% 	CF_A_kHz=.536;
% 	SR_A_sps=50;
% 	Cohc_A=1.0;
% 	Cihc_A=1.0;
% 
% 	%%%%%%%%%%%%%
% 	% Condition B
% 	stim_FileName_B='BBN_A.wav';
% 	CF_B_kHz=.467;
% 	SR_B_sps=SR_A_sps;
% 	Cohc_B=Cohc_A;
% 	Cihc_B=Cihc_A;

	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Stimuli generation
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[stim_A FsA_Hz] = wavread(fullfile(STIM_dir,stim_FileName_A));
	[stim_B FsB_Hz] = wavread(fullfile(STIM_dir,stim_FileName_B));

	dur_sec=min([length(stim_A)/FsA_Hz length(stim_B)/FsB_Hz]);  % Set duration to minimum of two stimuli

	%%%%%%%%%%%
	%% RESAMPLE to ANmodel_Fs_Hz, with original as Fs_Hz
	disp(sprintf('... resampling stim_A (%s) from %.f Hz to %.f Hz',stim_FileName_A,FsA_Hz,ANmodel_Fs_Hz))
	dBSPL_A_before=20*log10(sqrt(mean(stim_A.^2))/(20e-6));
	sfreq=FsA_Hz;	   sfreqNEW=ANmodel_Fs_Hz;
	P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
	if(P/Q*sfreq~=sfreqNEW) disp('Integer sfreq conversion NOT exact'); end
	Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
	stim_A_model=resample(stim_A,P,Q,Nfir);
	dBSPL_A_after=20*log10(sqrt(mean(stim_A_model.^2))/(20e-6));
	if abs(dBSPL_A_before-dBSPL_A_after)>2;
		error(sprintf('RESAMPLING CHANGED stim_A dBSPL by %f dB',dBSPL_A_after-dBSPL_A_before))
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('... resampling stim_B (%s) from %.f Hz to %.f Hz',stim_FileName_B,FsB_Hz,ANmodel_Fs_Hz))
	dBSPL_B_before=20*log10(sqrt(mean(stim_B.^2))/(20e-6));
	sfreq=FsB_Hz;     Q=round(sfreq/10);  %Integers used to up sample
	if(P/Q*sfreq~=sfreqNEW) disp('Integer sfreq conversion NOT exact'); end
	stim_B_model=resample(stim_B,P,Q,Nfir);
	dBSPL_B_after=20*log10(sqrt(mean(stim_B_model.^2))/(20e-6));
	if abs(dBSPL_B_before-dBSPL_B_after)>2;
		error(sprintf('RESAMPLING CHANGED stim_B dBSPL by %f dB',dBSPL_B_after-dBSPL_B_before))
	end

	%% Scale stimuli to correct OALevel_dBSPL
	stim_A_model=stim_A_model*10^((OALevel_dBSPL-dBSPL_A_after)/20);
	stim_B_model=stim_B_model*10^((OALevel_dBSPL-dBSPL_B_after)/20);
	%% REFIT and WINDOWwavefile at ANmodel_Fs
	% Repeat or truncate waveform to fit requested stimulus duration:
	stim_A_model = refit_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
	stim_B_model = refit_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
	% Window waveform using linear rise/fall:
	stim_A_model = window_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
	stim_B_model = window_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Spike generation (from AN model for all 4 conditions: A+,A-,B+,B-)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	disp(sprintf('... Processing CONDITION:A \n    ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_A="%s"\n    ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
		CF_A_kHz,SR_A_sps,Cohc_A,Cihc_A,Nreps,stim_FileName_A,OALevel_dBSPL,dur_sec))
	%% stim_A_plus
	[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
		= zbcatmodel(stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,SR_A_sps);
	[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
	% Format spikes into NEL spikes format then cell array
	NELspikes=ANmodelSTs2nel(sptimes,Nreps);
	SpikeTrainsA_plus=nelSTs2cell(NELspikes);
	%% stim_A_minus
	[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
		= zbcatmodel(-stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_A,Cihc_A,SR_A_sps);
	[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
	% Format spikes into NEL spikes format then cell array
	NELspikes=ANmodelSTs2nel(sptimes,Nreps);
	SpikeTrainsA_minus=nelSTs2cell(NELspikes);
	
	disp(sprintf('... Processing CONDITION:B \n    ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_B="%s"\n    ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
		CF_B_kHz,SR_B_sps,Cohc_B,Cihc_B,Nreps,stim_FileName_B,OALevel_dBSPL,dur_sec))
	%% stim_B_plus
	[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
		= zbcatmodel(stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,SR_B_sps);
	[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
	% Format spikes into NEL spikes format then cell array
	NELspikes=ANmodelSTs2nel(sptimes,Nreps);
	SpikeTrainsB_plus=nelSTs2cell(NELspikes);
	%% stim_B_minus
	[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth] ...
		= zbcatmodel(-stim_B_model.',CF_B_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc_B,Cihc_B,SR_B_sps);
	[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
	% Format spikes into NEL spikes format then cell array
	NELspikes=ANmodelSTs2nel(sptimes,Nreps);
	SpikeTrainsB_minus=nelSTs2cell(NELspikes);
	clear sout timeout meout c1filterout c2filterout c1vihc c2vihc vihc psth NELspikes sptimes
	clear stim_A stim_B 

	disp(sprintf('SAVING SPIKE TRAINS: %s_STs.mat',DATAfilename))
	eval(['save ' DATAfilename '_STs.mat'])
else
	disp(sprintf('LOADING EXISTING SPIKE TRAINS: %s_STs.mat',DATAfilename))
	eval(['load ' DATAfilename '_STs.mat'])
	disp(sprintf('... CONDITION:A \n ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n   ... Stim_A="%s"\n   ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
		CF_A_kHz,SR_A_sps,Cohc_A,Cihc_A,Nreps,stim_FileName_A,OALevel_dBSPL,dur_sec))
	disp(sprintf('... CONDITION:B \n ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n   ... Stim_B="%s"\n   ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
		CF_B_kHz,SR_B_sps,Cohc_B,Cihc_B,Nreps,stim_FileName_B,OALevel_dBSPL,dur_sec))
end


%% Organize variables for CCCanal
SpikeTrains=cell(2); % {condition (1,2), polarity (plus,minus)}
SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

% specify params to be used
clear paramsIN
paramsIN.durA_msec=dur_sec*1000;
paramsIN.durB_msec=dur_sec*1000;
paramsIN.CF_A_Hz=CF_A_kHz*1000;
paramsIN.CF_B_Hz=CF_B_kHz*1000;
% Need to include CF_A, CF_B for more generality
paramsIN.MAXspikes=2000;

paramsIN.PSD_LHfreqs_Hz=[0 64; 0 50];  %additional freq ranges to compute CCCenv for
% Can specify which SCpeak and CCCenv to use 
% (DEFAULTS: SCpeak='adj: 0-CF';CCCenv='10-300, subBIAS') 
% paramsIN.SCpeak_TOUSE='IFFTraw';
% paramsIN.CCCenv_TOUSE='10-300, subBIAS';

which_CCCanal=4;
if which_CCCanal==0
	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_0(SpikeTrains,paramsIN,1);
elseif which_CCCanal==1
	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_1(SpikeTrains,paramsIN,1);
elseif which_CCCanal==2
	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_2(SpikeTrains,paramsIN,1);
elseif which_CCCanal==3
	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_3(SpikeTrains,paramsIN,1);
elseif which_CCCanal==4
	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_4(SpikeTrains,paramsIN,1);
end

%% SAVE CCCanal Data
disp(sprintf('SAVING CCCanal Data: %s_SCCs.mat',DATAfilename))
eval(['save ' DATAfilename '_SCCs.mat'])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OUTPUT various metrics of interest for AVGed value
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCpeakLIST={'raw       ','adj: 0-300','adj: 0-CF ','IFFTraw'};
SCpeakLIST={'IFFTraw'};
for i=1:length(SCpeakLIST)
	disp(sprintf('SCpeaks("%s") = A: %.2f;  B:%.2f;  AB:%.2f',SCpeakLIST{i}, ...
		SACSCCmetrics.SCpeaks_A(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics.SCpeaks_legend))), ...
		SACSCCmetrics.SCpeaks_B(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics.SCpeaks_legend))), ...
		SACSCCmetrics.SCpeaks_AB(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics.SCpeaks_legend))) ))
end

% CCCenvLIST={'rawSC','10-300, subBIAS','10-300, withBIAS','0-300, subBIAS','0-50, subBIAS','0-CF, subBIAS','IFFTrawSC'};
CCCenvLIST={'0-300, subBIAS'};
for i=1:length(CCCenvLIST)
	disp(sprintf('CCCenv("%s") = %.2f',CCCenvLIST{i}, ...
		SACSCCmetrics.CCCenvs(find(strcmp(CCCenvLIST{i},SACSCCmetrics.CCCenvs_legend)))))
end


% plot_CCCanal_0(SACSCCfunctions,SACSCCmetrics,paramsOUT)





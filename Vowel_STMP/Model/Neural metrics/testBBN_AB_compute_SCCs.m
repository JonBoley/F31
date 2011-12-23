% File: testBBN_AB_compute_SCCs.m
% M. Heinz 
% May 21, 2008
%
% Tests CCCanal.m (general function) for BBN stimuli 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

LOADexistingDATA=0;
% DATAfilename='testBBN_AB_2';  % test any two stim
DATAfilenames={'JON_Fig2_CF250Hz','JON_Fig2_CF500Hz','JON_Fig2_CF1000Hz','JON_Fig2_CF2000Hz','JON_Fig2_CF4000Hz'};  % test JON Fig. 2
CFs_kHz=[250 500 1000 2000 4000]/1000;

cccrep=1;

% global STMP_dir 
% global MATLABfiles_dir

% % path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'ARLO'),path)
% path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'Jackson Spike Generator FAST'),path)
% path(strcat(MATLABfiles_dir,filesep,'AN models',filesep,'ZB0607_v2'),path)
% STIM_dir=fullfile(STMP_dir,'Data Analysis','STMP Mfiles','BBNstim');
% path(STIM_dir,path);
% 
% DATAdir='C:\Documents and Settings\Mike\My Documents\Work\Research\R03 Experiments\Data Analysis\STMP Mfiles\testBBN_data';
% eval(['cd ''' DATAdir ''''])

STIM_dir='C:\School\Purdue\Research\MATLAB\Temporal Coding\Neural metrics\BBNstim';
DATAdir='C:\School\Purdue\Research\MATLAB\Temporal Coding\Neural metrics\testBBN_data';

%for i=1:length(DATAfilenames)
for i=1:1

    
    i=2;
    
	DATAfilename=DATAfilenames{i}
	CF_kHz=CFs_kHz(i)

    for k=1:cccrep

	if ~LOADexistingDATA

		%% MODEL params
		ANmodel_Fs_Hz=100000;
		% 	CF_kHz=1;
		SR_sps=50;
		Cohc=1.0;
		Cihc=1.0;
		Nreps=120;

		%% STIM params
		OALevel_dBSPL=60;
		dur_sec=1.7;
		stim_FileName_A='BBN_A.wav';
		stim_FileName_B='BBN_B.wav';

		disp(sprintf('... Processing CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n   ... Stim_A="%s"; Stim_B="%s"\n   ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
			CF_kHz,SR_sps,Cohc,Cihc,Nreps,stim_FileName_A,stim_FileName_B,OALevel_dBSPL,dur_sec))


		%% Stim generation
		[stim_A FsA_Hz] = wavread(fullfile(STIM_dir,stim_FileName_A));
		[stim_B FsB_Hz] = wavread(fullfile(STIM_dir,stim_FileName_B));

		%% RESAMPLE to ANmodel_Fs_Hz, with original as Fs_Hz
		dBSPL_A_before=20*log10(sqrt(mean(stim_A.^2))/(20e-6));
		sfreq=FsA_Hz;
		sfreqNEW=ANmodel_Fs_Hz;
		P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
		if(P/Q*sfreq~=sfreqNEW) disp('Integer sfreq conversion NOT exact'); end
		Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
		stim_A_model=resample(stim_A,P,Q,Nfir);
		dBSPL_A_after=20*log10(sqrt(mean(stim_A_model.^2))/(20e-6));
		% 					t_old=(1:length(stim_A))/Fs_A_Hz;
		% 					t_new=(1:length(stim_A_model))/ANmodel_Fs_Hz;
		if abs(dBSPL_A_before-dBSPL_A_after)>2;
			error(sprintf('RESAMPLING CHANGED stim_A dBSPL by %f dB',dBSPL_A_after-dBSPL_A_before))
		end
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		dBSPL_B_before=20*log10(sqrt(mean(stim_B.^2))/(20e-6));
		sfreq=FsB_Hz;
		sfreqNEW=ANmodel_Fs_Hz;
		P=round(sfreqNEW/10); Q=round(sfreq/10);  %Integers used to up sample
		if(P/Q*sfreq~=sfreqNEW) disp('Integer sfreq conversion NOT exact'); end
		Nfir=30;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time
		stim_B_model=resample(stim_B,P,Q,Nfir);
		dBSPL_B_after=20*log10(sqrt(mean(stim_B_model.^2))/(20e-6));
		% 					t_old=(1:length(stim_B))/Fs_B_Hz;
		% 					t_new=(1:length(stim_B_model))/ANmodel_Fs_Hz;
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

            
		%% Generate spikes from AN model for all 4 conditions
		%% stim_A_plus
		[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
			= zbcatmodel(stim_A_model.',CF_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
		[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
		% Format spikes into NEL spikes format then cell array
		NELspikes=ANmodelSTs2nel(sptimes,Nreps);
		SpikeTrainsA_plus=nelSTs2cell(NELspikes);
		%% stim_A_minus
		[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
			= zbcatmodel(-stim_A_model.',CF_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
		[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
		% Format spikes into NEL spikes format then cell array
		NELspikes=ANmodelSTs2nel(sptimes,Nreps);
		SpikeTrainsA_minus=nelSTs2cell(NELspikes);
		%% stim_B_plus
		[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
			= zbcatmodel(stim_B_model.',CF_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
		[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
		% Format spikes into NEL spikes format then cell array
		NELspikes=ANmodelSTs2nel(sptimes,Nreps);
		SpikeTrainsB_plus=nelSTs2cell(NELspikes);
		%% stim_B_minus
		[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,sout,psth500k] ...
			= zbcatmodel(-stim_B_model.',CF_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.05,Cohc,Cihc,SR_sps);
		[sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
		% Format spikes into NEL spikes format then cell array
		NELspikes=ANmodelSTs2nel(sptimes,Nreps);
		SpikeTrainsB_minus=nelSTs2cell(NELspikes);
		clear sout timeout meout c1filterout c2filterout c1vihc c2vihc vihc psth500k NELspikes sptimes
		clear stim_A stim_A_model stim_B stim_B_model


		disp(sprintf('SAVING SPIKE TRAINS: %s_STs.mat',DATAfilename))
		eval(['save ' DATAfilename '_STs.mat'])
	else
		disp(sprintf('LOADING EXISTING SPIKE TRAINS: %s_STs.mat',DATAfilename))
		eval(['load STs_' DATAfilename '_STs.mat'])
		disp(sprintf('... Condition CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n   ... Stim_A="%s"; Stim_B="%s"\n   ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
			CF_kHz,SR_sps,Cohc,Cihc,Nreps,stim_FileName_A,stim_FileName_B,OALevel_dBSPL,dur_sec))
	end


	%% Organize variables for CCCanal
	SpikeTrains=cell(2); % {condition (1,2), polarity (plus,minus)}
	SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

	% specify params to be used
	clear paramsIN
	paramsIN.durA_msec=dur_sec*1000;
	paramsIN.durB_msec=dur_sec*1000;
	paramsIN.CF_Hz=CF_kHz*1000;

	[SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal(SpikeTrains,paramsIN,0);  % don't plot inside CCCanals
    
    cccenv(k)=SACSCCmetrics.CCCenv;
	cccenv5(k)=SACSCCmetrics.CCCenv5;
    ccctfs(k)=SACSCCmetrics.CCCtfs;
    
    plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsOUT)
    
    end

    SACSCCmetrics.cccenv=cccenv;
    SACSCCmetrics.cccenv5=cccenv5;
    SACSCCmetrics.ccctfs=ccctfs;

	% CALL GENERAL SCC_CCC code (with option to do only 1 or 3 columns)
	% LATER: TODO - setup to only COMPUTE COL 1
	% EG for SCC[CFi,CFi] - just do 1st column


	%% SAVE CCCanal Data
	disp(sprintf('SAVING CCCanal Data: %s_SCCs.mat',DATAfilename))
	eval(['save ' DATAfilename '_SCCs.mat'])


end


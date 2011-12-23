% File: testBBN_AB_compute_SCCs.m
% M. Heinz / J.Swaminathan
% May 21, 2008
%
% Tests CCCanal.m (general function) for BBN stimuli 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

LOADexistingDATA=0;

DATAfilenames(1,:)={'ama_SNR_20','ama_SNR_10','ama_SNR_5','ama_SNR_0','ama_SNR_-2','ama_SNR_-5','ama_SNR_-7','ama_SNR_-10','ama_SNR_-15','ama_SNR_-20'};  % test JON Fig. 2
DATAfilenames(2,:)={'ana_SNR_20','ana_SNR_10','ana_SNR_5','ana_SNR_0','ana_SNR_-2','ana_SNR_-5','ana_SNR_-7','ana_SNR_-10','ana_SNR_-15','ana_SNR_-20'};  % test JON Fig. 2
DATAfilenames(3,:)={'anga_SNR_20','anga_SNR_10','anga_SNR_5','anga_SNR_0','anga_SNR_-2','anga_SNR_-5','anga_SNR_-7','anga_SNR_-10','anga_SNR_-15','anga_SNR_-20'};  % test JON Fig. 2
DATAfilenames(4,:)={'imi_SNR_20','imi_SNR_10','imi_SNR_5','imi_SNR_0','imi_SNR_-2','imi_SNR_-5','imi_SNR_-7','imi_SNR_-10','imi_SNR_-15','imi_SNR_-20'};  % test JON Fig. 2
DATAfilenames(5,:)={'ini_SNR_20','ini_SNR_10','ini_SNR_5','ini_SNR_0','ini_SNR_-2','ini_SNR_-5','ini_SNR_-7','ini_SNR_-10','ini_SNR_-15','ini_SNR_-20'};  % test JON Fig. 2
DATAfilenames(6,:)={'umu_SNR_20','umu_SNR_10','umu_SNR_5','umu_SNR_0','umu_SNR_-2','umu_SNR_-5','umu_SNR_-7','umu_SNR_-10','umu_SNR_-15','umu_SNR_-20'};  % test JON Fig. 2
DATAfilenames(7,:)={'unu_SNR_20','unu_SNR_10','unu_SNR_5','unu_SNR_0','unu_SNR_-2','unu_SNR_-5','unu_SNR_-7','unu_SNR_-10','unu_SNR_-15','unu_SNR_-20'};  % test JON Fig. 2

%CFs_kHz=[250 500 1000 2000 4000]/1000;

% DATAfilenames(1,:)={'ama_SNR_20','ama_SNR_10','ama_SNR_5','ama_SNR_0','ama_SNR_-5','ama_SNR_-10','ama_SNR_-15','ama_SNR_-20'};  % test JON Fig. 2

cccrep=5;

CFs_kHz=550/1000;

	
STIM_dir='C:\Program Files\MATLAB\R2007a\work\Neural metrics\BBNstim';
DATAdir='C:\Program Files\MATLAB\R2007a\work\Neural metrics\testBBN_data';

for j=1:length(DATAfilenames(:,1))

DATAfilename=DATAfilenames(j,:);
    
for i=1:length(DATAfilenames(1,:))
% for i=1:1

    
%     i=2;
    
 	CF_kHz=CFs_kHz;
    
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
		dur_sec=0.3;
		stim_FileName_A=strcat(DATAfilename{1},'.wav');
		stim_FileName_B=strcat(DATAfilename{i},'.wav');

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
		% 					t_n ew=(1:length(stim_A_model))/ANmodel_Fs_Hz;
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


		disp(sprintf('SAVING SPIKE TRAINS: %s_STs.mat',DATAfilename{i}))
		eval(['save ' DATAfilename{i} '_STs.mat'])
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
    
%     plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsOUT)
    
    end

    SACSCCmetrics.cccenv=cccenv;
    SACSCCmetrics.cccenv5=cccenv5;
    SACSCCmetrics.ccctfs=ccctfs;

    cccenvTP(i)=mean(cccenv); %#ok<AGROW>
    cccenv5TP(i)=mean(cccenv5); %#ok<AGROW>
    ccctfsTP(i)=mean(ccctfs); %#ok<AGROW,AGROW>
    

	%% SAVE CCCanal Data
	disp(sprintf('SAVING CCCanal Data: %s_SCCs.mat',DATAfilename{i}))
	eval(['save ' DATAfilename{i} '_SCCs.mat'])


end

%%%Plotting CCCENV and CCCTFS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[20 10 5 0 -2 -5 -7 -10 -15 -20];
% x=[20 10 5 0 -5 -10 -15 -20];
figure;
plot(x,cccenvTP,'x-k','LineWidth', 2.5); hold on
plot(x,cccenv5TP,'x-b','LineWidth', 2.5); 
plot(x,ccctfsTP,'x-r','LineWidth', 2.5); hold off
set(gca,'XDir','rev')
ylim([0 1])
title(DATAfilenames(j))
ylabel('CCC')
xlabel('SNR (dB)')
legend('cccenv','cccenv5','ccctfs');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


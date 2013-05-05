clear; close; home;
DATAfilename='normal_hsr';   % test any two stim
fibertype=3; %[3,2,1] = {high, med, low} spont rate

cf=2250;%F2
octdiffs=-0.5:0.1:0.5;

STIM_dir='C:\school\Purdue\Research\MATLAB\PhaseModulation';
DATAdir='C:\school\Purdue\Research\MATLAB\PhaseModulation\Data\DeltaBF\2250Hz';
ANmodel_Fs_Hz=100e3;
Nreps=25;
OALevel_dBSPL=70;

% Condition A
stim_FileName_A='heed_100.wav';
CF_A_kHz=cf/1000;
SR_A_sps=50;
Cohc_A=1.0;
Cihc_A=1.0;

% Condition B
stim_FileName_B='heed_100.wav';
CF_B_kHz=(cf.*2.^(octdiffs/2))./1000;
SR_B_sps=SR_A_sps;
Cohc_B=Cohc_A;
Cihc_B=Cihc_A;


%% Stimuli generation
[stim_A FsA_Hz] = wavread(fullfile(STIM_dir,stim_FileName_A));
[stim_B FsB_Hz] = wavread(fullfile(STIM_dir,stim_FileName_B));
dur_sec=min([length(stim_A)/FsA_Hz length(stim_B)/FsB_Hz]);  % Set duration to minimum of two stimuli

% RESAMPLE to ANmodel_Fs_Hz, with original as Fs_Hz
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
% REFIT and WINDOWwavefile at ANmodel_Fs
% Repeat or truncate waveform to fit requested stimulus duration:
stim_A_model = refit_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = refit_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);
% Window waveform using linear rise/fall:
stim_A_model = window_waveform(stim_A_model,ANmodel_Fs_Hz,dur_sec*1000);
stim_B_model = window_waveform(stim_B_model,ANmodel_Fs_Hz,dur_sec*1000);

%% Spike generation (from AN model for all 4 conditions: A+,A-,B+,B-)
tic;
for k=1:length(CF_B_kHz)
    
    disp(sprintf('... Processing CONDITION:A \n    ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_A="%s"\n    ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
        CF_A_kHz,SR_A_sps,Cohc_A,Cihc_A,Nreps,stim_FileName_A,OALevel_dBSPL,dur_sec))
    % stim_A_plus
    vihc = catmodel_IHC(stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,Cohc_A,Cihc_A);
    [sout,psth] = catmodel_Synapse(vihc,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,fibertype,0);
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsA_plus=nelSTs2cell(NELspikes);
    % stim_A_minus
    vihc = catmodel_IHC(-stim_A_model.',CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,Cohc_A,Cihc_A);
    [sout,psth] = catmodel_Synapse(vihc,CF_A_kHz*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,fibertype,0);
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsA_minus=nelSTs2cell(NELspikes);

    disp(sprintf('... Processing CONDITION:B \n    ... CF=%.3f kHz; SR=%.1f sps; Cohc=%.2f; Cihc=%.2f; Nreps=%d;\n    ... Stim_B="%s"\n    ... Overall Sound Level: %.1f dB SPL; Duration = %.2f sec', ...
        CF_B_kHz(k),SR_B_sps,Cohc_B,Cihc_B,Nreps,stim_FileName_B,OALevel_dBSPL,dur_sec))
    % stim_B_plus
    vihc = catmodel_IHC(stim_B_model.',CF_B_kHz(k)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,Cohc_B,Cihc_B);
    [sout,psth] = catmodel_Synapse(vihc,CF_B_kHz(k)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,fibertype,0);
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsB_plus=nelSTs2cell(NELspikes);
    % stim_B_minus
    vihc = catmodel_IHC(-stim_B_model.',CF_B_kHz(k)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,Cohc_B,Cihc_B);
    [sout,psth] = catmodel_Synapse(vihc,CF_B_kHz(k)*1000,1,1/ANmodel_Fs_Hz,dur_sec+0.200,fibertype,0);
    [sptimes nspikes]= SGfast([1/ANmodel_Fs_Hz, Nreps],sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps);
    SpikeTrainsB_minus=nelSTs2cell(NELspikes);
    
    clear sout timeout meout c1filterout c2filterout c1vihc c2vihc vihc psth NELspikes sptimes
    clear stim_A stim_B

    DATAfilename1=strcat(DATAfilename,'_',num2str(octdiffs(k)));
    disp(sprintf('SAVING SPIKE TRAINS: %s_STs.mat',DATAfilename1));
    eval(['cd ''' DATAdir '''']);
    eval(['save ' DATAfilename1 '_STs.mat']);

    % Organize variables for CCCanal
    SpikeTrains=cell(2); % {condition (1,2), polarity (plus,minus)}
    SpikeTrains={SpikeTrainsA_plus,SpikeTrainsA_minus;SpikeTrainsB_plus,SpikeTrainsB_minus};

    % specify params to be used
    clear paramsIN
    paramsIN.durA_msec=dur_sec*1000;
    paramsIN.durB_msec=dur_sec*1000;
    paramsIN.CF_A_Hz=CF_A_kHz*1000;
    paramsIN.CF_B_Hz=CF_B_kHz(k)*1000;
    % Need to include CF_A, CF_B for more generality
    paramsIN.MAXspikes=3600;

    paramsIN.PSD_LHfreqs_Hz=[0 64; 0 50; 0 100];  %additional freq ranges to compute CCCenv for
    % Can specify which SCpeak and CCCenv to use
    % (DEFAULTS: SCpeak='adj: 0-CF';CCCenv='10-300, subBIAS')
    % paramsIN.SCpeak_TOUSE='IFFTraw';
    % paramsIN.CCCenv_TOUSE='10-300, subBIAS';

    which_CCCanal=3;
    if which_CCCanal==0
        [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_0(SpikeTrains,paramsIN,1);
    elseif which_CCCanal==1
        [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_1(SpikeTrains,paramsIN,1);
    elseif which_CCCanal==2
        [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_2(SpikeTrains,paramsIN,1);
    elseif which_CCCanal==3
        [SACSCCfunctions,SACSCCmetrics,paramsOUT] = CCCanal_3(SpikeTrains,paramsIN,1);
    end

    % SAVE CCCanal Data
    DATAfilename2=strcat(DATAfilename,'_',num2str(octdiffs(k)));
    disp(sprintf('SAVING CCCanal Data: %s_SCCs.mat',DATAfilename2))
    eval(['cd ''' DATAdir ''''])
    eval(['save ' DATAfilename2 '_SCCs.mat'])

    % OUTPUT various metrics of interest for AVGed value
    SCpeakLIST={'raw       ','adj: 0-300','adj: 0-CF ','IFFTraw'};
    for i=1:length(SCpeakLIST)
        disp(sprintf('SCpeaks("%s") = A: %.2f;  B:%.2f;  AB:%.2f',SCpeakLIST{i}, ...
            SACSCCmetrics{end}.SCpeaks_A(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))), ...
            SACSCCmetrics{end}.SCpeaks_B(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))), ...
            SACSCCmetrics{end}.SCpeaks_AB(find(strcmp(deblank(SCpeakLIST{i}),SACSCCmetrics{end}.SCpeaks_legend))) ))
    end

    CCCenvLIST={'rawSC','10-300, subBIAS','10-300, withBIAS','0-300, subBIAS','0-50, subBIAS','0-CF, subBIAS','IFFTrawSC'};
    for i=1:length(CCCenvLIST)
        disp(sprintf('CCCenv("%s") = %.2f',CCCenvLIST{i}, ...
            SACSCCmetrics{end}.CCCenvs(find(strcmp(CCCenvLIST{i},SACSCCmetrics{end}.CCCenvs_legend)))))
    end

    %%%Getting all the sumcor and difcor coeffs for plotting%%%%%%%%
    scpeak_A(k)=SACSCCmetrics{6}.SCpeaks_A(4)-1;
    scpeak_B(k)=SACSCCmetrics{6}.SCpeaks_B(4)-1;
    scpeak_AB(k)=SACSCCmetrics{6}.SCpeaks_AB(4)-1;

    dcpeak_A(k)=SACSCCmetrics{6}.DCpeak_A;
    dcpeak_B(k)=SACSCCmetrics{6}.DCpeak_B;
    dcpeak_AB(k)=SACSCCmetrics{6}.DCpeak_AB;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%Getting all the CCC's for plotting%%%%%%%%%%%%%%%%%%%%%%%
    cccenv_0_300(k)= SACSCCmetrics{6}.CCCenvAVGs(7);
    cccenv_10_300(k)= SACSCCmetrics{6}.CCCenvAVGs(13);
    cccenv_0_cf(k)= SACSCCmetrics{6}.CCCenvAVGs(9);
    cccenv_0_50(k)= SACSCCmetrics{6}.CCCenvAVGs(1);
    cccenv_ifft(k)= SACSCCmetrics{6}.CCCenvAVGs(19);
    ccctfs(k) = SACSCCmetrics{6}.CCCtfs;
end
toc


%%%%To plot all the sumcor coeffs%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(octdiffs,scpeak_A,'--ks','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,scpeak_B,'--ko','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,scpeak_AB,'--kv','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
xlabel('Octdiffs'); ylabel('SC peak');
legend('SC A','SC B','SC AB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%to plot all the difcorr coeffs%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(octdiffs,dcpeak_A,'--ks','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,dcpeak_B,'--ko','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,dcpeak_AB,'--kv','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
xlabel('Octdiffs'); ylabel('DC peak');
legend('DC A','DC B','DC AB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%to plot all the ccc env coeffs%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(octdiffs,cccenv_0_300,'--kv','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,cccenv_10_300,'--ks','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,cccenv_0_cf,'--ko','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,cccenv_0_50,'--kx','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,cccenv_ifft,'--kh','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
plot(octdiffs,ccctfs,'--rv','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',10); hold on
xlabel('Octdiffs'); ylabel('CCC ENV and TFS');
legend('0-300','10-300','0-cf','0-50','ifft','tfs');
title(sprintf('CF=%.3f Hz; SR=%.1f sp/sec; OAL = %.f dB',cf,SR_A_sps,OALevel_dBSPL))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save cccenv_deltaBF_500_R1.mat cccenv_0_300
% save ccctfs_deltaBF_500_R1.mat ccctfs
% save scpeakA_deltaBF_500_R1.mat scpeak_A
% save scpeakB_deltaBF_500_R1.mat scpeak_B
% save scpeakAB_deltaBF_500_R1.mat scpeak_AB
% save dcpeakA_deltaBF_500_R1.mat dcpeak_A
% save dcpeakB_deltaBF_500_R1.mat dcpeak_B
% save dcpeakAB_deltaBF_500_R1.mat dcpeak_AB
% save dcdelay_deltaBF_500_R1.mat dc_delay
% save scdelay_deltaBF_500_R1.mat sc_delay

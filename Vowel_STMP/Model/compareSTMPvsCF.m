% compareSTMPvsCF
% this script compares actual nearby fibers and predicted fibers via STMP
%
% Written by Jon Boley

% Add all subdirectories to the path
addpath(genpath('C:\Research\MATLAB\Vowel_STMP\Model'));

featureNum = 2; % [1 2 3 ...] = [F1 F2 F3 ...]
Levels = 65;%[60 70 80 90];
SNRs = Inf;

RunCFs = 0; % otherwise, load(CFfilename)
RunSTMP = 0; % otherwise, load(STMPfilename)
CFfilename = 'STMPvsCF_CF_2012-11-27_131415.mat'; % default file to load
STMPfilename = 'STMPvsCF_STMP_2012-11-27_141514.mat'; % default file to load

midCF_kHz = 1.7*2.^(0); %model this and surrounding CFs

%% Initialize model parameters
numCFs = 100; % in addition to midCF_kHz
numCFs = numCFs + mod(numCFs,2); % make this an even number

spread=2; % total number of octaves (half in each direction)
lowCF_kHz=midCF_kHz*2^(-spread/2);
CF_kHz=(lowCF_kHz*2.^(0:spread/(numCFs-1):spread));

% Make the first fiber the one centered at the feature
CF_kHz = [midCF_kHz CF_kHz];
numCFs = length(CF_kHz);
deltaCF = log2(CF_kHz'/midCF_kHz);

fibertype=3; %[1,2,3]=[low,med,high] spontaneous rate
Cohc=1.0; %outer hair cell health
Cihc=1.0; %inner hair cell health
Nreps=120; %number of repetitions
ANmodel_Fs_Hz=100e3; %100kHz sampling

% Initialize analysis variables
analysis_init;

% Get stimulus
FeaturesText={'F1','F2','F3','F4','F5'};
formants = [500    1700    2500    3300   3750];
BWs      = [60      90     200     250    200];
F0=100;
Fs=33000;
dur = 0.5;

ORIGfreq=formants(featureNum);
ShiftFact=midCF_kHz*1e3/ORIGfreq;

[time, vowel] = dovowl(formants*ShiftFact,BWs*ShiftFact,F0,dur,Fs);
vowel=vowel-mean(vowel);  % Remove DC
vowel=vowel./max(abs(vowel))*0.99; % normalize
signal = [-1; zeros(dur*Fs-1,1)];%vowel;%sin(2*pi*midCF_kHz*1e3*(1:dur*Fs)/Fs)';

%% Run model (actual CFs)
if RunCFs
    for LevelIndex = 1:length(Levels)
        for SNRindex=1:length(SNRs)
            % Set the level
            OALevel_dBSPL = Levels(LevelIndex);

            % Adjust the sample rate & level for the model
            refit_stim;

            % Adjust stimulus length and apply window
            signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
            signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);

            fprintf('Calculating %d fibers...\n',numCFs)
            for FiberNumber=1:numCFs
                fprintf('.');

                % Run model
                vihc = catmodel_IHC(signal_model.',midCF_kHz*2.^deltaCF(FiberNumber)*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
                [sout,psth]=catmodel_Synapse(vihc,midCF_kHz*2.^deltaCF(FiberNumber)*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
                [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
                NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
                SpikeTrains_plus = nelSTs2cell(NELspikes);
                if length(SpikeTrains_plus)<Nreps
                    SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
                end
                Spikes_plus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_plus;

                vihc = catmodel_IHC(-signal_model.',midCF_kHz*2.^deltaCF(FiberNumber)*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
                [sout,psth]=catmodel_Synapse(vihc,midCF_kHz*2.^deltaCF(FiberNumber)*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
                SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
                [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
                NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
                SpikeTrains_minus = nelSTs2cell(NELspikes);
                if length(SpikeTrains_minus)<Nreps
                    SpikeTrains_minus = [SpikeTrains_minus cell(1,Nreps-length(SpikeTrains_minus))];
                end
                Spikes_minus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_minus;
                
                calc_xCF;

                SCCdelays_usec{FiberNumber,LevelIndex,SNRindex} = SACSCCfunctions.delays_usec;
                SCC{FiberNumber,LevelIndex,SNRindex} = SACSCCfunctions.SCC_AB_avg;
                
                Rate(FiberNumber,LevelIndex,SNRindex)=...
                    mean(mean(SACSCCmetrics.AvgRate_sps));

                LocalSpikeTimes = [];
                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    LocalSpikeTimes = [LocalSpikeTimes;...
                        repmat(i,length(Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}),1), ...
                        Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}];
                end

                PSTHbinWidth_sec = 40e-6;
                PSTH{FiberNumber,LevelIndex,SNRindex}=...
                    histc(LocalSpikeTimes(:,2),0:PSTHbinWidth_sec:(dur_sec+1.00));
                PSTH{FiberNumber,LevelIndex,SNRindex}=...
                    PSTH{FiberNumber,LevelIndex,SNRindex}/...
                    length(unique(LocalSpikeTimes(:,1)))/PSTHbinWidth_sec; %spikes/sec
                
                [PERhist_sps{FiberNumber,LevelIndex,SNRindex},...
                    binWidth_sec,NumDrivenSpikes] = ...
                    PERhist(LocalSpikeTimes,dur_sec*1000,F0);

                if ~mod(FiberNumber,25)
                    fprintf('\n');
                end
            end
            fprintf('\n');

            Rho_vDeltaCF(LevelIndex,SNRindex,:) = [NaN; Rho(2:end,LevelIndex,SNRindex)];
            CD_vDeltaCF(LevelIndex,SNRindex,:) = [0; CD(2:end,LevelIndex,SNRindex)];
        end
    end
    
    CFfilename = ['STMPvsCF_CF_' datestr(datenum(now),'yyyy-mm-dd_HHMMSS')];
    eval(sprintf('save %s.mat',CFfilename));
end

%% Run model (STMP version)
if RunSTMP
    for LevelIndex = 1:length(Levels)
        for SNRindex=1:length(SNRs)
            % Set the level
            OALevel_dBSPL = Levels(LevelIndex);% adjust sample rate (STMP)

            %keyboard;
            fprintf('Calculating %d fibers...\n',numCFs)
            for FiberNumber=1:numCFs
                fprintf('.');

                % adjust sample rate (STMP)
                STMPfactor_freq = 2^deltaCF(FiberNumber);
                STMPfactor_time = 1/STMPfactor_freq;

                % Adjust the sample rate & level for the model
                % refit the stimulus for the model
                dBSPL_before=20*log10(sqrt(mean(signal.^2))/(20e-6));
                sfreq=Fs/STMPfactor_freq;
                sfreqNEW=ANmodel_Fs_Hz;
                P=round(sfreqNEW/10);
                Q=round(sfreq/10);
%                 if (P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion not exact'); end
                Nfir=30;
                signal_model=resample(signal,P,Q,Nfir);
                dBSPL_after=20*log10(sqrt(mean(signal_model.^2))/(20e-6));
                if abs(dBSPL_before-dBSPL_after)>3
                    error('RESAMPLING CHANGED input by %f dB',dBSPL_before-dBSPL_after);
                end
                dur_sec_STMP(FiberNumber)=length(signal_model)/ANmodel_Fs_Hz;

                signal_model = signal_model*10^((OALevel_dBSPL-dBSPL_after)/20);

                % Adjust stimulus length and apply window
                signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)*1000);
                signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)*1000);

                % Run model
                vihc = catmodel_IHC(signal_model.',midCF_kHz*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)+1.00,Cohc,Cihc);
                [sout,psth]=catmodel_Synapse(vihc,midCF_kHz*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)+1.00,fibertype,1);
                SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
                [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
                NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
                SpikeTrains_plus = nelSTs2cell(NELspikes);
                if length(SpikeTrains_plus)<Nreps
                    SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
                end                
%                 SpikeTrains_plus=cellfun(@(x) x*STMPfactor_time,SpikeTrains_plus,'UniformOutput',false);
                Spikes_plus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_plus;


                vihc = catmodel_IHC(-signal_model.',midCF_kHz*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)+1.00,Cohc,Cihc);
                [sout,psth]=catmodel_Synapse(vihc,midCF_kHz*1e3,1,...
                    1/ANmodel_Fs_Hz,dur_sec_STMP(FiberNumber)+1.00,fibertype,1);
                SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
                [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
                NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
                SpikeTrains_minus = nelSTs2cell(NELspikes);
                if length(SpikeTrains_minus)<Nreps
                    SpikeTrains_minus = [SpikeTrains_minus cell(1,Nreps-length(SpikeTrains_minus))];
                end
%                 SpikeTrains_minus=cellfun(@(x) x*STMPfactor_time,SpikeTrains_minus,'UniformOutput',false);
                Spikes_minus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_minus;
                
                
                LocalSpikeTimes = [];
                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    LocalSpikeTimes = [LocalSpikeTimes;...
                        repmat(i,length(Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}),1), ...
                        Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}];
                end
                PSTHbinWidth_sec = 40e-6;
                PSTH_preSTMP{FiberNumber,LevelIndex,SNRindex}=...
                    histc(LocalSpikeTimes(:,2),0:PSTHbinWidth_sec:(dur_sec_STMP+1.00));
                PSTH_preSTMP{FiberNumber,LevelIndex,SNRindex}=...
                    PSTH_preSTMP{FiberNumber,LevelIndex,SNRindex}/...
                    length(unique(LocalSpikeTimes(:,1)))/PSTHbinWidth_sec; %spikes/sec

                if ~mod(FiberNumber,25)
                    fprintf('\n');
                end
            end
            fprintf('\n');
            
            %%%%%%%%%%%%
            % Determine neural delay
            LocalSpikeTimes = [];
            for FiberNumber=1:numCFs %pool across simulated CF
                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    LocalSpikeTimes = [LocalSpikeTimes;...
                        repmat(i,length(Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}),1), ...
                        Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}];
                end
            end
            NeuralDelay_sec = 0.0032 - 0.0024*log10(midCF_kHz);%estimateLatency(LocalSpikeTimes);
            fprintf('[Using neural delay = %1.1f - %1.1f - %1.1fms]\n',min(NeuralDelay_sec)*1000,NeuralDelay_sec(1)*1000,max(NeuralDelay_sec)*1000);
            %%%%%%%%%%%%
            
            fprintf('... Calculating STMP adjustments for %d fibers...\n',numCFs)
            for FiberNumber=1:numCFs
                fprintf('.');
                %%%%%%%%%%%%
                % adjust spike times (STMP)
                STMPfactor_freq = 2^deltaCF(FiberNumber);
                STMPfactor_time = 1/STMPfactor_freq;
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%
                % shift spikes (by NeuralDelay_sec(1)/STMPfactor_time ?)
%                 if STMPfactor_time>1
%                     spikeOffset_sec = NeuralDelay_sec(1)/STMPfactor_time;
%                 else
%                     spikeOffset_sec = NeuralDelay_sec(1)*STMPfactor_time;
%                 end
                spikeOffset_sec = NeuralDelay_sec(1)-NeuralDelay_sec(1)*STMPfactor_time;
                fprintf('CF = %1.1fkHz; offset = %1.1fms\n',CF_kHz(FiberNumber),spikeOffset_sec*1e3);
                %%%%%%%%%%%%%%%%%%%%%%%%%

                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    % Scale all spikes times by x(FeatureFreq/BF) AFTER
                    % compensating for neural delay
                    NDcompORIGspikes=Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i};% -  NeuralDelay_sec;
                    NEGATIVEinds=find(NDcompORIGspikes<0);
                    NDcompSCALEDspikes = NDcompORIGspikes*STMPfactor_time;
                    % Leave spikes before NeuralDelay AS IS since they must be spontaneous
                    % spikes
                    if ~isempty(NEGATIVEinds)
                        NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
                    end
                    temp = NDcompSCALEDspikes + spikeOffset_sec;% + NeuralDelay_sec;
                    Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}=temp(temp>=0);
                end
                for i=1:length(Spikes_minus{FiberNumber,LevelIndex,SNRindex})
                    % Scale all spikes times by x(FeatureFreq/BF) AFTER
                    % compensating for neural delay
                    NDcompORIGspikes=Spikes_minus{FiberNumber,LevelIndex,SNRindex}{i};% - NeuralDelay_sec;
                    NEGATIVEinds=find(NDcompORIGspikes<0);
                    NDcompSCALEDspikes = NDcompORIGspikes*STMPfactor_time;
                    % Leave spikes before NeuralDelay AS IS since they must be spontaneous
                    % spikes
                    if ~isempty(NEGATIVEinds)
                        NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
                    end
                    temp = NDcompSCALEDspikes + spikeOffset_sec;% + NeuralDelay_sec;
                    Spikes_minus{FiberNumber,LevelIndex,SNRindex}{i}=temp(temp>=0);
                end
                %%%%%%%%%%%%
                dur_sec = dur_sec_STMP(FiberNumber) * STMPfactor_time;
                
                LocalSpikeTimes = [];
                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    LocalSpikeTimes = [LocalSpikeTimes;...
                        repmat(i,length(Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}),1), ...
                        Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}];
                end

                PSTHbinWidth_sec = 40e-6;
                PSTH_STMP{FiberNumber,LevelIndex,SNRindex}=...
                    histc(LocalSpikeTimes(:,2),0:PSTHbinWidth_sec:(dur_sec+1.00));
                PSTH_STMP{FiberNumber,LevelIndex,SNRindex}=...
                    PSTH_STMP{FiberNumber,LevelIndex,SNRindex}/...
                    length(unique(LocalSpikeTimes(:,1)))/PSTHbinWidth_sec; %spikes/sec
                
                calc_xCF;

                SCCdelays_usec{FiberNumber,LevelIndex,SNRindex} = SACSCCfunctions.delays_usec;
                SCC{FiberNumber,LevelIndex,SNRindex} = SACSCCfunctions.SCC_AB_avg;

                Rate(FiberNumber,LevelIndex,SNRindex)=...
                    mean(mean(SACSCCmetrics.AvgRate_sps));

                [PERhist_sps{FiberNumber,LevelIndex,SNRindex},...
                    binWidth_sec,NumDrivenSpikes] = ...
                    PERhist(LocalSpikeTimes,dur_sec*1000,F0);
                
                if ~mod(FiberNumber,25)
                    fprintf('\n');
                end
            end
            fprintf('\n');
            

            Rho_vDeltaCF(LevelIndex,SNRindex,:) = [NaN; Rho(2:end,LevelIndex,SNRindex)];
            CD_vDeltaCF(LevelIndex,SNRindex,:) = [0; CD(2:end,LevelIndex,SNRindex)];
        end
    end

    STMPfilename = ['STMPvsCF_STMP_' datestr(datenum(now),'yyyy-mm-dd_HHMMSS')];
    eval(sprintf('save %s.mat',STMPfilename));
end


%% Plot PerHist across CF
% CFfilename = 'STMPvsCF_CF_2012-09-22_141530';
% STMPfilename = 'STMPvsCF_STMP_2012-09-22_141713'; % 1ms
if ~exist('CFfilename')
    CFfiles=dir('STMPvsCF_CF_*.mat');
    STMPfiles=dir('STMPvsCF_STMP_*.mat');
    CFfilename = CFfiles(end).name;
    STMPfilename = STMPfiles(end).name;
end
 
load(CFfilename,'-regexp','[^CFfilename^STMPfilename]');
for i=1:length(CF_kHz)
    PerHist2(:,i) = PERhist_sps{i,1,1};
%     PSTH2(:,i) = PSTH{i,1,1};
end

[temp,CFindices]=sort(CF_kHz);
if length(CF_kHz)>=10
    figure(1), subplot(211),
    imagesc(1e3*binWidth_sec*(1:256),midCF_kHz*2.^deltaCF(CFindices),PerHist2(:,CFindices)');
    % imagesc(0:PSTHbinWidth_sec:dur_sec,midCF_kHz*2.^deltaCF,PSTH');
    xlabel('Time (ms)'); ylabel('CF (kHz)');
    title('Actual CFs'); %colormap('bone');
    % PerHist_center(:,1) = PERhist_sps{ceil(end/2),1,1};
else
    figure(1), plotIndex=1; scaleSize=1e-4;
    for i=CFindices
%         subplot(length(CF_kHz),1,plotIndex), plot(1e3*binWidth_sec*(1:256),PerHist2(:,i),'b-'); hold on;
        hPerHist(1)=plot(1e3*binWidth_sec*(1:256),CF_kHz(i)+scaleSize*PerHist2(:,i),'b-'); hold on;
%         title(sprintf('Period Histogram (CF %1.2fkHz)',midCF_kHz*2.^deltaCF(i)));
        plotIndex=plotIndex+1;
    end
    xlabel('Time (ms)'); ylabel('CF (kHz)');
end

figure(2)
hRho(1)=plot(deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),'b'); hold on;
plot([0 0],[min(min(Rho_vDeltaCF)) max(max(Rho_vDeltaCF))],'k:');
text(0,max(max(Rho_vDeltaCF)),FeaturesText{featureNum},...
    'HorizontalAlignment','Center'); hold off;
title('\rho vs \DeltaCF');
xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
ylabel(['\rho (re ' FeaturesText{featureNum} ')']);

%%% LOAD STMPfilename
load(STMPfilename,'-regexp','[^CFfilename^STMPfilename^hRho]');
for i=1:length(CF_kHz)
    PerHist3(:,i) = PERhist_sps{i,1,1};
%     PSTH3(:,i) = PSTH{i,1,1};
end

[temp,CFindices]=sort(CF_kHz);
if length(CF_kHz)>=10
    figure(1), subplot(212),
    imagesc(1e3*binWidth_sec*(1:256),midCF_kHz*2.^deltaCF(CFindices),PerHist3(:,CFindices)');
    % imagesc(0:PSTHbinWidth_sec:dur_sec,midCF_kHz*2.^deltaCF,PSTH2');
    xlabel('Time (ms)'); ylabel('CF (kHz)');
    title('STMP'); %colormap('bone');
    % PerHist_center(:,2) = PERhist_sps{ceil(end/2),1,1};
else
    figure(1), plotIndex=1; scaleSize=1e-4;
    for i=CFindices
%         subplot(length(CF_kHz),1,plotIndex), plot(1e3*binWidth_sec*(1:256),PerHist3(:,i),'r-'); hold off;
        hPerHist(2)=plot(1e3*binWidth_sec*(1:256),CF_kHz(i)+scaleSize*PerHist3(:,i),'r-'); hold on;
        plotIndex=plotIndex+1;
    end
    legend(hPerHist,{'Actual CF','STMP'});
    xlabel('Time (ms)'); ylabel('CF (kHz)');
%     xlabel('Time(ms)'); ylabel('Rate(sps)');
end

figure(2), hold on;
hRho(2)=plot(deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),'g'); hold off;
legend(hRho,{'Actual CF','STMP'});

% figure, plot(1e3*binWidth_sec*(1:256),PerHist_center);

%% Plot PSTH
load(CFfilename,'-regexp','[^CFfilename^STMPfilename]');
load(STMPfilename,'-regexp','[^CFfilename^STMPfilename^PSTH]');
if length(CF_kHz)<10
    [B,IX]=sort(CF_kHz);
    maxTime = 0.030; %sec
    count=0;
    figure(3),
    for index=IX
        count=count+1;
        subplot(length(CF_kHz),1,count), hold on;
        hPSTH(1)=area(PSTHbinWidth_sec:PSTHbinWidth_sec:maxTime,PSTH{index}(1:round(maxTime/PSTHbinWidth_sec)),'FaceColor','b','EdgeColor','none');
%         hPSTH(1)=area(PSTHbinWidth_sec:PSTHbinWidth_sec:maxTime,PSTH_STMP{1}(1:round(maxTime/PSTHbinWidth_sec)),'FaceColor','k','EdgeColor','none');
        hPSTH(2)=area(PSTHbinWidth_sec:PSTHbinWidth_sec:maxTime,PSTH_preSTMP{index}(1:round(maxTime/PSTHbinWidth_sec)),'FaceColor','g','EdgeColor','none');
        hPSTH(3)=area(PSTHbinWidth_sec:PSTHbinWidth_sec:maxTime,PSTH_STMP{index}(1:round(maxTime/PSTHbinWidth_sec)),'FaceColor','r','EdgeColor','none');
        xlim([0 0.03]); ylim([0 3e3]); ylabel(sprintf('%1.2fkHz',CF_kHz(index)));
    end
    % legend('CF','BF (@ feature)','pre STMP','STMP');
    legend(hPSTH,{'CF','pre STMP','STMP'});
end

%% plot CD error vs. CF
if ~exist('CFfilename')
    CFfiles=dir('STMPvsCF_CF_*.mat');
    STMPfiles=dir('STMPvsCF_STMP_*.mat');
    CFfilename = CFfiles(end-3).name;
    STMPfilename = STMPfiles(end-3).name;
end

clear CD; 
colors='gbr'; %[CF,STMP,error]
for i=1:2
    switch i
        case 1
            load(CFfilename,'CF_kHz','deltaCF','SCCdelays_usec','SCC','CD_vDeltaCF');
        case 2
            load(STMPfilename,'CF_kHz','deltaCF','SCCdelays_usec','SCC','CD_vDeltaCF');
    end %switch i
    
    manualCD = 1; % manually pick CD?
    if manualCD
        CD(:,i) = calcCDreBF_manual(CF_kHz,squeeze(SCCdelays_usec),squeeze(SCC),squeeze(CD_vDeltaCF));
    else
        CD(:,i) = squeeze(CD_vDeltaCF);
    end %if manualCD
    
end %for i=1:2

T = 1e6/(CF_kHz(1)*1e3); % period (us)
figure(99), hold on; 
plot(deltaCF(2:end),CD(2:end,1)/T,[colors(1) '.:']);
plot(deltaCF(2:end),CD(2:end,2)/T,[colors(2) '.:']);
plot(deltaCF(2:end),(CD(2:end,1)-CD(2:end,2))/T,[colors(3) '.-']);
hold off;
legend('CF','STMP','error');
xlabel('CF difference (octaves)'); ylabel('Characteristic delay (CF cycles)');

figure(999), hold on;
N= length(deltaCF);
isort = [2:floor(N/2) 1 ceil(N/2):N];
isort2 = [2:floor(N/2) 1 ceil(N/2):N-1];
isort3 = [3:floor(N/2) 1 ceil(N/2):N];
% xAxis = (deltaCF(isort2)+deltaCF(isort3))/2; % CFdiff (octaves)
% yAxis1 = -diff(CD(isort,1)/T)./diff(deltaCF(isort)); % (-CF cycles/octave)
% yAxis2 = -diff(CD(isort,2)/T)./diff(deltaCF(isort)); % (-CF cycles/octave)
xAxis = (2.^deltaCF(isort2)*CF_kHz(1)+2.^deltaCF(isort3)*CF_kHz(1))/2; % CF (kHz)
yAxis1 = -diff(CD(isort,1)*1e-6)./diff(2.^deltaCF(isort)*CF_kHz(1)*1e3)*1e6; % (us)
yAxis2 = -diff(CD(isort,2)*1e-6)./diff(2.^deltaCF(isort)*CF_kHz(1)*1e3)*1e6; % (us)
plot(xAxis,yAxis1,[colors(1) '.:']);
plot(xAxis,yAxis2,[colors(2) '.:']);
hold off;
legend('CF','STMP');
% xlabel('CF difference (octaves)'); ylabel('Group Delay (-CF cycles/octave)');
xlabel('CF (kHz)'); ylabel('Group Delay (\mus)');
set(gca,'XScale','log','XTick',[0.25 0.5 1 1.7 2 2.5 4],...
    'XTickLabel',{'0.25','0.5','1.0','1.7','2.0','2.5','4.0'});
title(sprintf('CF = %1.2fkHz',CF_kHz(1)));


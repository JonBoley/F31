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
RunSTMP = 1; % otherwise, load(STMPfilename)
CFfilename = 'STMPvsCF_CF_2012-08-03_100213.mat'; % default file to load
STMPfilename = 'STMPvsCF_STMP_2012-08-03_100355.mat'; % default file to load

midCF_kHz = 1.7*2.^(0);%(-2:0.5:2); %model this and surrounding CFs

%% Initialize model parameters
% model_init; % default params
% initialize model parameters
% midCF_kHz=formants(featureNum)/1000; %center on this feature

numCFs=4; % in addition to midCF_kHz
if mod(numCFs,2)==1, numCFs=numCFs+1; end % make this an even number

spread=1; % total number of octaves (half in each direction)
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
signal = vowel;

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
                    histc(LocalSpikeTimes(:,2),0:PSTHbinWidth_sec:dur_sec);
                
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
                if abs(dBSPL_before-dBSPL_after)>2
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
            NeuralDelay_sec = estimateLatency(LocalSpikeTimes);
            fprintf('[Using neural delay = %1.1fms]\n',NeuralDelay_sec*1000);
            %%%%%%%%%%%%
            
            fprintf('... Calculating STMP adjustments for %d fibers...\n',numCFs)
            for FiberNumber=1:numCFs
                fprintf('.');
                %%%%%%%%%%%%
                % adjust spike times (STMP)
                STMPfactor_freq = 2^deltaCF(FiberNumber);
                STMPfactor_time = 1/STMPfactor_freq;
                for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
                    % Scale all spikes times by x(FeatureFreq/BF) AFTER
                    % compensating for neural delay
                    NDcompORIGspikes=Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}-NeuralDelay_sec;
                    NEGATIVEinds=find(NDcompORIGspikes<0);
                    NDcompSCALEDspikes = NDcompORIGspikes*STMPfactor_time;
                    % Leave spikes before NeuralDelay AS IS since they must be spontaneous
                    % spikes
                    if ~isempty(NEGATIVEinds)
                        NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
                    end
                    Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}=NDcompSCALEDspikes+NeuralDelay_sec;
                end
                for i=1:length(Spikes_minus{FiberNumber,LevelIndex,SNRindex})
                    % Scale all spikes times by x(FeatureFreq/BF) AFTER
                    % compensating for neural delay
                    NDcompORIGspikes=Spikes_minus{FiberNumber,LevelIndex,SNRindex}{i}-NeuralDelay_sec;
                    NEGATIVEinds=find(NDcompORIGspikes<0);
                    NDcompSCALEDspikes = NDcompORIGspikes*STMPfactor_time;
                    % Leave spikes before NeuralDelay AS IS since they must be spontaneous
                    % spikes
                    if ~isempty(NEGATIVEinds)
                        NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
                    end
                    Spikes_minus{FiberNumber,LevelIndex,SNRindex}{i}=NDcompSCALEDspikes+NeuralDelay_sec;
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
                PSTH{FiberNumber,LevelIndex,SNRindex}=...
                    histc(LocalSpikeTimes(:,2),0:PSTHbinWidth_sec:1+dur_sec);
                PSTH{FiberNumber,LevelIndex,SNRindex}=...
                    PSTH{FiberNumber,LevelIndex,SNRindex}/...
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
% CFfilename = 'STMPvsCF_CF_2012-07-27_085342'; % F1 +/- 1oct
% STMPfilename = 'STMPvsCF_STMP_2012-07-27_100701'; % F1 +/- 1oct
% CFfilename = 'STMPvsCF_CF_2012-07-26_213815'; % F2 +/- 1oct
% STMPfilename = 'STMPvsCF_STMP_2012-07-26_225211'; % F2 +/- 1oct

% CFfilename = 'STMPvsCF_CF_2012-07-31_101231.mat'; % F2(@1700) +/- 1oct
% STMPfilename = 'STMPvsCF_STMP_2012-07-31_112836.mat'; % F2(@1700) +/- 1oct
 
load(CFfilename,'-regexp','[^CFfilename^STMPfilename]');
for i=1:length(CF_kHz)
    PerHist2(:,i) = PERhist_sps{i,1,1};
%     PSTH2(:,i) = PSTH{i,1,1};
end

figure(1), subplot(211),
imagesc(1e3*binWidth_sec*(1:256),midCF_kHz*2.^deltaCF,PerHist2');
% imagesc(0:PSTHbinWidth_sec:dur_sec,midCF_kHz*2.^deltaCF,PSTH');
xlabel('Time (ms)'); ylabel('CF (kHz)');
title('Actual CFs'); %colormap('bone');
% PerHist_center(:,1) = PERhist_sps{ceil(end/2),1,1};

if length(CF_kHz)<10
    figure(2),
    for i=1:length(CF_kHz)
        subplot(length(CF_kHz),1,i), plot(1e3*binWidth_sec*(1:256),PerHist2(:,i),'b-'); hold on;
        title(sprintf('Period Histogram (CF %1.2fkHz)',midCF_kHz*2.^deltaCF(i)));
    end
end

figure(3)
hRho(1)=plot(deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),'b'); hold on;
plot([0 0],[min(min(Rho_vDeltaCF)) max(max(Rho_vDeltaCF))],'k:');
text(0,max(max(Rho_vDeltaCF)),FeaturesText{featureNum},...
    'HorizontalAlignment','Center'); hold off;
title('\rho vs \DeltaCF');
xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
ylabel(['\rho (re ' FeaturesText{featureNum} ')']);

%%% LOAD STMPfilename
load(STMPfilename,'-regexp','[^CFfilename^STMPfilename]');
for i=1:length(CF_kHz)
    PerHist3(:,i) = PERhist_sps{i,1,1};
%     PSTH3(:,i) = PSTH{i,1,1};
end

figure(1), subplot(212),
imagesc(1e3*binWidth_sec*(1:256),midCF_kHz*2.^deltaCF,PerHist2');
% imagesc(0:PSTHbinWidth_sec:dur_sec,midCF_kHz*2.^deltaCF,PSTH2');
xlabel('Time (ms)'); ylabel('CF (kHz)');
title('STMP'); %colormap('bone');
% PerHist_center(:,2) = PERhist_sps{ceil(end/2),1,1};

if length(CF_kHz)<10
    figure(2),
    for i=1:length(CF_kHz)
        subplot(length(CF_kHz),1,i), plot(1e3*binWidth_sec*(1:256),PerHist3(:,i),'r-'); hold off;
        legend('Actual CF','STMP');
    end
    xlabel('Time(ms)'); ylabel('Rate(sps)');
end

figure(3), hold on;
hRho(2)=plot(deltaCF(2:end),squeeze(Rho_vDeltaCF(:,:,2:end)),'g'); hold off;
legend(hRho,{'Actual CF','STMP'});

% figure, plot(1e3*binWidth_sec*(1:256),PerHist_center);


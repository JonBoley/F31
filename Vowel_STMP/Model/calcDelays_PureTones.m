% calcDelays_PureTones.m  (Measure f/CF delay)
% 1) plot PST histograms for tones at CF (CF ranging form 0.5-8kHz)
% 2) plot PST histograms for CF=1kHz (tones ranging +/- 1 octave)
% 3) plot PST histograms for CF=2kHz (tones ranging +/- 1 octave)
% 4) plot PST histograms for CF=4kHz (tones ranging +/- 1 octave)
% 5) do the same with impairment (does f/CF delay change?)

%% setup

% Add all subdirectories to the path
addpath(genpath(fileparts(mfilename('fullpath'))));

midCF_kHz = 0.5*2.^(0.01:0.01:1.99); 

numCFs = 100; % in addition to midCF_kHz
numCFs = numCFs + mod(numCFs,2); % make this an even number
spread=2; % total number of octaves (half in each direction)

fibertype=3; %[1,2,3]=[low,med,high] spontaneous rate
[Cohc,Cihc]=fitaudiogram(midCF_kHz*1e3,0*ones(size(midCF_kHz))); %flat loss
Nreps=120; %number of repetitions
ANmodel_Fs_Hz=100e3; %100kHz sampling
OALevel_dBSPL = 65;

Fs=33000;
dur = 0.25;

% Initialize analysis variables
SynOut=cell(length(midCF_kHz),numCFs);
Spikes_plus=cell(length(midCF_kHz),numCFs);

%figure setup
figure(1); hold on;
ylabel('Frequency (kHz)'); xlabel('Time (sec)');
set(gca,'YScale','log','YTick',[0.25 0.5 1 1.7 2 2.5 4 8],...
    'YTickLabel',{'0.25','0.5','1.0','1.7','2.0','2.5','4.0','8.0'});

%% run conditions
for CFindex = 1:length(midCF_kHz)
%     fprintf('CF = %1.1fkHz\n',midCF_kHz(CFindex));
    
    % if 1kHz, 2kHz, or 4kHz
    if (midCF_kHz(CFindex)==1) || (midCF_kHz(CFindex)==2) || (midCF_kHz(CFindex)==4)
        lowF_kHz=midCF_kHz(CFindex)*2^(-spread/2);
        F_kHz=(lowF_kHz*2.^(0:spread/(numCFs-1):spread));
    else
        F_kHz = midCF_kHz(CFindex);
    end

    for ToneIndex = 1:length(F_kHz)
%         fprintf('... f = %1.1fkHz\n',CF_kHz(ToneIndex));
        signal = sin(2*pi*F_kHz(ToneIndex)*1e3*(1:dur*Fs)/Fs)';
        
        % Adjust the sample rate & level for the model
        refit_stim;

        % Adjust stimulus length and apply window
        signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
        signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);

        % run model
        vihc = catmodel_IHC(signal_model.',midCF_kHz(CFindex)*1e3,1,...
            1/ANmodel_Fs_Hz,dur_sec+0.100,Cohc(CFindex),Cihc(CFindex));
        [sout,psth]=catmodel_Synapse(vihc,midCF_kHz(CFindex)*1e3,1,...
            1/ANmodel_Fs_Hz,dur_sec+0.100,fibertype,1);
        SynOut{CFindex,ToneIndex}=sout; % save the synapse output
        [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
        NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
        SpikeTrains_plus = nelSTs2cell(NELspikes);
        if length(SpikeTrains_plus)<Nreps
            SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
        end
        Spikes_plus{CFindex,ToneIndex} = SpikeTrains_plus;

        % calc PST histograms
        PSTHbinWidth_sec = 40e-6;
        LocalSpikeTimes = [];
        for i=1:length(Spikes_plus{CFindex,ToneIndex})
            LocalSpikeTimes = [LocalSpikeTimes;...
                repmat(i,length(Spikes_plus{CFindex,ToneIndex}{i}),1), ...
                Spikes_plus{CFindex,ToneIndex}{i}];
        end

        nreps = length(unique(LocalSpikeTimes(:,1)));
        startIndices = 1+[0; find(max(0,diff(LocalSpikeTimes(:,1)==1)))];
        endIndices = [find(max(0,diff(LocalSpikeTimes(:,1)==1))); length(LocalSpikeTimes(:,1))];
        mreps = length(startIndices); %may have reps of reps

        clear localPSTH;
        for i=1:mreps
            localPSTH(:,i) = histc(LocalSpikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(LocalSpikeTimes(:,2)));
            localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
        end

        % plot PST histograms
        switch midCF_kHz(CFindex)
            case 1, color='r'; plotYES=1;
            case 2, color='g'; plotYES=1;
            case 4, color='b'; plotYES=1;
            otherwise, color='k'; plotYES=0;
        end
        if plotYES
            sizeFactor = 50/max(max(localPSTH.^2));
            thresh = 500^2;
            for i=1:size(localPSTH,2)
                localData = localPSTH(1:750,i).^2;
                localData(localData<=thresh)=NaN;
                scatter((1:size(localData,1))*PSTHbinWidth_sec,...
                    F_kHz(ToneIndex)*ones(1,size(localData,1)),sizeFactor*localData,...
                    color);%,'MarkerFaceColor',color);
                drawnow;
            end
        end

    end


    F_kHz = [1];
    for ToneIndex = 1:length(F_kHz)
%         fprintf('... f = %1.1fkHz\n',CF_kHz(ToneIndex));
        signal = sin(2*pi*F_kHz(ToneIndex)*1e3*(1:dur*Fs)/Fs)';
        
        % Adjust the sample rate & level for the model
        refit_stim;

        % Adjust stimulus length and apply window
        signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
        signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);

        % run model
        vihc = catmodel_IHC(signal_model.',midCF_kHz(CFindex)*1e3,1,...
            1/ANmodel_Fs_Hz,dur_sec+0.100,Cohc(CFindex),Cihc(CFindex));
        [sout,psth]=catmodel_Synapse(vihc,midCF_kHz(CFindex)*1e3,1,...
            1/ANmodel_Fs_Hz,dur_sec+0.100,fibertype,1);
        SynOut{CFindex,ToneIndex}=sout; % save the synapse output
        [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
        NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
        SpikeTrains_plus = nelSTs2cell(NELspikes);
        if length(SpikeTrains_plus)<Nreps
            SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
        end
        Spikes_plus{CFindex,ToneIndex} = SpikeTrains_plus;

        % calc PST histograms
        PSTHbinWidth_sec = 40e-6;
        LocalSpikeTimes = [];
        for i=1:length(Spikes_plus{CFindex,ToneIndex})
            LocalSpikeTimes = [LocalSpikeTimes;...
                repmat(i,length(Spikes_plus{CFindex,ToneIndex}{i}),1), ...
                Spikes_plus{CFindex,ToneIndex}{i}];
        end

        nreps = length(unique(LocalSpikeTimes(:,1)));
        startIndices = 1+[0; find(max(0,diff(LocalSpikeTimes(:,1)==1)))];
        endIndices = [find(max(0,diff(LocalSpikeTimes(:,1)==1))); length(LocalSpikeTimes(:,1))];
        mreps = length(startIndices); %may have reps of reps

        clear localPSTH;
        for i=1:mreps
            localPSTH(:,i) = histc(LocalSpikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(LocalSpikeTimes(:,2)));
            localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
        end

        % plot PST histograms
        switch F_kHz(ToneIndex)
            case {1,2,4}, color='b'; plotYES=1;
%             case 2, color='b'; plotYES=1;
%             case 4, color='r'; plotYES=1;
            otherwise, color='k'; plotYES=0;
        end
        if plotYES
            sizeFactor = 50/max(max(localPSTH.^2));
            thresh = 500^2;
            for i=1:size(localPSTH,2)
                localData = localPSTH(1:750,i).^2;
                localData(localData<=thresh)=NaN;
                scatter((1:size(localData,1))*PSTHbinWidth_sec,...
                    midCF_kHz(CFindex)*ones(1,size(localData,1)),sizeFactor*localData,...
                    color);%,'MarkerFaceColor',color);
                drawnow;
            end
        end

    end

end


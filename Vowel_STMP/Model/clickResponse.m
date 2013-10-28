% clickResponse.m

%% setup

% Add all subdirectories to the path
addpath(genpath(fileparts(mfilename('fullpath'))));

baseCF_kHz = 1.0;

numCFs = 100; % in addition to midCF_kHz
numCFs = numCFs + mod(numCFs,2); % make this an even number
spread=6; % total number of octaves (half in each direction)

midCF_kHz = baseCF_kHz*2.^(-spread/2:spread/numCFs:spread/2);

fibertype=3; %[1,2,3]=[low,med,high] spontaneous rate
[Cohc,Cihc]=fitaudiogram(midCF_kHz*1e3,0*ones(size(midCF_kHz))); %flat loss
Nreps=120; %number of repetitions
ANmodel_Fs_Hz=100e3; %100kHz sampling

Fs=33000;
dur_sec = 0.500;
dur_plot = 0.010;

% 0.1ms condensation click at 50dBSPL, in accordance with:
% Wang, G. I., & Delgutte, B. (2012). Sensitivity of Cochlear Nucleus
% Neurons to Spatio-Temporal Changes in Auditory Nerve Activity. Journal of
% neurophysiology. doi:10.1152/jn.00160.2012
clickDur = 0.1e-3;
OALevel_dBSPL = 50;
signal_model=[ones(round(clickDur*ANmodel_Fs_Hz),1);...
    zeros(round(dur_sec*ANmodel_Fs_Hz)-round(clickDur*ANmodel_Fs_Hz),1)]; %click

% Initialize analysis variables
SynOut=cell(length(midCF_kHz),1);
Spikes_plus=cell(length(midCF_kHz),1);
PSTHbinWidth_sec = 40e-6;
bigPSTH = NaN*ones(length(midCF_kHz),round(dur_plot/PSTHbinWidth_sec));

% Greenwood Function for Chinchilla:
A = 163.5; k = 0.85; a = 2.1;
x = 0:0.001:1; %proportion of cochlear length
F_Hz = A * (10.^(a*x) - k);
dist = interp1(F_Hz,x,midCF_kHz*1e3);
fitLatency_sec = 0.005228*dist.^2 - 0.01203*dist + 0.008404; % max response
fitLatency_sec = -0.001121*dist.^2 - 0.00250*dist + 0.003349; % 1st peak

%figure setup
% figure(1); hold on;
ylabel('CF (kHz)'); xlabel('time (sec)');
set(gca,'YScale','log','YTick',[0.125 0.25 0.5 1 2 3 4 8]);
% xlim([0 0.010]); %ylim([1.5 4.5]);

latencies = NaN*ones(length(midCF_kHz),1);

%% run conditions
for CFindex = 1:length(midCF_kHz)
    % run model
    vihc = catmodel_IHC(signal_model.',midCF_kHz(CFindex)*1e3,1,...
        1/ANmodel_Fs_Hz,dur_sec+0.100,Cohc(CFindex),Cihc(CFindex));
    [sout,psth]=catmodel_Synapse(vihc,midCF_kHz(CFindex)*1e3,1,...
        1/ANmodel_Fs_Hz,dur_sec+0.100,fibertype,1);
    SynOut{CFindex}=sout; % save the synapse output
    [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
    SpikeTrains_plus = nelSTs2cell(NELspikes);
    if length(SpikeTrains_plus)<Nreps
        SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
    end
    if (midCF_kHz(CFindex)==baseCF_kHz)
        for i=1:length(SpikeTrains_plus)
            SpikeTrains_plus{i}=SpikeTrains_plus{i}-0.001;
        end
    end
    Spikes_plus{CFindex} = SpikeTrains_plus;

    % calc PST histograms
    LocalSpikeTimes = [];
    for i=1:length(Spikes_plus{CFindex})
        LocalSpikeTimes = [LocalSpikeTimes;...
            repmat(i,length(Spikes_plus{CFindex}{i}),1), ...
            Spikes_plus{CFindex}{i}];
    end
    
    temp = estimateClickLatency(LocalSpikeTimes);
    if ~isempty(temp), latencies(CFindex)=temp; end

    nreps = length(unique(LocalSpikeTimes(:,1)));
    startIndices = 1+[0; find(max(0,diff(LocalSpikeTimes(:,1)==1)))];
    endIndices = [find(max(0,diff(LocalSpikeTimes(:,1)==1))); length(LocalSpikeTimes(:,1))];
    mreps = length(startIndices); %may have reps of reps

    clear localPSTH;
    for i=1:mreps
        localPSTH(:,i) = histc(LocalSpikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(LocalSpikeTimes(:,2)));
        localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec

        bigPSTH(CFindex,:) = localPSTH(1:round(dur_plot/PSTHbinWidth_sec),i);
    end

end

% sizeFactor = 50/max(max(localPSTH.^2));
% thresh = 500^2;
% localData=bigPSTH;
% localData(localData<=thresh)=NaN;
% 
% scatter((1:size(localData,1))*PSTHbinWidth_sec,...
%     midCF_kHz(CFindex)*ones(1,size(localData,1)),sizeFactor*localData,...
%     'k');%,'MarkerFaceColor',color);

% surf((1:size(bigPSTH,2))*PSTHbinWidth_sec,midCF_kHz,bigPSTH,'EdgeColor','none');
% colormap(flipud(gray)); view([0 0 1]);
% % hold on; line(latencies,midCF_kHz,max(max(bigPSTH))*ones(length(midCF_kHz),1)); hold off;
% hold on; line(fitLatency_sec,midCF_kHz,max(max(bigPSTH))*ones(length(midCF_kHz),1)); hold off;

%% find new latencies...
% Greenwood Function for Chinchilla:
A = 163.5; k = 0.85; a = 2.1;
x = 0:0.001:1; %proportion of cochlear length
F_Hz = A * (10.^(a*x) - k);
dist = interp1(F_Hz,x,midCF_kHz*1e3);

figure, plot(dist,latencies); 
% manually picked points:
y_latency = [0.001	0.00116 0.0014  0.00168 0.00196 0.00232 0.00264 0.00300];
x_dist    = [0.8	0.7     0.6     0.5     0.4     0.3     0.2     0.1];
fitLatency_sec = interp1(x_dist,y_latency,dist);

figure,
surf((1:size(bigPSTH,2))*PSTHbinWidth_sec,midCF_kHz,bigPSTH,'EdgeColor','none');
colormap(flipud(gray)); view([0 0 1]);
% hold on; line(latencies,dist,max(max(bigPSTH))*ones(length(dist),1)); hold off;
hold on; line(fitLatency_sec,midCF_kHz,max(max(bigPSTH))*ones(length(midCF_kHz),1)); hold off;


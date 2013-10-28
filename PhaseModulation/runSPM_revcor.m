%% SPM - spatiotemporal phase modulation
% revcor to find best frequency, then apply SPM

% Add all subdirectories to the path
addpath(genpath(fileparts(mfilename('fullpath'))));

impaired = 0; % yes/no
Levels = 65;
SNRs = Inf;

% GenEh;              % create vowel
% actually, just use noise for now
dur = 0.5;
F0=100;
Fs=24414;
rng(1); % always use the same seed for random noise generator
vowel = randn(round(dur*Fs),1);
vowel=vowel./max(abs(vowel))*0.99; % normalize

midCF_kHz = 1; %center on this feature
model_init;         % Initialize model parameters
impairment_init;    % Set up impairment
analysis_init;      % Initialize analysis variables

% create filters
fPeak = midCF_kHz*1e3; %Hz
numFilters = 5;
numStages = 1;
GD=0.005*Fs; % group delay (samples)
k0=(GD-2)/(GD+2);
k2=1;
k1=-cos(2*pi*fPeak/Fs);
B=[k0 k1*(1+k0*k2) k2];
A=fliplr(B);
% figure(999), plot(1:Fs/2,grpdelay(B,A,Fs/2,Fs));
H=dfilt.df2t(B,A);
Hcas=dfilt.scalar;
for i=1:numFilters
    Hcas(i)=dfilt.scalar;
    for j=1:numStages
        Hcas(i)=dfilt.cascade(Hcas(i),H);
    end
end
for i=2:numFilters
    Hcas(i)=dfilt.cascade(Hcas(i-1),Hcas(i));
end

binSize_sec = 10e-6; % for PSTH

% get spikes
Rho_vDeltaCF = NaN*ones(numel(Hcas)+1,numel(Levels),numel(SNRs),numCFs);
CD_vDeltaCF = NaN*ones(numel(Hcas)+1,numel(Levels),numel(SNRs),numCFs);
for LevelIndex = 1:numel(Levels)
    for SNRindex=1:numel(SNRs)
        for FilterIndex=1:(numel(Hcas)+1)
            % Set the level
            OALevel_dBSPL = Levels(LevelIndex);
            fprintf('Level = %ddBSPL    ',Levels(LevelIndex))
            
            % Set the SNR
            if SNRs(SNRindex)~=Inf
                [signal,delta_dB]=AddWhiteNoise(vowel,SNRs(SNRindex));
            else
                signal=vowel;
            end
            fprintf('SNR = %ddB\n',SNRs(SNRindex))
            
            if FilterIndex>1
                signal = filter(Hcas(FilterIndex-1),signal);
            end
                
            
            % Adjust the sample rate & level for the model
            refit_stim;
            
            % Adjust stimulus length and apply window
            signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
            signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
            
            fprintf('Calculating %d fibers...\n',numCFs)
%             numCFs = 1; % for testing subset of fibers
            for FiberNumber=1:numCFs
                fprintf('.');
                
                % Set impairment
                if impaired
                    Cohc=interp1(Audiogram_freq,Cohc_impaired,CF_kHz(1)*1e3);
                    Cihc=interp1(Audiogram_freq,Cihc_impaired,CF_kHz(1)*1e3);
                end
                
                get_spikes;
                if numCFs>1
                    calc_xCF;
                end
                
                
                Spikes{LevelIndex,SNRindex,FilterIndex,FiberNumber} = ...
                    SpikeTrains_plus;
                
                PSTH{LevelIndex,SNRindex,FilterIndex,FiberNumber} = ...
                    getPSTH(SpikeTrains_plus,binSize_sec);
                
                Revcors{LevelIndex,SNRindex,FilterIndex,FiberNumber} = ...
                    revcor(SpikeTrains_plus, vowel, Fs);
                
            end
            fprintf('\n');
            
            if numCFs>1
                [Rho,CD] = calcCDreBF_manual(SACSCCfunctions,SACSCCmetrics,CF_kHz);
                
                Rho_vDeltaCF(FilterIndex,LevelIndex,SNRindex,:) = ...
                    [Rho(2:midCF_index), NaN, Rho(midCF_index+1:end)];
                CD_vDeltaCF(FilterIndex,LevelIndex,SNRindex,:) = ...
                    [CD(2:midCF_index), 0, CD(midCF_index+1:end)];
            end
            
            figure(1)
            plot(deltaCF,squeeze(CD_vDeltaCF),'.-');
            title('CD vs \DeltaCF');
            xlabel(['\DeltaCF (octaves re ' sprintf('%1.1fkHz)',fPeak/1e3)]);
            ylabel(['CD (re 0 octaves)']);
            drawnow;
        end
    end
end

%% lots of plots

if numCFs>1
    % plot CD vs. CF cycles
    figure,
    cycles = (squeeze(CD_vDeltaCF)/1e6)*fPeak;
    plot(deltaCF,cycles,'.-');
    title('CD vs \DeltaCF');
    xlabel(['\DeltaCF (octaves re ' sprintf('%1.1fkHz)',fPeak/1e3)]);
    ylabel(['CD (cycles re ' sprintf('%1.1fkHz)',fPeak/1e3)]);
    
    figure,
    plot(0:numFilters,interp1(deltaCF,cycles',-0.5),'o-');
    xlabel('Number of all-pass filters');
    ylabel(['CD @ -0.5 octaves (cycles re ' sprintf('%1.1fkHz)',fPeak/1e3)]);
end

% plot response of filters
maxPhi = NaN*ones(numel(Hcas),1);
figure, hold on;
lineWidth=numel(Hcas);
for i=1:numel(Hcas)
    [gd,w_gd] = grpdelay(Hcas(i)); % group delay
    [phid,w_phid] = phasedelay(Hcas(i)); % phase delay
    [phi,w_phi] = phasez(Hcas(i)); % phase
    maxPhi(i) = max(phid)/Fs*fPeak;
    plot(Fs*w_gd/(2*pi),gd/Fs*fPeak,'Linewidth',lineWidth,'Color','b');
    plot(Fs*w_phid/(2*pi),phid/Fs*fPeak,'Linewidth',lineWidth,'Color','g');
    plot(Fs*w_phi/(2*pi),phi/Fs*fPeak,'Linewidth',lineWidth,'Color','k');
    lineWidth = lineWidth - 1;
end
xlabel('Frequency (Hz)');
ylabel('Delay (CF cycles)');
set(gca,'XLim',[100 10e3],...
        'XTick',[100 250 500 1e3 2500 5e3 10e3],...
        'XScale','log');
legend('Group Delay','Phase Delay','Phase');

figure,
plot(mod(maxPhi,1),'o-');
xlabel('Number of all-pass filters');
ylabel(['max phase delay (F2 cycles)']);

if numCFs>1
    figure,
    plot([0; mod(maxPhi,1)],interp1(deltaCF,cycles',-0.5),'o');
    xlabel('filter delay (cycles)');
    ylabel('characteristic delay (cycles)');
end


%% plot revcor
doSmoothing = 1;
CFindex = 2;%numel(CF_kHz);
Npoints=numel(Revcors{1,1,1,1}); 
freqs = (1:Npoints)/Npoints*Fs/2;
lineWidth=numel(Hcas)+1;
figure, 
% plot revcor magnitude
for ii=1:(numel(Hcas)+1)
    if doSmoothing
        plot(freqs,20*log10(ThirdOctSmoothing(...
            abs(Revcors{1,1,ii,CFindex}),freqs)),'b','LineWidth',lineWidth);
    else
        plot(freqs,20*log10(abs(Revcors{1,1,ii,CFindex})),'b',...
            'LineWidth',lineWidth);
    end
    hold on;
    lineWidth = lineWidth - 1;
end
hold off;
ax1 = gca;
set(ax1,'XColor','b','YColor','b','XScale','log','xlim',[100 10e3]);
xlabel('Frequency (Hz)');
ylabel('Revcor Magnitude (dB)');
ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'XScale','log',...
           'xlim',[100 10e3],...
           'Color','none',...
           'XColor','k','YColor','k');
lineWidth=numel(Hcas)+1;
% plot revcor phase
for ii=1:(numel(Hcas)+1)
    if doSmoothing
        line(freqs,...
            ThirdOctSmoothing(unwrap(angle(Revcors{1,1,ii,CFindex}))/(2*pi),freqs),...
            'Color','k','Parent',ax2,'LineWidth',lineWidth);
    else
        line(freqs,unwrap(angle(Revcors{1,1,ii,CFindex}))/(2*pi),...
            'Color','k','Parent',ax2,'LineWidth',lineWidth);
    end
    % phase delay
%    line(freqs,...
%        ThirdOctSmoothing(unwrap(angle(Revcors{1,1,ii,1}))./...
%        (2*pi*freqs)'*1e3,freqs),...
%        'Color','k','Parent',ax2,'LineWidth',lineWidth);
    hold(ax2,'on');
    lineWidth = lineWidth - 1;
end
hold(ax2,'off');
ylabel('Revcor Phase (cycles)');


%% plot phase difference (re no filter) as a function of CF
% (phases at each CF)
doSmoothing = 0;
Npoints=numel(Revcors{1,1,1,1}); 
freqs = (1:Npoints)/Npoints*Fs/2;
[sortedCF_kHz,ix] = sort(CF_kHz);

refPhase = zeros(numCFs,1);
for FiberNumber=1:numCFs
    if doSmoothing
        phaseArray = ThirdOctSmoothing(...
            unwrap(angle(Revcors{1,1,1,FiberNumber})),freqs);
    else
        phaseArray = unwrap(angle(Revcors{1,1,1,FiberNumber}));
    end
    refPhase(FiberNumber) = phaseArray(find(freqs<=1e3*CF_kHz(FiberNumber),1,'last'));
end

lineWidth=numel(Hcas)+1;
figure, 
phaseDiff = zeros(numCFs,numel(Hcas));
for ii=1:numel(Hcas)
    for FiberNumber=1:numCFs
        if doSmoothing
            phaseArray = ThirdOctSmoothing(...
                unwrap(angle(Revcors{1,1,ii,FiberNumber})),freqs);
        else
            phaseArray = unwrap(angle(Revcors{1,1,ii,FiberNumber}));
        end
        phaseDiff(FiberNumber,ii) = ...
            phaseArray(find(freqs<=1e3*CF_kHz(FiberNumber),1,'last'));
    end   
    
    phaseToPlot = unwrap(phaseDiff(ix,ii))/(2*pi);
    phaseToPlot = phaseToPlot - phaseToPlot(ceil(end/2));
    plot(sortedCF_kHz,phaseToPlot,'LineWidth',lineWidth);
    hold on;
    lineWidth = lineWidth - 1;
end
title('phase difference (re no filter) at CF');
xlabel('characteristic frequency (kHz)');
ylabel('\Delta\phi (cycles)');


%% plot phase difference (re no filter) as a function of CF
% (phase at common frequency)
doSmoothing = 0;
Npoints=numel(Revcors{1,1,1,1}); 
freqs = (1:Npoints)/Npoints*Fs/2;
[sortedCF_kHz,ix] = sort(CF_kHz);

refPhase = zeros(numCFs,1);
for FiberNumber=1:numCFs
    if doSmoothing
        phaseArray = ThirdOctSmoothing(...
            unwrap(angle(Revcors{1,1,1,FiberNumber})),freqs);
    else
        phaseArray = unwrap(angle(Revcors{1,1,1,FiberNumber}));
    end
    refPhase(FiberNumber) = phaseArray(find(freqs<=1e3*CF_kHz(1),1,'last'));
end

lineWidth=numel(Hcas)+1;
figure, 
phaseDiff = zeros(numCFs,numel(Hcas));
for ii=1:numel(Hcas)
    for FiberNumber=1:numCFs
        if doSmoothing
            phaseArray = ThirdOctSmoothing(...
                unwrap(angle(Revcors{1,1,ii,FiberNumber})),freqs);
        else
            phaseArray = unwrap(angle(Revcors{1,1,ii,FiberNumber}));
        end
        phaseDiff(FiberNumber,ii) = ...
            phaseArray(find(freqs<=1e3*CF_kHz(1),1,'last'));
    end   
    
%     plot(sortedCF_kHz,unwrap(phaseDiff(ix,ii)),'LineWidth',lineWidth);
%     plot(sortedCF_kHz,unwrap(mod(unwrap(phaseDiff(ix,ii))/(2*pi),1),-0.5),'LineWidth',lineWidth);

    plotVals = mod(unwrap(phaseDiff(ix,ii))/(2*pi),1);
    jumpIndx = 1+find(diff(plotVals)<-0.5);
    plotVals(jumpIndx:end) = plotVals(jumpIndx:end) + 1;
    if min(plotVals)>0.5
        plotVals = plotVals - 0.5;
    end
    plot(sortedCF_kHz,plotVals,'LineWidth',lineWidth);
    hold on;
    lineWidth = lineWidth - 1;
end
title('phase difference (re no filter) at 1kHz');
xlabel('characteristic frequency (kHz)');
ylabel('\Delta\phi (cycles)');


%% plot coincidence detection
cWinLen_sec = 0.010e-3; % (0.1msec in Wang & Delgutte 2012)
L = 2; % min input spikes for output spike
w = 0.5; % CF range (octaves)
N = numCFs; % number of input fibers

maxBins = (dur_sec+1.00)/binSize_sec;

% get PSTH for each CF
CDprob = zeros(numel(Hcas)+1,maxBins);
for FilterIndex=1:(numel(Hcas)+1)
    PSTHmat = zeros(numCFs,maxBins);
    for FiberNumber=2:2:N
        % probability of a spike in each bin
        PSTHmat(FiberNumber,1:numel(PSTH{1,1,FilterIndex,FiberNumber})) = ...
            PSTH{1,1,FilterIndex,FiberNumber}/Nreps;
    end
    
    % calculate coincidence probability
    % sliding window of cWinLen_sec
    CDhist = filter(ones(round(cWinLen_sec/binSize_sec),1),1,PSTHmat);
    % count # of CFs with spikes
%     CFcount = sum(CDhist>0,1); 
    % threshold (any bin with <L spikes gets no output)
    CDprob(FilterIndex,:) = sum(CDhist,1);
    CDprob(FilterIndex,CDprob(FilterIndex,:)<(L/numCFs))=0;
%     CDprob(FilterIndex,:) = CDprob(FilterIndex,:)/numCFs;
end
figure, plot((1:maxBins)*binSize_sec,CDprob');
% 

totalCDprob = mean(CDprob(:,1:round(dur_sec/binSize_sec)),2)'/binSize_sec;
% figure, plot(0:numel(Hcas),totalCDprob,':o');
% % ylim([0 1]);
% xlabel('number of all-pass filters');
% ylabel('coincidence spikes/sec');


% figure, 
% for FilterIndex=1:(numel(Hcas)+1)
%     subplot((numel(Hcas)+1),1,FilterIndex)
%     temp = abs(fft(CDprob(FilterIndex,:)));
%     temp = temp(2:(end/2));
%     semilogx(temp);
% end
%

% period histogram of CD cell
ignoreDur = 20e-3;
perBins = round(1/F0/binSize_sec);
numPer = round((dur_sec-ignoreDur)/binSize_sec/perBins);
CDperHist = NaN*ones(1,perBins);
for FilterIndex=1:(numel(Hcas)+1)
    temp = CDprob(FilterIndex,1+round(ignoreDur/binSize_sec):round(dur_sec/binSize_sec));
    CDperHist(FilterIndex,:) = sum(reshape(temp,numPer,perBins),1);
end
figure,plot((0:(1/F0/binSize_sec))*binSize_sec,CDperHist);

%% FFT of CD cell PSTH
CDfft = NaN*ones(1,round((dur_sec-ignoreDur)/binSize_sec)/2);
for FilterIndex=1:(numel(Hcas)+1)
    temp = CDprob(FilterIndex,1+round(ignoreDur/binSize_sec):round(dur_sec/binSize_sec));
    temp = abs(fft(temp));
    temp = temp(1:end/2);
    CDfft(FilterIndex,:) = temp;
end
figure,plot(CDfft(:,2:end)');



%% more coincidence detection
cWinLen_sec = 0.1e-3; % 0.1msec (as in Wang & Delgutte 2012)
L = 2; % min input spikes for output spike
w = 0.5; % CF range (octaves)
N = numCFs; % number of input fibers

cWinDiv = 10; % divider for PSTH
maxBins = (dur_sec+1.00)/(cWinLen_sec/cWinDiv);

CDinputs = zeros(numel(Hcas)+1,maxBins);
CDoutputs = zeros(size(CDinputs));
for FilterIndex=1:(numel(Hcas)+1)
    for repNum=1:Nreps
        spikes = cell(1,numCFs);
        for FiberNumber=1:(round(numCFs/N)):numCFs
            spikes{FiberNumber} = Spikes{1,1,FilterIndex,FiberNumber}{repNum};
        end
        
        % inputPSTH = numFibers firing per bin
        % (assumes binSize << refractory period)
        inputPSTH = getPSTH(spikes,cWinLen_sec/cWinDiv)';
        CDhist = zeros(1,size(CDinputs,2));
        CDhist(1:numel(inputPSTH)) = ...
            filter(ones(cWinDiv,1),1,inputPSTH);
        CDinputs(FilterIndex,:) = CDinputs(FilterIndex,:) + CDhist;
        
        % threshold (any bin with <L spikes gets no output)
        CDoutputs(FilterIndex,CDinputs(FilterIndex,:)>=L) = ...
            1 + CDoutputs(FilterIndex,CDinputs(FilterIndex,:)>=L);
    end
    
end
% figure, plot(CDoutputs'/Nreps/binSize_sec);

totalCDprob = mean(CDoutputs(:,1:round(dur_sec/(cWinLen_sec/cWinDiv))),2)'/Nreps;
figure, plot(0:numel(Hcas),totalCDprob,':o');
ylim([0 1]);
xlabel('number of all-pass filters');
ylabel('coincidence probability');


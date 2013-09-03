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
            numCFs = 1; % for testing subset of fibers
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
Npoints=numel(Revcors{1,1,1,1}); 
freqs = (1:Npoints)/Npoints*Fs/2;
lineWidth=numel(Hcas)+1;
figure, 
% plot revcor magnitude
for ii=1:(numel(Hcas)+1)
    if doSmoothing
        plot(freqs,20*log10(ThirdOctSmoothing(...
            abs(Revcors{1,1,ii,1}),freqs)),'b','LineWidth',lineWidth);
    else
        plot(freqs,20*log10(abs(Revcors{1,1,ii,1})),'b',...
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
            ThirdOctSmoothing(unwrap(angle(Revcors{1,1,ii,1}))/(2*pi),freqs),...
            'Color','k','Parent',ax2,'LineWidth',lineWidth);
    else
        line(freqs,unwrap(angle(Revcors{1,1,ii,1}))/(2*pi),...
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


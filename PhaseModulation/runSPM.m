%% SPM - spatiotemporal phase modulation

% Add all subdirectories to the path
addpath(genpath('C:\Research\MATLAB\PhaseMod'));

impaired = 1; % yes/no
featureNum = 2; % [1 2 3 ...] = [F1 F2 F3 ...]
Levels = 65;
SNRs = Inf;

GenEh;              % create vowel
model_init;         % Initialize model parameters
impairment_init;    % Set up impairment
analysis_init;      % Initialize analysis variables

% create filters
fPeak = formants(featureNum); %Hz
numFilters = 5;
numStages = 2;
GD=0.005*Fs; % group delay (samples)
k0=(GD-2)/(GD+2);
k2=1;
% note that phase delay (~0.2ms/section) is independent of group delay
% also note that max phase delay may not be at center freq
k1=-cos(2*pi*fPeak/Fs);
B=[k0 k1*(1+k0*k2) k2];
A=fliplr(B);
% figure(999), plot(1:Fs/2,grpdelay(B,A,Fs/2,Fs));
H=dfilt.df2t(B,A);
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
            for FiberNumber=1:numCFs
                fprintf('.');
                
                % Set impairment
                if impaired
                    Cohc=interp1(Audiogram_freq,Cohc_impaired,CF_kHz(1)*1e3);
                    Cihc=interp1(Audiogram_freq,Cihc_impaired,CF_kHz(1)*1e3);
                end
                
                get_spikes;
                calc_xCF;
                blah=1;
                
            end
            fprintf('\n');
            
            [Rho,CD] = calcCDreBF_manual(SACSCCfunctions,SACSCCmetrics,CF_kHz);
            
            Rho_vDeltaCF(FilterIndex,LevelIndex,SNRindex,:) = ...
                [Rho(2:midCF_index), NaN, Rho(midCF_index+1:end)];
            CD_vDeltaCF(FilterIndex,LevelIndex,SNRindex,:) = ...
                [CD(2:midCF_index), 0, CD(midCF_index+1:end)];
            
            figure(1)
            plot(deltaCF,squeeze(CD_vDeltaCF),'.-');
            title('CD vs \DeltaCF');
            xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
            ylabel(['CD (re ' FeaturesText{featureNum} ')']);
            drawnow;
        end
    end
end

% plot CD vs. CF cycles
figure,
cycles = (squeeze(CD_vDeltaCF)/1e6)*formants(featureNum);
plot(deltaCF,cycles,'.-');
title('CD vs \DeltaCF');
xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
ylabel(['CD (cycles re ' FeaturesText{featureNum} ')']);

figure, 
plot(0:numFilters,interp1(deltaCF,cycles',-0.5),'o-');
xlabel('Number of all-pass filters');
ylabel(['CD @ -0.5 octaves (cycles re ' FeaturesText{featureNum} ')']);


% plot response of filters
maxPhi = NaN*ones(numel(Hcas),1);
figure, hold on;
for i=1:numel(Hcas)
    [gd,w_gd] = grpdelay(Hcas(i)); % group delay
    [phi,w_phi] = phasedelay(Hcas(i)); % phase delay
    maxPhi(i) = max(phi)/Fs*formants(featureNum);
    plot(Fs*w_gd/(2*pi),gd/Fs*formants(featureNum),'Linewidth',3,'Color','b');
    plot(Fs*w_phi/(2*pi),phi/Fs*formants(featureNum),'Linewidth',3,'Color','k');
end
xlabel('Frequency (Hz)');
ylabel('Delay (F2 cycles)');
set(gca,'XLim',[100 10e3],...
        'XTick',[100 250 500 1e3 2500 5e3 10e3],...
        'XScale','log');
legend('Group Delay','Phase Delay');

figure,
plot(mod(maxPhi,1),'o-');
xlabel('Number of all-pass filters');
ylabel(['max phase delay (F2 cycles)']);

figure,
plot([0 mod(maxPhi,1)],interp1(deltaCF,cycles',-0.5),'o');
xlabel('filter delay (cycles)');
ylabel('characteristic delay (cycles)');


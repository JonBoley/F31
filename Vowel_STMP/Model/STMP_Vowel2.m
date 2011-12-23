% This file analyzes the spatio-temporal coding of a vowel
% (using an auditory nerve model)
%
% Written by Jon Boley

% Add all subdirectories to the path
addpath(genpath('C:\Research\MATLAB\Vowel_STMP\Model'));

impaired = 1; % yes/no
amplified = 1; % turn on hearing aid?
featureNum = 1; % [1 2 3 ...] = [F1 F2 F3 ...]
Levels = 65;%[60 70 80 90]; % Should these actually be based on the RLF or BML?
SNRs = [Inf 6 0 -6]; % Is the noise actually getting loud enough in the speech band?

% Get stimulus
GenEh;

% Initialize model parameters
model_init; % default params

% Set up impairment for mild hearing loss (Bruce, ISAAR 2007)
impairment_init; % default params

% Set up hearing aid parameters
amplify_init;

% Initialize analysis variables
analysis_init;

for LevelIndex = 1:length(Levels)
    for SNRindex=1:length(SNRs)
        % Set the level
        OALevel_dBSPL = Levels(LevelIndex);
        fprintf('Level = %ddBSPL    ',Levels(LevelIndex))
        
        % Set the SNR
        [signal,delta_dB]=AddWhiteNoise(vowel,SNRs(SNRindex));
        fprintf('SNR = %ddB\n',SNRs(SNRindex))
        
        % Apply hearing aid gain
        if amplified
%             % Linear gain:
%             signal = ApplyFFTgain(signal,Fs,Audiogram_freq,NAL_IG);
            
            % Nonlinear gain:
            if OALevel_dBSPL==65
                % Based on DSL[i/o] pure tones in noise
                DSL_REAG=[-13 -2 10 25 30 31 25];
            else
                errordlg('No Gain specified for this level');
            end
            signal = ApplyFFTgain(signal,Fs,Audiogram_freq,DSL_REAG);
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
            
            % Run model
            get_spikes;
            
            % Calculate ALSR
            %  - Calculate Period Histogram
            %  - Calculate Synchronized Rate
            
            % Calculate cross-fiber correlations
            calc_xCF;
            
        end
        fprintf('\n');
        
        Rho_vDeltaCF(LevelIndex,SNRindex,:) = [Rho(2:midCF_index,LevelIndex,SNRindex); NaN; Rho(midCF_index+1:end,LevelIndex,SNRindex)];
        CD_vDeltaCF(LevelIndex,SNRindex,:) = [CD(2:midCF_index,LevelIndex,SNRindex); 0; CD(midCF_index+1:end,LevelIndex,SNRindex)];
        
        figure(1)
        plot(deltaCF,squeeze(Rho_vDeltaCF)); hold on;
        plot([0 0],[min(min(Rho_vDeltaCF)) max(max(Rho_vDeltaCF))],'k:');
        text(0,max(max(Rho_vDeltaCF)),FeaturesText{featureNum},...
            'HorizontalAlignment','Center'); hold off;
        title('\rho vs \DeltaCF');
        xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
        ylabel(['\rho (re ' FeaturesText{featureNum} ')']);
        pause(0.1);
    end
end

if length(SNRs)>1
    % Plot SNRs
    lines = {'bo-','rs-','b^-','gv-'};
    figure(2), subplot(2,1,1),
    for i=1:size(squeeze(Rho_vDeltaCF),1)
        plot(squeeze(deltaCF),squeeze(Rho_vDeltaCF(1,i,:)),lines{i}); hold on;
    end
    plot([0 0],[min(min(Rho_vDeltaCF)) max(max(Rho_vDeltaCF))],'k:');
    text(0,max(max(Rho_vDeltaCF)),FeaturesText{featureNum},...
        'HorizontalAlignment','Center'); hold off;
    if impaired
        if amplified
            title('Impaired + Amplified');
        else
            title('Impaired');
        end
    else
        if amplified
            title('Normal + Amplified');
        else
            title('Normal');
        end
    end
    xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
    ylabel(['\rho (re ' FeaturesText{featureNum} ')']);
    legend([num2str(SNRs') repmat('dB SNR',4,1)]);
    axis([-0.5 0.5 0.5 1.5]);
    
    figure(2), subplot(2,1,2),
    for i=1:size(squeeze(Rho_vDeltaCF),1)
        plot(squeeze(deltaCF),squeeze(CD_vDeltaCF(1,i,:)),lines{i}); hold on;
    end
    xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
    ylabel(['Characteristic Delay (\mus re ' FeaturesText{featureNum} ')']);
    legend([num2str(SNRs') repmat('dB SNR',4,1)]);
    axis([-0.5 0.5 -1000 1000]);
end

if length(Levels)>1
    % Plot Levels
    lines = {'bo-','rs-','b^-','gv-'};
    figure(2), subplot(2,1,1),
    for i=1:size(squeeze(Rho_vDeltaCF),1)
        plot(squeeze(deltaCF),squeeze(Rho_vDeltaCF(i,1,:)),lines{i}); hold on;
    end
    plot([0 0],[min(min(Rho_vDeltaCF)) max(max(Rho_vDeltaCF))],'k:');
    text(0,max(max(Rho_vDeltaCF)),FeaturesText{featureNum},...
        'HorizontalAlignment','Center'); hold off;
    if impaired
        title('Impaired');
    else
        title('Normal');
    end
    xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
    ylabel(['\rho (re ' FeaturesText{featureNum} ')']);
    legend([num2str(Levels') repmat('dB SPL',4,1)]);
    axis([-0.5 0.5 0.5 1.5]);
    
    figure(2), subplot(2,1,2),
    for i=1:size(squeeze(Rho_vDeltaCF),1)
        plot(squeeze(deltaCF),squeeze(CD_vDeltaCF(i,1,:)),lines{i}); hold on;
    end
    xlabel(['\DeltaCF (octaves re ' sprintf('%s)',FeaturesText{featureNum})]);
    ylabel(['Characteristic Delay (\mus re ' FeaturesText{featureNum} ')']);
    legend([num2str(Levels') repmat('dB SPL',4,1)]);
    axis([-0.5 0.5 -1000 1000]);
end


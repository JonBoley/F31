% This file just plots which gains were used

addpath(genpath('C:\Research\MATLAB\Vowel_STMP\Model'));

impaired = 1; % yes/no
amplified = 2; % [0 1 2] = [off linear nonlinear]
featureNum = 3; % [1 2 3 ...] = [F1 F2 F3 ...]
Levels = 65;%[60 70 80 90]; %65; % Should these actually be based on the RLF or BML?
SNRs = [Inf 6 0 -6]; % Is the noise actually getting loud enough in the speech band?

lines = {'go-','rs-','b^-','bv:'};

% Get stimulus
GenEh;

% Initialize model parameters
model_init; % default params

% Set up impairment for mild hearing loss (Bruce, ISAAR 2007)
impairment_init; % default params

% Set up hearing aid parameters
amplify_init;

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
            if amplified==1 
                % Linear gain:
                IG=NAL_IG;
            else
                % Nonlinear gain:
                % Based on DSL[i/o] in noise (monaural)
                if (OALevel_dBSPL+delta_dB)<=60
                    DSL_REAG=[-14 -3 7 14 16 17 18]; % speech REIG
                elseif (OALevel_dBSPL+delta_dB)<=65
                    DSL_REAG=[-14 -3 7 13 15 16 18]; % speech REIG
%                     DSL_REAG=[-13 -2 10 25 30 31 25]; % pure tone REAG
                elseif (OALevel_dBSPL+delta_dB)<=70
                    DSL_REAG=[-14 -3 7 12 14 15 17]; % speech REIG
                elseif (OALevel_dBSPL+delta_dB)<=80
                    DSL_REAG=[-14 -3 7 9 11 13 15]; % speech REIG
                elseif (OALevel_dBSPL+delta_dB)<=90 %(DSL only specifies up to 84dB)
                    DSL_REAG=[-14 -3 7 7 10 12 14]; % speech REIG
                else
                    errordlg('No Gain specified for this level');
                end
                IG=DSL_REAG;
            end
        else
            IG=[0 0 0 0 0 0 0];
        end
        
        
        for myformant=1:3
            Gain(SNRindex,myformant) = interp1(Audiogram_freq,IG,formants(myformant),'linear','extrap');
        end
    end
    for myformant=1:3
        figure(1), subplot(LevelIndex,3,myformant);
        plot([12 SNRs(2:end)],Gain(:,myformant),lines{impaired+amplified+1},'LineWidth',3);
        set(gca,'XTick',[fliplr(SNRs(2:end)) 12]);
        set(gca,'XTickLabel','-6|0|6|Quiet');
        set(gca,'Xdir','reverse');
        hold on;
    end
end


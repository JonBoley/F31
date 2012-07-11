function [vowel,noise,Fs,filename_vowel,filename_noise,newVowelAtten,...
    newNoiseAtten,dBreTONE,NEWforms_Hz]=...
    synth_BASELINE_ehLTASS_aid(TargetFreq_Hz,NEWF0_Hz,TargetFeature,...
    Fix2Harms,mode,aidParams)
% function [vowel,noise,Fs,filename_vowel,filename_noise,newVowelAtten,...
%     newNoiseAtten,dBreTONE,NEWforms_Hz]=...
%     synth_BASELINE_ehLTASS_aid(TargetFreq_Hz,NEWF0_Hz,TargetFeature,...
%     Fix2Harms,mode,aidParams)
% 
% aidParams.vowelAtten
% aidParams.noiseAttens
% aidParams.max_dBSPL
% aidParams.strategy
% aidParams.Fs_Hz
%
% J. Boley 07Mar2012 (from synth_BASELINE_eh.m)
% M. Heinz 30Jun2004 (from makeEH.m 04Nov2003)
% Make an /eh/ like in Conley and Keilson (1995)
%
% Creates a BASELINE vowel (based on ORIG eh) with F2 placed at the TargetFrequency (TF)
% and F0=75 Hz.  Expand BWs and formants based on shift, as would have occured with a Fs shift,
% but here F0 is not shifted, it is set to the desired value
%
% Also applies hearing aid gain, according to parameters supplied
%
% Original EH:
% F0=100Hz, Fs=10kHz;
%           F1     F2      F3      F4     F5
%  f (Hz)  500    1700    2500    3300   3750
% bw (Hz)   60      90     200     250    200
%
% Mode: 1: generate vowel and show plot
%       2: generate vowel only, don't show
%       3: only return filename and NEWforms_Hz
% strategy = 0 for 'none'
%         or 1 for 'linear'
%         or 2 for 'nonlinear_quiet'
%         or 3 for 'nonlinear_noise'
%

if ~exist('mode','var')
    mode=1;
end

ORIGforms_Hz=[500    1700    2500    3300   3750];
ORIGbws_Hz  =[60      90     200     250    200];
NUMforms=5;  % could use fewer formants
ORIGF0=100;
Fs=33000;

Features={'T0','F1','T1','F2','T2','F3','T3'};
if strcmp(TargetFeature,'F2')
    ORIGfreq=ORIGforms_Hz(2);
else
    error('Only set up for F2 as Target Feature')
end
ShiftFact=TargetFreq_Hz/ORIGfreq;

% Translate both the formant frequencies and bandwidths, as would occur with a shift in Fs
NEWforms_Hz=ORIGforms_Hz*ShiftFact;
NEWbws_Hz=ORIGbws_Hz*ShiftFact;

%%%%%%%%%%%% Fix formants to harmonics if option is set
FHtext='';
AdjustedHarms=[];
if Fix2Harms
    FHtext='_FH';
    for ind=1:NUMforms
        if rem(NEWforms_Hz(ind),NEWF0_Hz)
            AdjustedHarms=[AdjustedHarms ind];
            NEWforms_Hz(ind)=round(NEWforms_Hz(ind)/NEWF0_Hz)*NEWF0_Hz;
        end
    end
end
filename_vowel=sprintf('baseEH_F2at%05.f_F0at%03.f%s.wav',TargetFreq_Hz,NEWF0_Hz,FHtext);

% Mode 3 only returns info; Modes 1 and 2 return the synthesized vowel
if mode<3
    disp('***New vowel synthesized')
    if Fix2Harms
        disp(sprintf('     - Formants: %s adjusted to be at harmonics',mat2str(AdjustedHarms)))
    end
    %% returned vowel will be 1 cycle, but any analysis and/or plotting will be done with a longer vowel
    FULLdur=0.6;   % This will be made to an integer number of periods (<=FULLdur)
    dur=1/NEWF0_Hz;
    Nreps=ceil(FULLdur/dur);  % Create full-duration waveform

    [time, vowel] = dovowl(NEWforms_Hz(1:NUMforms), NEWbws_Hz(1:NUMforms), NEWF0_Hz, dur, Fs);

    vowel=vowel'-mean(vowel);  % Remove DC
    vowel=vowel/max(abs(vowel))*.999;  % Need amplitude to fit in wavwrite

    noise=GenLTASS(dur,Fs); % Generate noise
    noise=noise/max(abs(noise))*.999;  % Need amplitude to fit in wavwrite

    % Mode 2 uses getVowelParams to calculate params; mode 2 also will view the vowel
    if mode==1
        viewVowel=1;
        %       FULLvowel=repmat(vowel,1,Nreps);
        %       sound(FULLvowel,Fs)
        %       disp(sprintf('\nPlaying New Synthesized Vowel'))
        %       pause(1)
        %    [vowelORIG,FsORIG]=wavread('vow17_MH10k.wav');
        %        sound(vowelORIG,FsORIG)
        %    disp('Playing Resampled ORIGINAL Vowel (EH) from CN exps')
    else
        viewVowel=0;
    end
    [feat_freqs_Hz,feat_levs_dB,dBreTONE]=getVowelParams(vowel,Fs,FULLdur,NEWforms_Hz,viewVowel);
else
    vowel=[]; Fs=[]; dBreTONE=[]; noise=[];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY HEARING AID GAIN
if aidParams.strategy
    switch aidParams.strategy
        case 1, strategyText='_linearAid';
        case 2, strategyText='_nonlinearAid';
        case 3, strategyText='_nonlinearNoiseAid';
        otherwise, strategyText='';
    end
    filename_vowel=sprintf('baseEH_F2at%05.f_F0at%03.f%s%s_vowelat%03.fdB_noiseat%03.fdB.wav',...
        TargetFreq_Hz,NEWF0_Hz,FHtext,strategyText,aidParams.vowelAtten,aidParams.noiseAtten);
    filename_noise=sprintf('baseLTASS%s_vowelat%03.fdB_noiseat%03.fdB.wav',...
        strategyText,aidParams.vowelAtten,aidParams.noiseAtten);

    if ~isempty(vowel) && ~isempty(noise)
        % 0) Combine vowel & noise
        input = vowel + (10^(aidParams.vowelAtten-aidParams.noiseAtten)*noise);

        % 1) Calculate hearing aid amplification
        audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
        freqs_Hz = [500,1000,2000,4000,6000];
        gainStruct = calcgain(input,aidParams.Fs_Hz,...
            aidParams.max_dBSPL,aidParams.vowelAtten,...
            audiogram,freqs_Hz,aidParams.strategy);

        % 2) Apply amplification to signal & noise separately
        MaxLevel0 = max(abs(vowel));
        MaxLevel1 = max(abs(noise));

        % frequency shift of CF (@F2) relative to base F2
        freqShift = TargetFreq_Hz / 1700;
        fc = freqShift * 500*2^(0.5*log2(1700/500)); % cross-over frequency (between F1 & F2)

        switch aidParams.strategy
            case 1 % linear
                vowel = FIRgain(vowel,gainStruct.linear.gain_dB,...
                    gainStruct.linear.f_Hz*freqShift,aidParams.Fs_Hz);
                noise = FIRgain(noise,gainStruct.linear.gain_dB,...
                    gainStruct.linear.f_Hz*freqShift,aidParams.Fs_Hz);
            case 2 % nonlinear (quiet)
                vowel = TwoBandGain(vowel,fc,aidParams.Fs_Hz,...
                    gainStruct.nonlinear.quiet.gains);
                noise = TwoBandGain(noise,fc,aidParams.Fs_Hz,...
                    gainStruct.nonlinear.quiet.gains);
            case 3 % nonlinear (noise)
                vowel = TwoBandGain(vowel,fc,aidParams.Fs_Hz,...
                    gainStruct.nonlinear.noise.gains);
                noise = TwoBandGain(noise,fc,aidParams.Fs_Hz,...
                    gainStruct.nonlinear.noise.gains);
            otherwise
                % no amplification
        end

        % 3) Calculate how much the peak changed (e.g., increased)
        MaxLevelRatio0 = max(abs(vowel))/MaxLevel0;
        MaxLevelRatio1 = max(abs(noise))/MaxLevel1;

        % 4) Change scale back to original dynamic range (to prevent clipping)
        vowel = vowel / MaxLevelRatio0;
        noise = noise / MaxLevelRatio1;

        % 5) Adjust attenuation to get the intended level (e.g., decrease)
        newVowelAtten = aidParams.vowelAtten - 20*log10(MaxLevelRatio0);
        newNoiseAtten = aidParams.noiseAtten - 20*log10(MaxLevelRatio1);
    end %if ~isempty(vowel)
end %if aidParams.strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





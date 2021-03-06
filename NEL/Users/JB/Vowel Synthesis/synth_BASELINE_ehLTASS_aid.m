function [mixed,Fs,filename_mixed,newMixedAtten,dBreTONE,NEWforms_Hz]=...
    synth_BASELINE_ehLTASS_aid(TargetFreq_Hz,NEWF0_Hz,TargetFeature,...
    Fix2Harms,mode,stimDur,aidParams)
% function [vowel,noise,Fs,filename_vowel,filename_noise,newVowelAtten,...
%     newNoiseAtten,dBreTONE,NEWforms_Hz]=...
%     synth_BASELINE_ehLTASS_aid(TargetFreq_Hz,NEWF0_Hz,TargetFeature,...
%     Fix2Harms,mode,stimDur,aidParams)
% 
% aidParams.vowelAtten
% aidParams.noiseAtten
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

global signals_dir
EHsignals_dir=strcat(signals_dir,'JB\EHvowels');

if ~exist('mode','var')
    mode=1;
end

ORIGforms_Hz=[500    1700    2500    3300   3750];
ORIGbws_Hz  =[60      90     200     250    200];
NUMforms=5;  % could use fewer formants
ORIGF0=100;
Fs=33000;

Features={'T0','F1','T1','F2','T2','F3','T3'};
if strcmp(TargetFeature,'F1')
    ORIGfreq=ORIGforms_Hz(1);
elseif strcmp(TargetFeature,'F2')
    ORIGfreq=ORIGforms_Hz(2);
else
    error('Only set up for F1 or F2 as Target Feature')
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
filename_vowel=sprintf('baseEH_%sat%05.f_F0at%03.f%s.wav',TargetFeature,TargetFreq_Hz,NEWF0_Hz,FHtext);

% Mode 3 only returns info; Modes 1 and 2 return the synthesized vowel
if mode<3
%     disp('***New vowel synthesized')
    if Fix2Harms
        disp(sprintf('     - Formants: %s adjusted to be at harmonics',mat2str(AdjustedHarms)))
    end
%     FULLdur=0.6;   % This will be made to an integer number of periods (<=FULLdur)
%     dur=1/NEWF0_Hz;
%     Nreps=ceil(FULLdur/dur);  % Create full-duration waveform

    [time, vowel] = dovowl(NEWforms_Hz(1:NUMforms), NEWbws_Hz(1:NUMforms), NEWF0_Hz, stimDur, Fs);

    vowel=vowel'-mean(vowel);  % Remove DC
    vowel=vowel/max(abs(vowel))*.999;  % Need amplitude to fit in wavwrite

    noise=GenLTASS(stimDur,Fs); % Generate noise
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
    [feat_freqs_Hz,feat_levs_dB,dBreTONE]=getVowelParams(vowel,Fs,stimDur,NEWforms_Hz,viewVowel);
else
    vowel=[]; Fs=[]; dBreTONE=[]; noise=[];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY HEARING AID GAIN
CalcSeparately = 0; % Calculate signal & noise separately?

mixed = vowel + (10^((aidParams.vowelAtten-aidParams.noiseAtten)/20)*noise);
MaxLevelRatio2 = max(abs(mixed))/max(abs(vowel)); %this one relative to vowel alone
mixed = mixed / MaxLevelRatio2;
newMixedAtten = aidParams.vowelAtten - 20*log10(MaxLevelRatio2);

switch aidParams.strategy
case 1, strategyText='_linearAid';
case 2, strategyText='_nonlinearAid';
case 3, strategyText='_nonlinearNoiseAid';
otherwise, strategyText='';
end

filename_vowel=sprintf('baseEH_%sat%05.f_F0at%03.f%s%s_%03.fdBEH_%03.fdBLTASS.wav',...
    TargetFeature,TargetFreq_Hz,NEWF0_Hz,FHtext,strategyText,aidParams.vowelAtten,aidParams.noiseAtten);
filename_noise=sprintf('baseLTASS%s_%03.fdBEH_%03.fdBLTASS.wav',...
    strategyText,aidParams.vowelAtten,aidParams.noiseAtten);
filename_mixed=sprintf('baseEHLTASS_%sat%05.f_F0at%03.f%s%s_%03.fdBEH_%03.fdBLTASS.wav',...
    TargetFeature,TargetFreq_Hz,NEWF0_Hz,FHtext,strategyText,aidParams.vowelAtten,aidParams.noiseAtten);


if length(dir(fullfile(EHsignals_dir,[filename_mixed(1:end-4) '*.wav'])))
    existingFile = dir(fullfile(EHsignals_dir,[filename_mixed(1:end-4) '*.wav']));
    if length(existingFile)>1
        error(['Multiple file matches found matching ' filename_mixed(1:end-4) '*.wav']);
    end
    newMixedAtten=str2num(existingFile.name(strfind(existingFile.name,'dBWAV')-3:strfind(existingFile.name,'dBWAV')-1));
    fprintf('File found: %s\n... Using Atten=%1.fdB\n',existingFile.name,newMixedAtten);
    filename_mixed = [filename_mixed(1:end-4) sprintf('_%03.fdBWAV',newMixedAtten) '.wav'];
    return; % if this file already exists, just return
elseif aidParams.strategy & ~isempty(vowel) & ~isempty(noise)
    % 0) Combine vowel & noise
    mixed = vowel + (10^((aidParams.vowelAtten-aidParams.noiseAtten)/20)*noise);
    
    % 1) Calculate hearing aid amplification
    audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
    freqs_Hz = [500,1000,2000,4000,6000];
    gainStruct = calcgain(mixed,aidParams.Fs_Hz,...
        aidParams.max_dBSPL,aidParams.vowelAtten,...
        audiogram,freqs_Hz,aidParams.strategy);
    
    % 2) Apply amplification to signal & noise separately
    MaxLevel0 = max(abs(vowel));
    MaxLevel1 = max(abs(noise));
    MaxLevel2 = max(abs(mixed));
    
    % frequency shift of CF (@F2) relative to base F2
    if strcmp(TargetFeature,'F1')
        freqShift = TargetFreq_Hz / ORIGforms_Hz(1);
    elseif strcmp(TargetFeature,'F2')
        freqShift = TargetFreq_Hz / ORIGforms_Hz(2);
    else
        error('Only F1 and F2 formants are currently supported!');
    end
    fc = freqShift * 500*2^(0.5*log2(1700/500)); % cross-over frequency (between F1 & F2)
    
    switch aidParams.strategy
    case 1 % linear
        if CalcSeparately
            vowel = FIRgain(vowel,gainStruct.linear.gain_dB,...
                gainStruct.linear.f_Hz*freqShift,aidParams.Fs_Hz);
            noise = FIRgain(noise,gainStruct.linear.gain_dB,...
                gainStruct.linear.f_Hz*freqShift,aidParams.Fs_Hz);
        else
            mixed = FIRgain(mixed,gainStruct.linear.gain_dB,...
                gainStruct.linear.f_Hz*freqShift,aidParams.Fs_Hz);
        end
    case 2 % nonlinear (quiet)
        if CalcSeparately
            vowel = TwoBandGain(vowel,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.quiet.gains);
            noise = TwoBandGain(noise,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.quiet.gains);
        else
            mixed = TwoBandGain(mixed,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.quiet.gains);
        end
    case 3 % nonlinear (noise)
        if CalcSeparately
            vowel = TwoBandGain(vowel,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.noise.gains);
            noise = TwoBandGain(noise,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.noise.gains);
        else
            mixed = TwoBandGain(mixed,fc,aidParams.Fs_Hz,...
                gainStruct.nonlinear.noise.gains);
        end
    otherwise
        % no amplification
    end
    
    % 3) Calculate how much the peak changed (e.g., increased)
    MaxLevelRatio0 = max(abs(vowel))/MaxLevel0;
    MaxLevelRatio1 = max(abs(noise))/MaxLevel1;
    MaxLevelRatio2 = max(abs(mixed))/MaxLevel0; %this one relative to vowel alone
    
    % 4) Change scale back to original dynamic range (to prevent clipping)
    vowel = vowel / MaxLevelRatio0;
    noise = noise / MaxLevelRatio1;
    mixed = mixed / MaxLevelRatio2;
    
    % 5) Adjust attenuation to get the intended level (e.g., decrease)
    newVowelAtten = aidParams.vowelAtten - 20*log10(MaxLevelRatio0);
    newNoiseAtten = aidParams.noiseAtten - 20*log10(MaxLevelRatio1);
    newMixedAtten = aidParams.vowelAtten - 20*log10(MaxLevelRatio2);
    
    filename_vowel = [filename_vowel(1:end-4) sprintf('_%03.fdBWAV',newVowelAtten) '.wav'];
    filename_noise = [filename_noise(1:end-4) sprintf('_%03.fdBWAV',newNoiseAtten) '.wav'];
    
end %if ~isempty(vowel) (& filename doesn't exist)

filename_mixed = [filename_mixed(1:end-4) sprintf('_%03.fdBWAV',newMixedAtten) '.wav'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





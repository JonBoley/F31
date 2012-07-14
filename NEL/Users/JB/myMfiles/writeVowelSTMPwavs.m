function newMixedAtten = writeVowelSTMPwavs(CalibPicNum,BF_kHz,...
    Vowel_dBSPL,noiseAttenMid,strategy)
% function newMixedAtten = writeVowelSTMPwavs(CalibPicNum,BF_kHz,...
%     Vowel_dBSPL,noiseAttenMid,strategy)
%
% strategy = 0 for 'none'
%         or 1 for 'linear'
%         or 2 for 'nonlinear_quiet'
%         or 3 for 'nonlinear_noise'

global signals_dir
EHsignals_dir=strcat(signals_dir,'JB\EHvowels');

if isempty(dir(sprintf('p%04d_calib.m',CalibPicNum)))
    error('Cannot find calib pic!');
else
    x=loadpic(CalibPicNum);
    CalibData=x.CalibData(:,1:2);
    CalibData(:,2)=trifilt(CalibData(:,2)',5)';
    aidParams.max_dBSPL = CalibInterp(BF_kHz,CalibData);
end
aidParams.strategy = strategy;
aidParams.Fs_Hz = 33000;

F0=100;
Fix2Harms=0;

%%% generate temp stimuli to get dBreTONE
ORIGforms_Hz=[500    1700    2500    3300   3750];
ORIGbws_Hz  =[60      90     200     250    200];
[time, vowel] = dovowl(ORIGforms_Hz, ORIGbws_Hz, F0, 1, aidParams.Fs_Hz);
vowel=vowel'-mean(vowel);  % Remove DC
vowel=vowel/max(abs(vowel))*.999;  % Need amplitude to fit in wavwrite

noise=GenLTASS(1,aidParams.Fs_Hz); % Generate noise
noise=noise-mean(noise);           % Remove DC
noise=noise/max(abs(noise))*.999;  % Need amplitude to fit in wavwrite

dBreTONE=20*log10(sqrt(mean(vowel.^2))/.707);
dBreTONE_noise=20*log10(sqrt(mean(noise.^2))/.707);
clear vowel noise;
%%%

aidParams.vowelAtten = max_dBSPL-Vowel_dBSPL+dBreTONE;
noiseAttens(1) = 120; % quiet
noiseAttens(2) = noiseAttenMid; % Equal SL
noiseAttens(3) = max_dBSPL-Vowel_dBSPL+dBreTONE_noise; % Equal SPL

newNoiseAtten = NaN*ones(size(noiseAttens));
for AttenIndex = 1:length(noiseAttens)
    aidParams.noiseAtten = noiseAttens(AttenIndex);
    
    for TargetFeature = {'F1','F2'}
        [mixed,Fs,filename_mixed,newMixedAtten,dBreTONE,NEWforms_Hz]=...
            synth_BASELINE_ehLTASS_aid(BF_kHz*1e3,F0,TargetFeature,...
            Fix2Harms,2,aidParams);
        
        disp(sprintf('Writing %s ...',filename_mixed));
        if ~exist(fullfile(EHsignals_dir,filename_mixed),'file')
            wavwrite(mixed,Fs,fullfile(EHsignals_dir,filename_mixed));
        end
        
        %%% NEED TO ADD FILENAMES TO LIST
        % <signals_dir>\Lists\JB\VowelLTASS\VowelLTASS.m
        
    end
end



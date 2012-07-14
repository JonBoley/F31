function [newMixedAtten,filename_mixed] = writeVowelSTMPwavs(CalibPicNum,BF_kHz,...
    Vowel_dBSPL,noiseAttenMid,stimDur,strategy)
% function newMixedAtten = writeVowelSTMPwavs(CalibPicNum,BF_kHz,...
%     Vowel_dBSPL,noiseAttenMid,stimDur,strategy)
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
aidParams.Fs_Hz = 33000;

F0=100;
Fix2Harms=0;

%%% generate temp stimuli to get dBreTONE
ORIGforms_Hz=[500    1700    2500    3300   3750];
ORIGbws_Hz  =[60      90     200     250    200];
[time, vowel] = dovowl(ORIGforms_Hz, ORIGbws_Hz, F0, 1, aidParams.Fs_Hz);
vowel=vowel'-mean(vowel);  % Remove DC
vowel=vowel/max(abs(vowel))*.999;  % Need amplitude to fit in wavwrite
dBreTONE=20*log10(sqrt(mean(vowel.^2))/.707);

if isempty(which('GenLTASS')), addpath(genpath('C:\NEL\Users\JB\LTASS')); end
% noise=GenLTASS(1,aidParams.Fs_Hz); % Generate noise
% noise=noise-mean(noise);           % Remove DC
% noise=noise/max(abs(noise))*.999;  % Need amplitude to fit in wavwrite
dBreTONE_noise = -9.6039; %20*log10(sqrt(mean(noise.^2))/.707);

clear vowel noise;
%%%

aidParams.vowelAtten = aidParams.max_dBSPL-Vowel_dBSPL+dBreTONE;
noiseAttens(1) = 120; % quiet
noiseAttens(2) = noiseAttenMid; % Equal SL
noiseAttens(3) = aidParams.max_dBSPL-Vowel_dBSPL+dBreTONE_noise; % Equal SPL
strNoiseConds = {'_quiet','_equalSL','_equalSPL'};

newMixedAtten = NaN*ones(2,length(noiseAttens),length(strategy));
for strategyNum = 1:length(strategy)
    aidParams.strategy = strategy(strategyNum);
    
    for AttenIndex = 1:length(noiseAttens)
        aidParams.noiseAtten = noiseAttens(AttenIndex);
        
        featIndex=1;
        for TargetFeature = {'F1','F2'}
            [mixed,Fs,filename_mixed{featIndex,AttenIndex,strategyNum},...
                    newMixedAtten(featIndex,AttenIndex,strategyNum),dBreTONE,NEWforms_Hz]=...
                synth_BASELINE_ehLTASS_aid(BF_kHz*1e3,F0,TargetFeature{1},...
                Fix2Harms,2,stimDur,aidParams);
            
            % add {'_quiet','_equalSL','_equalSPL'} to end of filename_mixed
            filename_mixed{featIndex,AttenIndex,strategyNum} = ...
                [filename_mixed{featIndex,AttenIndex,strategyNum}(1:end-4),...
                    strNoiseConds{AttenIndex}, '.wav'];
            
            if ~exist(fullfile(EHsignals_dir,filename_mixed{featIndex,AttenIndex,strategyNum}),'file')
                disp(sprintf('Writing %s ...',filename_mixed{featIndex,AttenIndex,strategyNum}));
                if max(abs(mixed))>0.9995
                    fprintf('... FILE MAY BE CLIPPED! (max=%1.5f)',max(abs(mixed))); 
                end
                wavwrite(mixed,Fs,fullfile(EHsignals_dir,filename_mixed{featIndex,AttenIndex,strategyNum}));
            end
            
            featIndex=featIndex+1;
        end
    end
end



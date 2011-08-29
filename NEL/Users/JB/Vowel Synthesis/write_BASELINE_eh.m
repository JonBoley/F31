% File: write_BASELINE_eh.m
%

TargetFreq_Hz=2000;
NEWF0_Hz=75;
TargetFeature='F2';
Fix2Harms=0;
EHsignals_dir='C:\Signals\MH\EHvowels';

[stim,Fs,filename,dBreTONE]=synth_BASELINE_eh(TargetFreq_Hz,NEWF0_Hz,TargetFeature,Fix2Harms);

wavwrite(stim,Fs,fullfile(EHsignals_dir,filename))

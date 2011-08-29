function [vowel,dBreTONE,filename]= sfs_synth(forms_Hz,F0,dur,Fs)
% FIle: synth_ORIG_eh.m
%
% M. Heinz 30Jun2004 (from makeEH.m 04Nov2003)
% Make an /eh/ like in Conley and Keilson (1995)
%
% Original EH:
% F0=100Hz, Fs=10kHz;
%           F1     F2      F3      F4     F5
%  f (Hz)  500    1700    2500    3300   3750
% bw (Hz)   60      90     200     250    200

%clear

% forms_Hz=[500    1700    2500    3300   3750];
% bws_Hz  =[60      90     200     250    200];
%forms_Hz=[1500];
%forms_Hz=[1700];
%bws_Hz  =[60];
%F0=100;
%Fs=48828.125;
%Fs=48828
%dur=0.5;
% dur=1/F0;
% NUMforms=5;
NUMforms=1;
bws_Hz  =[60];

[time, vowel] = dovowl(forms_Hz(1:NUMforms), bws_Hz(1:NUMforms), F0, dur, Fs);

vowel=vowel-mean(vowel);
vowel=vowel/max(vowel);

sound(vowel,Fs)
%disp(sprintf('\nPlaying Synthesized Vowel'))

% viewVowel(vowel,Fs,forms_Hz)
[dBreTONE,RMS]=viewvowel_mod(vowel,Fs,forms_Hz)

vowel=vowel';
filename=sprintf('SFS_CFat%0.1fkHZ.wav',forms_Hz(1)/1000);
%wavwrite(vowel,Fs,'EH.wav')



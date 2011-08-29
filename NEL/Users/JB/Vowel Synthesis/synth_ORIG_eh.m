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

clear

forms_Hz=[500    1700    2500    3300   3750];
bws_Hz  =[60      90     200     250    200];
F0=100;
Fs=33000;
dur=0.6;
% dur=1/F0;
NUMforms=5;

[time, vowel] = dovowl(forms_Hz(1:NUMforms), bws_Hz(1:NUMforms), F0, dur, Fs);

vowel=vowel-mean(vowel);
vowel=vowel/max(vowel);

sound(vowel,Fs)
disp(sprintf('\nPlaying Synthesized Vowel'))
pause(1)
[stimORIG,FsORIG]=wavread('vow17_MH10k.wav');
sound(stimORIG,FsORIG)
disp('Playing Resampled Vowel from CN exps')

viewVowel(vowel,Fs,forms_Hz)




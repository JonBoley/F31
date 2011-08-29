% File: write_BASELINE_noise.m
%
% 11-Apr-2005 M. Heinz
% generating white noise for use with EH-vowel studies
%
% Needs to be the same SAMPLING RATE (Fs) as BASELINE EH (33000)
% Want to have at least 400 ms of independent noise (non-repeated) when played at highest 
%    possible Fs for NI board with 2 signals (167 kHz), i.e., 167000*.4 = 66800 (==> USE 70000)
% 1) generate 70000 samples of random normally distributed noise,
% 2) Normalize to max value =1
% 3) Save as wavfile with Fs=33000 Hz
% 4) Use in templates, and compute dBreTONE as 20*log10(RMSnoise/.707)

Nsamps=70000;
Fs=33000;

noise=randn(1,Nsamps);
noise=noise/max(abs(noise))*0.9999;  % avoid warning you sometimes get with a real max of 1

[min(noise) max(noise)]

wavwrite(noise,Fs,fullfile(signals_dir,'MH','EHvowels','baseNOISE.wav'))

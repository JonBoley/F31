function out=NoiseGen(len, depth, freq, Fs)
% Generate AM White Noise
% 
% out = NoiseGen(len, depth, freq, Fs)
% len = length of noise sequence (in sec)
% depth = amplitude modulation depth (in percent)
% freq = modulation frequency
% Fs = sample rate (in Hz)

% generate noise
Y = randn(1,len*Fs);  % Gaussian White Noise (zero mean, unit variance)

% generate sinusoidal modulation tone
t = 1/Fs:1/Fs:len;   % "len" seconds @ "Fs" sample rate
tone = sin(2*pi*freq*t);
tone = 0.5*depth*tone/100;  % adjust modulation depth
tone = tone + 1-max(tone);  % bring max up to one

% apply modulation
out = tone.*Y;

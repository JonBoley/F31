threshold = -20;   % dBFS
ratio     = 4;    % dB/dB
attack    = 0.005; % sec
release   = 0.005; % sec
gain      = 20;     % dB

ModFreq = 4;
ModDepth = 100;
SNR = -6;  % dB


[FileName,PathName,FilterIndex] = uigetfile('*.wav','Select Input Audio File','A_Boy_Fell.wav');
% [FileName2,PathName2,FilterIndex2] = uiputfile('*.wav','Select Output Audio File','out.wav');

[input,Fs,bitdepth] = wavread([PathName FileName]);

% Generate Noise
noise = NoiseGen(length(input)/Fs, ModDepth, ModFreq, Fs);
RMS1 = norm(input)/sqrt(length(input));
RMS2 = norm(noise)/sqrt(length(noise));
noise = noise * RMS1/RMS2 * 10^(-SNR/20);  % Normalize noise for long-term SNR

% Mix signal and noise
input = input' + noise;

output  = compressor(input,threshold,ratio,gain,release,attack,Fs,bitdepth);

% wavwrite(output,Fs,[PathName2 FileName2]);

subplot(2,1,1); plot(input);
title('Pre');
subplot(2,1,2); plot(output);
title('Post');


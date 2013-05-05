% Test Instantaneous Frequency

% Fs=48000;
% t = 0:1/Fs:2;
% y = chirp(t,80,t(end),300);%,'logarithmic');

[y,Fs]=wavread('heed_100.wav');

% filter
d = fdesign.bandpass(40,80,300,600,60,1,60,Fs);
Hd = design(d,'ellip');
y = filter(Hd,y);

xa = hilbert(y);
phi = unwrap(angle(xa));
freq = diff(phi)/(2*pi)*Fs;
freq = min(max(freq,80),300);
plot(freq);
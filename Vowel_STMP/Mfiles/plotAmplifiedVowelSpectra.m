%%
Fs=22e3;
MaxSPL = 105;
atten = 50;
audiogram = [16   18   20    9    9]; % based on avg animal data (500OBN exposure)
freqs_Hz = [500,1000,2000,4000,6000];
strategy = 2; %'nonlinear_quiet';
plotYes=1;
nLPC = 50;

forms_Hz=[500    1700    2500    3300   3750];
bws_Hz  =[60      90     200     250    200];
F0=100;
dur=0.6;
% dur=1/F0;
NUMforms=5;

addpath('C:\Research\MATLAB\NEL\Users\JB\Vowel Synthesis');
[time, vowel] = dovowl(forms_Hz(1:NUMforms), bws_Hz(1:NUMforms), F0, dur, Fs);

vowel=vowel-mean(vowel);
vowel=vowel/max(vowel);

[a,g] = lpc(vowel,nLPC);
[h,w]=freqz(1,a);

figure, hold on;
semilogx(Fs*w/(2*pi),20*log10(abs(h)),'k','LineWidth',3);
set(gca,'XScale','log'); set(gca,'XTick',1e3*[0.25 0.5 1 2 4]);
xlim([250 4000]); ylim([0 60]);

addpath('C:\Research\MATLAB\NEL\Users\JB\Amplification');
output_lin = ApplyGain(vowel,Fs,MaxSPL,atten,audiogram,freqs_Hz,1);
[a,g] = lpc(output_lin,nLPC);
[h,w]=freqz(1,a);
semilogx(Fs*w/(2*pi),20*log10(abs(h)),'m','LineWidth',3);

output_nonlin = ApplyGain(vowel,Fs,MaxSPL,atten,audiogram,freqs_Hz,2);
[a,g] = lpc(output_nonlin,nLPC);
[h,w]=freqz(1,a);
semilogx(Fs*w/(2*pi),15+20*log10(abs(h)),'g','LineWidth',3);

legend('No Gain','Linear','NonLinear');

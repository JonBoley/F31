% File: resampleEH_CNexps.m
% From: readEH.m (CNexps)
%
% M. Heinz 04Nov03
% Take vow17.wav file from Eric's directory and down-sample to 10 kHz
% Original Fs: 51281, dur=1 cycle
% 1) Make full-duration waveform (600 ms)
% 2) Down-sample
% 3) Find formant/trough features
%
% USED to make EH for CN Exps with GE/MH

clear

set(0,'DefaultTextInterpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur=0.6;
[stim,Fs]=wavread('EDY - EH Sounds/vow17.wav');
stim=stim';
durN=length(stim)/Fs;
Nreps=ceil(dur/durN);  % Create full-duration waveform
stimFull=repmat(stim,1,Nreps);
time=(1:length(stimFull))/Fs;

sound(stimFull,Fs)

pause(dur)
%wavwrite(stimFull,Fs,'vow17_MH.wav')

RMS=sqrt(mean(stimFull.^2));
RMS_TONE=1/sqrt(2);
dBreTONE=20*log10(RMS/RMS_TONE);

figure(1);
subplot(211)
plot(time,stimFull)
ylabel('Amplitude')
xlabel('Time (ms)')
title(sprintf('Vowel time waveform (vow17_MH.wav: %.2f dB re: TONE)',dBreTONE))

%%% FFT
vowel_fft=fft(stimFull);
Nfft=length(vowel_fft);
MagSpect=20*log10(abs(vowel_fft));
MagSpect=MagSpect-max(MagSpect);

freq=(0:Nfft-1)*Fs/Nfft/1000;
ind1=find(freq<=5);   % Find all power in 1st 50 harmonics
ind2=find(freq<=(Fs/2/1000)); % All positive freqs

feat_inds=[184 306 733 1038 1343 1526 1831];
feat_freqs=freq(feat_inds)*1000;
disp(sprintf('Feature Frequencies (T0,F1,T1,F2,T2,F3,T3):'))
for find=1:length(feat_inds)
   disp(sprintf('%.f Hz',feat_freqs(find)))
end
   
subplot(212)
plot(freq(ind2),MagSpect(ind2),'g')
hold on
plot(freq(ind1),MagSpect(ind1),'b')
plot(freq(feat_inds),MagSpect(feat_inds),'rx')
ylabel('Magnitude (dB)')
xlabel('Frequency (kHz)')
title('Vowel Spectrum')
xlim([0 5])
ylim([-50 0])

text(.75,.9,sprintf('Fs = %.2f Hz',Fs),'Units','norm')
text(.65,.8,sprintf('Feature Frequencies (T0,F1,T1,F2,T2,F3,T3):'),'Units','norm')
for find=1:length(feat_inds)
   text(.75,.8-.1*find,sprintf('%.f Hz',feat_freqs(find)),'Units','norm')
end

hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New sampling frequency
% resample command applies a anti-aliasing filter during the conversion!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs_newNom=10e3;
%Resample down to 10-kHz (newFs=oldFs*P/Q)
% Old=51281, new =10000 ==> P=10000; Q=51281;  %% Errors when I run using high Ps&Qs like this
P=round(Fs_newNom/10); Q=round(Fs/10);  %Integers used to up sample: Fsold*P/Q=Fsnew: (51281)*1000/5128=10.2k
Nfir=100;  % proportional to FIR filter length used for resampling: higher Nfir, better accuracy & longer comp time

if round(P)~=P,   error(sprintf('Resampling error: P (=%.2f) is not an integer!!!',P)), end
if round(Q)~=Q,   error(sprintf('Resampling error: Q (=%.2f) is not an integer!!!',Q)), end
if(P/Q*Fs==Fs_newNom) 
   disp('Integer Fs conversion exact')
   Fs_new=Fs_newNom;
else 
   disp(sprintf('Actual Fs_new = %.2f Hz',P/Q*Fs)) 
   Fs_new=P/Q*Fs;
end
disp(sprintf('Resampling with Nfir=%d .........',Nfir))

stimNew=resample(stimFull,P,Q,Nfir);

time=(1:length(stimNew))/Fs_new;

sound(stimNew,Fs_new)
%wavwrite(stimNew,Fs_new,'vow17_MH10k.wav')

RMS=sqrt(mean(stimNew.^2));
RMS_TONE=1/sqrt(2);
dBreTONE=20*log10(RMS/RMS_TONE);

figure(2);
subplot(211)
plot(time,stimNew)
ylabel('Amplitude')
xlabel('Time (ms)')
title(sprintf('Vowel time waveform (vow17_MH10k.wav: %.2f dB re: TONE)',dBreTONE))

%%% FFT
vowel_fft=fft(stimNew);
Nfft=length(vowel_fft);
MagSpect=20*log10(abs(vowel_fft));
MagSpect=MagSpect-max(MagSpect);

freq=(0:Nfft-1)*Fs_new/Nfft/1000;
% ind1=find(freq<=5);   % Find all power in 1st 50 harmonics
% ind2=find(freq<=(Fs_new/2/1000)); % All positive freqs

ind1=1:length(freq);
ind2=ind1;

% feat_inds=[184 306 733 1038 1343 1526 1831];
% feat_freqs=freq(feat_inds)*1000;
% disp(sprintf('Feature Frequencies (T0,F1,T1,F2,T2,F3,T3):'))
% for find=1:length(feat_inds)
%    disp(sprintf('%.f Hz',feat_freqs(find)))
% end
   
subplot(212)
plot(freq(ind2),MagSpect(ind2),'g')
hold on
plot(freq(ind1),MagSpect(ind1),'b')
plot(freq(feat_inds),MagSpect(feat_inds),'rx')
ylabel('Magnitude (dB)')
xlabel('Frequency (kHz)')
title('Vowel Spectrum')
xlim([0 5])
ylim([-50 0])

text(.75,.9,sprintf('Fs_new = %.2f Hz',Fs_new),'Units','norm')
text(.65,.8,sprintf('Feature Frequencies (T0,F1,T1,F2,T2,F3,T3):'),'Units','norm')
for find=1:length(feat_inds)
   text(.75,.8-.1*find,sprintf('%.f Hz',feat_freqs(find)),'Units','norm')
end

hold off
 
% print -depsc vow17_MH10k.ps

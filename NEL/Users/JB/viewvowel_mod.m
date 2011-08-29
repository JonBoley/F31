function [RMS ,dBreTONE]=viewvowel_mod (stim,Fs,forms_Hz)
% File: viewVowel.m
% From: resampleEH_CNexps.m (CNexps)
%
% M. Heinz 30Jun2004
% Views any vowel stimulus in time and spectral domain, and calculates formant trough frequencies
%

set(0,'DefaultTextInterpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur=length(stim)/Fs;
time=(1:length(stim))/Fs;

% sound(stim,Fs)

RMS=sqrt(mean(stim.^2));
RMS_TONE=1/sqrt(2);
dBreTONE=20*log10(RMS/RMS_TONE);

disp(sprintf('\nEH: Calibration'));
disp(sprintf('Vrms = %.3f  = %.1f dB re TONE',RMS,dBreTONE))

%figure(1); clf
% % % figure;
% % % subplot(211)
% % % plot(time,stim)
% % % ylabel('Amplitude')
% % % xlabel('Time (ms)')
% % % title(sprintf('Vowel time waveform (%.2f dB re: TONE)',dBreTONE))
% % % xlim([0 dur])

%%% FFT
vowel_fft=fft(stim);
Nfft=length(vowel_fft);
warning off
MagSpect=20*log10(abs(vowel_fft));
warning backtrace
MagSpect=MagSpect-max(MagSpect);  %Normalize to 0 dB at max
MagSpect(find(isinf(MagSpect)))=-500; % Set -Inf values to finite value

freqs_Hz=(0:Nfft-1)*Fs/Nfft;
ind1=find(freqs_Hz<=5000);   % Find all power in 1st 50 harmonics
ind2=find(freqs_Hz<=(Fs/2)); % All positive freqs

% %%%%%%%%%%%%%%%
% % Find features
% harm_inds=find(MagSpect(ind2)>-100);
% harmonics_Hz=freqs_Hz(harm_inds);
% %%%%%%%%% TODO HERE
% % *1) auto-find features and levels (formants are passed, so place in between harmonics if necessary
% % 2) test our needed resynthesized stimuli
% % 3) setup templates, 
% % 4) cleanup
% 
% % Set formants to actual formant freqs
% feat_inds=NaN*ones(1,8);
% for formIND=1:4
%    [y,i]=min(abs(freqs_Hz-(forms_Hz(formIND)+.001)));
%    feat_inds((formIND-1)*2+2)=i;
% end
% % Find T0
% temp_harm_inds=harm_inds(find(harmonics_Hz<freqs_Hz(feat_inds(2))));
% [y,Tind]=min(MagSpect(temp_harm_inds));
% feat_inds(1)=temp_harm_inds(Tind);
% % Find T1-T3
% for troughIND=find(isnan(feat_inds))
%    temp_harm_inds=harm_inds(find((harmonics_Hz<freqs_Hz(feat_inds(troughIND+1)))&(harmonics_Hz>freqs_Hz(feat_inds(troughIND-1)))));
%    [y,Tind]=min(MagSpect(temp_harm_inds));
%    feat_inds(troughIND)=temp_harm_inds(Tind);
% end
% feat_inds=feat_inds(1:end-1);
% 
% feat_freqs_Hz=freqs_Hz(feat_inds);
% feat_levs=MagSpect(feat_inds);
% 
% % Check if any formants not at harmonic
% for formIND=1:3
%    if isempty(find(harm_inds==feat_inds((formIND-1)*2+2)))
%       disp(sprintf('Form %d not at HARMONIC',formIND))
%       feat_levs((formIND-1)*2+2)=interp1(harmonics_Hz(max(find(harm_inds<feat_inds((formIND-1)*2+2)))+[0 1]), ...
%          MagSpect(harm_inds(max(find(harm_inds<feat_inds((formIND-1)*2+2)))+[0 1])),feat_freqs_Hz((formIND-1)*2+2));
%    end
% end
% 
% Features={'T0','F1','T1','F2','T2','F3','T3'};
% disp(sprintf('Feature\t\tFreq (Hz)\tLevel (dB re F1)'))
% for ind=1:length(feat_inds)
%    disp(sprintf('%s\t\t\t%4.f\t\t\t%4.1f',Features{ind},feat_freqs_Hz(ind),feat_levs(ind)))
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% % % subplot(212)
% % % plot(freqs_Hz(ind2)/1000,MagSpect(ind2),'g')
% % % hold on
% % % plot(freqs_Hz(ind1)/1000,MagSpect(ind1),'b')
% % % % plot(feat_freqs_Hz/1000,feat_levs,'rx')
% % % ylabel('Magnitude (dB)')
% % % xlabel('Frequency (kHz)')
% % % title('Vowel Spectrum')
% xlim([0 ceil(feat_freqs_Hz(end)/5000)*5])


%ylim([-50 0])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % TextFontSize=8;
% % % % text(.95,.93,sprintf('Fs = %.2f Hz;  F0 = %.2fHz',Fs,harmonics_Hz(1)),'Units','norm','Horiz','right','FontSize',TextFontSize)
% % % % text(.68,.83,'Feature   Freq(Hz)  Level(dB re F1)','Units','norm','Horiz','left','FontSize',TextFontSize)
% % % % for ind=1:length(feat_inds)
% % % %    text(.73,.83-.09*ind,sprintf('  %s       %4.f      %4.1f',Features{ind},feat_freqs_Hz(ind),feat_levs(ind)), ...
% % % %       'Units','norm','Horiz','left','FontSize',TextFontSize)
% % % % end

% % % hold off

% subplot(211)
% xlim([0 3/harmonics_Hz(1)])

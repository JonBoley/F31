function [feat_freqs_Hz,feat_levs_dB,dBreTONE]=getVowelParams(vowel,Fs,FULLdur,forms_Hz,viewVowel)
% File: getVowelParams.m
% From: viewvowel
%
% M. Heinz 01July2004
% Returns feature frequencies and levels, and overall dBreTONE
% takes a short version of the vowel (e.g., 1 cycle; MUST be periodic), and does analysis 
% and plotting on longer (FULLdur) version

if ~exist('viewVowel','var')
   viewVowel=1;
end
% Verify vowel is a row vector
Vsize=size(vowel);
if Vsize(1)>Vsize(2)
   vowel=vowel';
end

%% Compute Calibration re TONE
RMS=sqrt(mean(vowel.^2));
RMS_TONE=1/sqrt(2);
dBreTONE=20*log10(RMS/RMS_TONE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dur=length(vowel)/Fs;
Nreps=ceil(FULLdur/dur);  % Create full-duration waveform
FULLvowel=repmat(vowel,1,Nreps);


%%% FFT
FULLvowel_fft=fft(FULLvowel);
Nfft=length(FULLvowel_fft);
warning off  % Suppress warnings of "log of 0"
MagSpect=20*log10(abs(FULLvowel_fft));
warning backtrace
MagSpect=MagSpect-max(MagSpect);  %Normalize to 0 dB at max
MagSpect(find(isinf(MagSpect)))=-500; % Set -Inf values to finite value

freqs_Hz=(0:Nfft-1)*Fs/Nfft;
% ind1=find(freqs_Hz<=5000);   % Find all power in 1st 50 harmonics
ind2=find(freqs_Hz<=(Fs/2)); % All positive freqs

%%%%%%%%%%%%%%%
% Find features
harm_inds=find(MagSpect(ind2)>-100);
harmonics_Hz=freqs_Hz(harm_inds);

% Set formants to actual formant freqs
feat_inds=NaN*ones(1,8);
for formIND=1:4
   [y,i]=min(abs(freqs_Hz-(forms_Hz(formIND)+.001)));
   feat_inds((formIND-1)*2+2)=i;
end
% Find T0
temp_harm_inds=harm_inds(find(harmonics_Hz<freqs_Hz(feat_inds(2))));
[y,Tind]=min(MagSpect(temp_harm_inds));
if ~isempty(Tind)
   feat_inds(1)=temp_harm_inds(Tind);  % If not found, remains set at NaN
end
% Find T1-T3
for troughIND=setdiff(find(isnan(feat_inds)),1)  % Don't let it be 1, i.e., if T0 not found, it will show up as 1
   temp_harm_inds=harm_inds(find((harmonics_Hz<freqs_Hz(feat_inds(troughIND+1)))&(harmonics_Hz>freqs_Hz(feat_inds(troughIND-1)))));
   [y,Tind]=min(MagSpect(temp_harm_inds));
   if ~isempty(Tind)
      feat_inds(troughIND)=temp_harm_inds(Tind);   % If not found, remains set at NaN
   end
end
feat_inds=feat_inds(1:end-1);

feat_freqs_Hz=ones(size(feat_inds))*NaN;  % Need to handle possible NaNs for different features
feat_levs_dB=ones(size(feat_inds))*NaN;
feat_freqs_Hz(find(~isnan(feat_inds)))=freqs_Hz(feat_inds(find(~isnan(feat_inds))));
feat_levs_dB(find(~isnan(feat_inds)))=MagSpect(feat_inds(find(~isnan(feat_inds))));

% Check if any formants not at harmonic
for formIND=1:3
   if isempty(find(harm_inds==feat_inds((formIND-1)*2+2)))
      %       disp(sprintf('Form %d not at HARMONIC',formIND))
      feat_levs_dB((formIND-1)*2+2)=interp1(harmonics_Hz(max(find(harm_inds<feat_inds((formIND-1)*2+2)))+[0 1]), ...
         MagSpect(harm_inds(max(find(harm_inds<feat_inds((formIND-1)*2+2)))+[0 1])),feat_freqs_Hz((formIND-1)*2+2));
   end
end

if viewVowel
   time=(1:length(FULLvowel))/Fs; 
   sound(FULLvowel,Fs)
   
   Features={'T0','F1','T1','F2','T2','F3','T3'};
   
   figure(round(viewVowel)); clf
   set(gcf,'pos',[464   179   560   518])
   subplot(211)
   plot(time,FULLvowel)
   ylabel('Amplitude')
   xlabel('Time (ms)')
   title(sprintf('vowel time waveform (%.2f dB re: TONE)',dBreTONE))
   xlim([0 dur])
   
   subplot(212)
   plot(freqs_Hz(ind2)/1000,MagSpect(ind2),'b')
   hold on
   %    plot(freqs_Hz(ind1)/1000,MagSpect(ind1),'g')
   plot(feat_freqs_Hz/1000,feat_levs_dB,'rx')
   ylabel('Magnitude (dB)')
   xlabel('Frequency (kHz)')
   title('Vowel Spectrum')
   xlim([0 ceil(max(feat_freqs_Hz(find(~isnan(feat_freqs_Hz))))/5000)*5])
   ylim([-50 0])
   
   TextFontSize=8;
   text(.95,.93,sprintf('Fs = %.2f Hz;  F0 = %.2fHz',Fs,harmonics_Hz(1)),'Units','norm','Horiz','right','FontSize',TextFontSize)
   text(.6,.83,'Feature   Freq(Hz)  Level(dB re F1)','Units','norm','Horiz','left','FontSize',TextFontSize)
   for ind=1:length(feat_inds)
      text(.63,.83-.06*ind,sprintf('  %s       %4.f           %4.1f',Features{ind},feat_freqs_Hz(ind),feat_levs_dB(ind)), ...
         'Units','norm','Horiz','left','FontSize',TextFontSize)
   end 
   subplot(211)
   xlim([0 3/harmonics_Hz(1)])
   hold off
end
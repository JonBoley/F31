function [vowel,Fs,filename,dBreTONE,NEWforms_Hz,attenMod]=synth_BASELINE_eh_aid(TargetFreq_Hz,NEWF0_Hz,TargetFeature,Fix2Harms,mode)
% FIle: synth_BASELINE_eh_aid.m
%
% J. Boley 07Mar2012 (from synth_BASELINE_eh.m)
% M. Heinz 30Jun2004 (from makeEH.m 04Nov2003)
% Make an /eh/ like in Conley and Keilson (1995)
%
% Creates a BASELINE vowel (based on ORIG eh) with F2 placed at the TargetFrequency (TF) 
% and F0=75 Hz.  Expand BWs and formants based on shift, as would have occured with a Fs shift, 
% but here F0 is not shifted, it is set to the desired value
% 
% Also applies hearing aid gain, according to parameters supplied
%
% Original EH:
% F0=100Hz, Fs=10kHz;
%           F1     F2      F3      F4     F5
%  f (Hz)  500    1700    2500    3300   3750
% bw (Hz)   60      90     200     250    200
%
% Mode: 1: generate vowel and show plot
%       2: generate vowel only, don't show
%       3: only return filename and NEWforms_Hz
%

if ~exist('mode','var')
   mode=1;
end

ORIGforms_Hz=[500    1700    2500    3300   3750];
ORIGbws_Hz  =[60      90     200     250    200];
NUMforms=5;  % could use fewer formants
ORIGF0=100;
Fs=33000;

Features={'T0','F1','T1','F2','T2','F3','T3'};
if strcmp(TargetFeature,'F2')
   ORIGfreq=ORIGforms_Hz(2);
else
   error('Only set up for F2 as Target Feature')
end
ShiftFact=TargetFreq_Hz/ORIGfreq;

% Translate both the formant frequencies and bandwidths, as would occur with a shift in Fs
NEWforms_Hz=ORIGforms_Hz*ShiftFact;
NEWbws_Hz=ORIGbws_Hz*ShiftFact;

%%%%%%%%%%%% Fix formants to harmonics if option is set
FHtext='';
AdjustedHarms=[];
if Fix2Harms
   FHtext='_FH';
   for ind=1:NUMforms
      if rem(NEWforms_Hz(ind),NEWF0_Hz)
         AdjustedHarms=[AdjustedHarms ind];
         NEWforms_Hz(ind)=round(NEWforms_Hz(ind)/NEWF0_Hz)*NEWF0_Hz;
      end
   end
end
filename=sprintf('baseEH_F2at%05.f_F0at%03.f%s.wav',TargetFreq_Hz,NEWF0_Hz,FHtext);
         
% Mode 3 only returns info; Modes 1 and 2 return the synthesized vowel
if mode<3
   disp('***New vowel synthesized')
   if Fix2Harms
      disp(sprintf('     - Formants: %s adjusted to be at harmonics',mat2str(AdjustedHarms)))   
   end
   %% returned vowel will be 1 cycle, but any analysis and/or plotting will be done with a longer vowel
   FULLdur=0.6;   % This will be made to an integer number of periods (<=FULLdur)
   dur=1/NEWF0_Hz;
   Nreps=ceil(FULLdur/dur);  % Create full-duration waveform
   
   [time, vowel] = dovowl(NEWforms_Hz(1:NUMforms), NEWbws_Hz(1:NUMforms), NEWF0_Hz, dur, Fs);
   
   vowel=vowel'-mean(vowel);  % Remove DC
   vowel=vowel/max(vowel)*.999;  % Need amplitude to fit in wavwrite
   
   % Mode 2 uses getVowelParams to calculate params; mode 2 also will view the vowel 
   if mode==1
      viewVowel=1;
      %       FULLvowel=repmat(vowel,1,Nreps);
      %       sound(FULLvowel,Fs)
      %       disp(sprintf('\nPlaying New Synthesized Vowel'))
      %       pause(1)
      %    [vowelORIG,FsORIG]=wavread('vow17_MH10k.wav');
      %        sound(vowelORIG,FsORIG)
      %    disp('Playing Resampled ORIGINAL Vowel (EH) from CN exps')
   else
      viewVowel=0;
   end
   [feat_freqs_Hz,feat_levs_dB,dBreTONE]=getVowelParams(vowel,Fs,FULLdur,NEWforms_Hz,viewVowel);
else
   vowel=[]; Fs=[]; dBreTONE=[];
end

   
   

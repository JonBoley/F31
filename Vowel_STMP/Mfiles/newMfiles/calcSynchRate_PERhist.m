function PIC=calcSynchRate_PERhist(PIC)
% M.Heinz 13Sep2004.  Taken from PSTview.m
% Calculates Synchronized Rate vs frequency from Period histogram in PIC.
% Calls calcPERhist.m if needed.  Uses window_msec=[20,dur].
%
% Usage:PIC=calcSynchRate_PERhist(PIC)
% Input PIC: has PIC.x with stimulus parameters and spikes
% Output PIC: Stores calcs in PIC.SynchRate_PERhist, with paramaters used

% Need PERhist to calculate Sychronized Rate
if ~isfield(PIC,'PERhist')
   PIC=calcPERhist(PIC);
end

%%%% Calculate DFT from PERIOD histogram - to get Synch and phase
SynchRate_PERhist=fft(PIC.PERhist.PERhist)/length(PIC.PERhist.PERhist);
FFTfreqs=(0:length(SynchRate_PERhist)-1)*(1/PIC.PERhist.params.binWidth_sec)/length(SynchRate_PERhist);

Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
RayleighCRIT=chi2inv(1-Rayleigh_P,2);
%%%% Compute Synchs and Phases for each feature
FeatureINDs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureSynchs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeaturePhases=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureRaySig=zeros(size(PIC.FeatureFreqs_Hz));

for i=1:length(PIC.FeatureFreqs_Hz)
   %    FEATind=find(round(FFTfreqs*1000)==round(PIC.FeatureFreqs_Hz(i)*1000));
   %    FEATind=find(round(FFTfreqs)==round(PIC.FeatureFreqs_Hz(i)));  % OK to have less digits to match for PERhist, since DFT freqs are sparse
   %%%% Because the FFTfreqs are not exactly equal to the harmonics, we need to pick the closest one
   [yyy,FEATind]=min(abs(FFTfreqs-PIC.FeatureFreqs_Hz(i)));
   if ~isempty(FEATind)
      FeatureINDs(i)=FEATind;
      FeatureSynchs(i)=abs(SynchRate_PERhist(FeatureINDs(i)))/SynchRate_PERhist(1);
      FeaturePhases(i)=angle(SynchRate_PERhist(FeatureINDs(i)));
		% NOTE: for sim_FF conditions, the SCALED PERhist is used, but the
		% original number of spikes is used for statistics!!
		RayleighStat=2*PIC.PERhist.NumDrivenSpikes*FeatureSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
      if RayleighStat>RayleighCRIT
         FeatureRaySig(i)=1;
      end
   end
end

% calculate synch to each harmonic
PIC.HarmonicFreqs_Hz = PIC.FundamentalFreq_Hz*(1:floor(max(PIC.FeatureFreqs_Hz)/PIC.FundamentalFreq_Hz));
HarmonicINDs=NaN+zeros(size(PIC.HarmonicFreqs_Hz));
HarmonicSynchs=NaN+zeros(size(PIC.HarmonicFreqs_Hz));
HarmonicPhases=NaN+zeros(size(PIC.HarmonicFreqs_Hz));
HarmonicRaySig=zeros(size(PIC.HarmonicFreqs_Hz));

for i=1:length(PIC.HarmonicFreqs_Hz)
    %%%% Because the FFTfreqs are not exactly equal to the harmonics, we need to pick the closest one
   [yyy,HARMind]=min(abs(FFTfreqs-PIC.HarmonicFreqs_Hz(i)));
   if ~isempty(HARMind)
      HarmonicINDs(i)=HARMind;
      HarmonicSynchs(i)=abs(SynchRate_PERhist(HarmonicINDs(i)))/SynchRate_PERhist(1);
      HarmonicPhases(i)=angle(SynchRate_PERhist(HarmonicINDs(i)));
		% NOTE: for sim_FF conditions, the SCALED PERhist is used, but the
		% original number of spikes is used for statistics!!
		RayleighStat=2*PIC.PERhist.NumDrivenSpikes*HarmonicSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
      if RayleighStat>RayleighCRIT
         HarmonicRaySig(i)=1;
      end
   end
end

% Store calcs
PIC.SynchRate_PERhist.SynchRate_PERhist=SynchRate_PERhist;
PIC.SynchRate_PERhist.FFTfreqs=FFTfreqs;

PIC.SynchRate_PERhist.FeatureINDs=FeatureINDs;
PIC.SynchRate_PERhist.FeatureSynchs=FeatureSynchs;
PIC.SynchRate_PERhist.FeaturePhases=FeaturePhases;
PIC.SynchRate_PERhist.FeatureRaySig=FeatureRaySig;

PIC.SynchRate_PERhist.HarmonicINDs=HarmonicINDs;
PIC.SynchRate_PERhist.HarmonicSynchs=HarmonicSynchs;
PIC.SynchRate_PERhist.HarmonicPhases=HarmonicPhases;
PIC.SynchRate_PERhist.HarmonicRaySig=HarmonicRaySig;

% Store parameters used for calcs
PIC.SynchRate_PERhist.params.Rayleigh_P=Rayleigh_P;

return;

function PIC=calcSynchRate_PST(PIC)
% M.Heinz 13Sep2004.  Taken from PSTview.m
% Calculates Synchronized Rate vs frequency from PST histogram in PIC.
% Calls calcPST.m if needed.  Uses window=[20,dur].
%
% Usage:PIC=calcSynchRate_PST(PIC)
% Input PIC: has PIC.x with stimulus parameters and spikes
% Output PIC: Stores calcs in PIC.SynchRate_PST, with paramaters used


% Need PST to calculate Sychronized Rate
if ~isfield(PIC,'PST')
   PIC=calcPST(PIC);
end

SynchR_window_sec=[20 PIC.x.Hardware.Trigger.StmOn]/1000; % matches Wong et al 1998, Miller et al 1999a,b

%%%% Window for Calculating Synch Rate
drivenPSTinds = find( (PIC.PST.pst_X_sec >= SynchR_window_sec(1)) & (PIC.PST.pst_X_sec <= SynchR_window_sec(2)) );
NumDrivenSpikes=sum(PIC.PST.pst_Y_sps(drivenPSTinds)*PIC.x.Stimuli.fully_presented_lines*PIC.PST.params.binWidth_sec);

%%%% Calculate DFT from PST - to get Synch and phase
SynchRate_PST=fft(PIC.PST.pst_Y_sps(drivenPSTinds))/length(drivenPSTinds);
FFTfreqs=(0:length(SynchRate_PST)-1)*(1/PIC.PST.params.binWidth_sec)/length(SynchRate_PST);

%%%% Compute Synchs and Phases for each feature
Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
RayleighCRIT=chi2inv(1-Rayleigh_P,2);
FeatureINDs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureSynchs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeaturePhases=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureRaySig=zeros(size(PIC.FeatureFreqs_Hz));
for i=1:length(PIC.FeatureFreqs_Hz)
   %%%% Here, the FFTfreqs are not equal to the harmonics, and so we need to pick the closest one
   [yyy,FEATind]=min(abs(FFTfreqs-PIC.FeatureFreqs_Hz(i)));
   if ~isempty(FEATind)
      FeatureINDs(i)=FEATind;
      FeatureSynchs(i)=abs(SynchRate_PST(FeatureINDs(i)))/SynchRate_PST(1);
      FeaturePhases(i)=angle(SynchRate_PST(FeatureINDs(i)));
      RayleighStat=2*NumDrivenSpikes*FeatureSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
      if RayleighStat>RayleighCRIT
         FeatureRaySig(i)=1;
      end
   end
end

% Store calcs
PIC.SynchRate_PST.NumDrivenSpikes=NumDrivenSpikes;
PIC.PST.drivenPSTinds=drivenPSTinds;
PIC.SynchRate_PST.SynchRate_PST=SynchRate_PST;
PIC.SynchRate_PST.FFTfreqs=FFTfreqs;

PIC.SynchRate_PST.FeatureINDs=FeatureINDs;
PIC.SynchRate_PST.FeatureSynchs=FeatureSynchs;
PIC.SynchRate_PST.FeaturePhases=FeaturePhases;
PIC.SynchRate_PST.FeatureRaySig=FeatureRaySig;

% Store parameters used for calcs
PIC.SynchRate_PST.params.SynchR_window_sec=SynchR_window_sec;
PIC.SynchRate_PST.params.Rayleigh_P=Rayleigh_P;

return;
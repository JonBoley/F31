function [STMPunit_BF_Hz,STMPstim_FeatureTarget_Hz,STMP_SpikeTrains,NeuralDelay_sec,FreqFact] = ...
   STMP_spikeconversion(ORIGunit_BF_Hz,ORIGstim_FeatureTarget_Hz,ORIG_SpikeTrains);
% File: STMP_spikeconversion.m
% M. Heinz;   June 8, 2007
%
% FROM: newPIC=simFF_PICshift(PIC)%
%
% Converts spike times collected in the STMP scheme:
% Data collected from single unit with BForig:
%     - stimulus feature shifted above (or below) BForig by 1/FreqFact by a
%     change in sampling rate
% Predicted data:
%     - effective stimulus is centered at BForig, predicted unit has
%     BFeffective, which is shifted below (or above) BForig by FreqFact 
%     - spike times are scaled by 1/FreqFact to correct for change in
%     temporal waveform corresponding to change in stimulus sampling rate
%     by 1/FreqFact
%
% A fixed NEURAL DELAY is assumed, and is removed prior to time-scaling
% spike times.
%
% This simulation method allows the responses of a single neuron to N
% shifted stimuli to be used to simulate the responses of N neurons with
% different BFs to the same stimulus (at the nominal BF). 

NeuralDelay_sec = 0.001;   % Neural Delay to compensate for before scaling
FreqFact=ORIGunit_BF_Hz/ORIGstim_FeatureTarget_Hz;
TimeFact=1/FreqFact;

%%%% Scale all frequency parameters by x(BF/FeatureFreq)
STMPstim_FeatureTarget_Hz=ORIGstim_FeatureTarget_Hz*FreqFact;
STMPunit_BF_Hz=ORIGunit_BF_Hz*FreqFact;

%%%% Scale all spikes times by x(FeatureFreq/BF) AFTER
%%%% compensating for assumed neural delay
ORIGspikes=ORIG_SpikeTrains(:,2);
NDcompORIGspikes=ORIGspikes-NeuralDelay_sec;
NEGATIVEinds=find(NDcompORIGspikes<0);
NDcompSCALEDspikes = NDcompORIGspikes*TimeFact;
% Leave spikes before NeuralDelay ASIS since they must be spontaneous
% spikes
if ~isempty(NEGATIVEinds)
   NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
end
SCALEDspikes=NDcompSCALEDspikes+NeuralDelay_sec;
STMP_SpikeTrains(:,1)=ORIG_SpikeTrains(:,1);
STMP_SpikeTrains(:,2)=SCALEDspikes;

return;

function newPIC=simFF_PICshift(PIC)
% M.Heinz 28Sep2004.
%
% Modified 27Apr2006 - Added parameter - neural delay to take out before
% time scaling spikes
%
% converts PIC by shifting spike times and all frequencies in the opposite direction
% of where the stimulus was originally shifted relative to BF.
% - reBF data is a stimulus that has been shifted above (or below) BF.
% - simulated FF data is translated backwards, such that it models a neuron with BF shifted below (or above)
%   the nominal BF, responding to a stimulus at the nominal BF.
%
% This simulation method allows the responses of a single neuron to two shifted stimuli to be used to 
% simulate the responses of two neurons with different BFs to the same stimulus (at the nominal BF).

newPIC=PIC;

NeuralDelay_sec = 0.001;   % Neural Delay to compensate for before scaling

if strcmp(newPIC.TEMPLATE,'TrBF')
   
   FreqFact=PIC.BF_Hz/PIC.FeatureFreqs_Hz;
   TimeFact=1/FreqFact;
    
   %%%% Shift all frequency parameters by x(BF/FeatureFreq)
   newPIC.FeatureFreqs_Hz=PIC.FeatureFreqs_Hz*FreqFact;
   newPIC.FundamentalFreq_Hz=PIC.FundamentalFreq_Hz*FreqFact;   
   newPIC.BF_Hz=PIC.BF_Hz*FreqFact;
   
   newPIC.x.Stimuli.Condition.BaseFrequency_kHz=PIC.x.Stimuli.Condition.BaseFrequency_kHz*FreqFact;
   newPIC.x.Stimuli.main.tone.freq=PIC.x.Stimuli.main.tone.freq*FreqFact;   
   newPIC.x.Stimuli.main.freqs=PIC.x.Stimuli.main.freqs*FreqFact;
   
   
   %%%% Shift all time parameters by x(FeatureFreq/BF)
   newPIC.x.General.spike_res=PIC.x.General.spike_res*TimeFact;
   
   newPIC.x.Stimuli.rise_fall=PIC.x.Stimuli.rise_fall*TimeFact;
   
   newPIC.x.Hardware.Trigger.StmOn=PIC.x.Hardware.Trigger.StmOn*TimeFact;
   newPIC.x.Hardware.Trigger.StmOff=PIC.x.Hardware.Trigger.StmOff*TimeFact;
   
	%%% COMPENSATE for Neural Delay
	ORIGspikes=PIC.x.spikes{1}(:,2);
	NDcompORIGspikes=ORIGspikes-NeuralDelay_sec;
	NEGATIVEinds=find(NDcompORIGspikes<0);
	NDcompSCALEDspikes = NDcompORIGspikes*TimeFact;
	% Leave spikes before NeuralDelay ASIS since they must be spontaneous
	% spikes
	if ~isempty(NEGATIVEinds)
		NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
	end
	SCALEDspikes=NDcompSCALEDspikes+NeuralDelay_sec;
	%	OLD_SCALEDspikes = PIC.x.spikes{1}(:,2)*TimeFact;
	newPIC.x.spikes{1}(:,2)=SCALEDspikes;
   
elseif strcmp(newPIC.TEMPLATE,'EHrBF')
   
   FreqFact=PIC.BF_Hz/PIC.x.Stimuli.featureTargetFreq_Hz;
   TimeFact=1/FreqFact;
    
   %%%% Shift all frequency parameters by x(BF/FeatureFreq)
   newPIC.FeatureFreqs_Hz=PIC.FeatureFreqs_Hz*FreqFact;
   newPIC.FundamentalFreq_Hz=PIC.FundamentalFreq_Hz*FreqFact;   
   newPIC.BF_Hz=PIC.BF_Hz*FreqFact;

   newPIC.x.Stimuli.Condition.BaseFrequency_kHz=PIC.x.Stimuli.Condition.BaseFrequency_kHz*FreqFact;
   
   newPIC.x.Stimuli.BASELINE.F0_Hz=PIC.x.Stimuli.BASELINE.F0_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.TargetFreq_Hz=PIC.x.Stimuli.BASELINE.TargetFreq_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.FormFreqs_Hz=PIC.x.Stimuli.BASELINE.FormFreqs_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.Fs_Hz=PIC.x.Stimuli.BASELINE.Fs_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.FeatFreqs_Hz=PIC.x.Stimuli.BASELINE.FeatFreqs_Hz*FreqFact;

   newPIC.x.Stimuli.Computed.FeatureTarget_Hz=PIC.x.Stimuli.Computed.FeatureTarget_Hz*FreqFact;
   newPIC.x.Stimuli.Computed.UpdateRate_Hz=PIC.x.Stimuli.Computed.UpdateRate_Hz*FreqFact;
   
   newPIC.x.Stimuli.origstimUpdateRate_Hz=PIC.x.Stimuli.origstimUpdateRate_Hz*FreqFact;
   newPIC.x.Stimuli.featureFreq_Hz=PIC.x.Stimuli.featureFreq_Hz*FreqFact;
   newPIC.x.Stimuli.featureTargetFreq_Hz=PIC.x.Stimuli.featureTargetFreq_Hz*FreqFact;
   newPIC.x.Stimuli.updateRate_Hz=PIC.x.Stimuli.updateRate_Hz*FreqFact;
   newPIC.x.Stimuli.playback_sampling_rate=PIC.x.Stimuli.playback_sampling_rate*FreqFact;
   
   
   %%%% Shift all time parameters by x(FeatureFreq/BF)
   newPIC.x.General.spike_res=PIC.x.General.spike_res*TimeFact;
     
   newPIC.x.Hardware.Trigger.StmOn=PIC.x.Hardware.Trigger.StmOn*TimeFact;
   newPIC.x.Hardware.Trigger.StmOff=PIC.x.Hardware.Trigger.StmOff*TimeFact;
   
	%%% COMPENSATE for Neural Delay
	ORIGspikes=PIC.x.spikes{1}(:,2);
	NDcompORIGspikes=ORIGspikes-NeuralDelay_sec;
	NEGATIVEinds=find(NDcompORIGspikes<0);
	NDcompSCALEDspikes = NDcompORIGspikes*TimeFact;
	% Leave spikes before NeuralDelay ASIS since they must be spontaneous
	% spikes
	if ~isempty(NEGATIVEinds)
		NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
	end
	SCALEDspikes=NDcompSCALEDspikes+NeuralDelay_sec;
	%	OLD_SCALEDspikes = PIC.x.spikes{1}(:,2)*TimeFact;
	newPIC.x.spikes{1}(:,2)=SCALEDspikes;
   
elseif sum(strcmp(newPIC.TEMPLATE,{'EHrBFi','EHvNrBFi','EHvLTASSrBFi','WAVreBFi','TrBFi'}))
    
   FreqFact=PIC.BF_Hz/PIC.x.Stimuli.Used.FeatureTarget_Hz_List(PIC.CONDind);
   TimeFact=1/FreqFact;
    
   %%%% Shift all frequency parameters by x(BF/FeatureFreq)
   newPIC.FeatureFreqs_Hz=PIC.FeatureFreqs_Hz*FreqFact;
   newPIC.FundamentalFreq_Hz=PIC.FundamentalFreq_Hz*FreqFact;   
   newPIC.BF_Hz=PIC.BF_Hz*FreqFact;

   newPIC.x.Stimuli.Condition.BaseFrequency_kHz=PIC.x.Stimuli.Condition.BaseFrequency_kHz*FreqFact;
   
   newPIC.x.Stimuli.BASELINE.F0_Hz=PIC.x.Stimuli.BASELINE.F0_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.TargetFreq_Hz=PIC.x.Stimuli.BASELINE.TargetFreq_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.FormFreqs_Hz=PIC.x.Stimuli.BASELINE.FormFreqs_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.Fs_Hz=PIC.x.Stimuli.BASELINE.Fs_Hz*FreqFact;
   newPIC.x.Stimuli.BASELINE.FeatFreqs_Hz=PIC.x.Stimuli.BASELINE.FeatFreqs_Hz*FreqFact;

   newPIC.x.Stimuli.Computed.FeatureTarget_Hz_List=PIC.x.Stimuli.Computed.FeatureTarget_Hz_List*FreqFact;
   newPIC.x.Stimuli.Computed.UpdateRate_Hz_List=PIC.x.Stimuli.Computed.UpdateRate_Hz_List*FreqFact;

   newPIC.x.Stimuli.Used.FeatureTarget_Hz_List=PIC.x.Stimuli.Used.FeatureTarget_Hz_List*FreqFact;
   newPIC.x.Stimuli.Used.UpdateRate_Hz_List=PIC.x.Stimuli.Used.UpdateRate_Hz_List*FreqFact;
   newPIC.x.Stimuli.Used.OctShifts_List=PIC.x.Stimuli.Used.OctShifts_List;

   %    newPIC.x.Stimuli.origstimUpdateRate_Hz=PIC.x.Stimuli.origstimUpdateRate_Hz*FreqFact;
   %    newPIC.x.Stimuli.featureFreq_Hz=PIC.x.Stimuli.featureFreq_Hz*FreqFact;
   %    newPIC.x.Stimuli.featureTargetFreq_Hz=PIC.x.Stimuli.featureTargetFreq_Hz*FreqFact;
   newPIC.x.Stimuli.updateRate_Hz=PIC.x.Stimuli.updateRate_Hz*FreqFact;   
   %    newPIC.x.Stimuli.playback_sampling_rate=PIC.x.Stimuli.playback_sampling_rate*FreqFact;
   
   %%%% Shift all time parameters by x(FeatureFreq/BF)
   newPIC.x.General.spike_res=PIC.x.General.spike_res*TimeFact;
     
   newPIC.x.Hardware.Trigger.StmOn=PIC.x.Hardware.Trigger.StmOn*TimeFact;
   newPIC.x.Hardware.Trigger.StmOff=PIC.x.Hardware.Trigger.StmOff*TimeFact;
   
   
    % calculate NeuralDelay_sec
%     NeuralDelay_sec = 0.0032 - 0.0024*log10(PIC.BF_Hz/1e3);
%     spikeOffset_sec = NeuralDelay_sec-NeuralDelay_sec*TimeFact;
    % Greenwood Function:
    A = 163.5; k = 0.85; a = 2.1; % Chinchilla
    x = 0:0.001:1; %proportion of cochlear length
    F_Hz = A * (10.^(a*x) - k);
    [B,IX] = sort(PSTH_BF_kHz);
    dist = interp1(F_Hz,x,B*1e3);
    NeuralDelay_sec = 0.005228*dist.^2 - 0.01203*dist + 0.008404; %from clickResponse.m
    
    % apply offset based on result from findDelay.m
    % normal=0.014; impaired=0.0135; amplified=0.0115;
    NeuralDelay_sec = NeuralDelay_sec + 0.014; 
	
	%%% COMPENSATE for Neural Delay
	ORIGspikes=PIC.x.spikes{1}(:,2);
	NDcompORIGspikes=ORIGspikes-NeuralDelay_sec;
	NEGATIVEinds=find(NDcompORIGspikes<0);
	NDcompSCALEDspikes = NDcompORIGspikes*TimeFact;
	% Leave spikes before NeuralDelay ASIS since they must be spontaneous
	% spikes
	if ~isempty(NEGATIVEinds)
		NDcompSCALEDspikes(NEGATIVEinds)=NDcompORIGspikes(NEGATIVEinds);
	end
	SCALEDspikes=NDcompSCALEDspikes+NeuralDelay_sec;
	%	OLD_SCALEDspikes = PIC.x.spikes{1}(:,2)*TimeFact;
	newPIC.x.spikes{1}(:,2)=SCALEDspikes;
	  
else
   error('NOT a ''reBF'' TEMPLATE!!!!')
end

newPIC.simFF_PICshift.FreqFact=FreqFact;   
newPIC.simFF_PICshift.TimeFact=TimeFact;   
newPIC.simFF_PICshift.NeuralDelay_sec=NeuralDelay_sec;

return;
   

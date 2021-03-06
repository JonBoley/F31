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

% calculate NeuralDelay_sec
% Greenwood Function:
A = 163.5; k = 0.85; a = 2.1; % Chinchilla
x = 0:0.001:1; %proportion of cochlear length
F_Hz = A * (10.^(a*x) - k);
[B,IX] = sort(PIC.BF_Hz/1e3);%PSTH_BF_kHz);
dist = interp1(F_Hz,x,B*1e3);
% manually picked latencies (from clickResponse.m)
y_latency = [0.001	0.00116 0.0014  0.00168 0.00196 0.00232 0.00264 0.00300];
x_dist    = [0.8	0.7     0.6     0.5     0.4     0.3     0.2     0.1];
NeuralDelay_sec = interp1(x_dist,y_latency,dist);

% switch PIC.x.General.date
%     values from findDelay.m
%     case {'18-Apr-2011 ','27-Jun-2011 ','20-Jul-2011 '} % normal
%         NeuralDelay_sec = NeuralDelay_sec + 0.014;
%     case {'23-Jun-2011 ','21-Jul-2011 ','01-Aug-2011 ','09-Aug-2011 '} % impaired
%         NeuralDelay_sec = NeuralDelay_sec + 0.0135;
%     case {'04-May-2011 '} % amplified
%         NeuralDelay_sec = NeuralDelay_sec + 0.0115;
%     otherwise
%         NeuralDelay_sec = 0.001;   % Neural Delay to compensate for before scaling
% end
% fprintf('Neural delay = %1.1fms\n',NeuralDelay_sec*1e3);

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
   

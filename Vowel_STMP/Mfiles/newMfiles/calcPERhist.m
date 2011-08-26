function PIC=calcPERhist(PIC)
% M.Heinz 13Sep2004.  Taken from PSTview.m
%
% Modified: 27Apr2006 (M.Heinz)
% Adds scaling of PERhist for simFF cases.  If simFF_PICshift has been run,
% the PERhist is multiplied by TimeFact, to compensate. The NumDrivenSpikes
% is stored as the actual number after spike-time scaling, and therefore
% will not equal sum(PERhist) after PERhist is scaled.
% SCALED_NumDrivenSpikes is new value = sum(sacledPERhist)
%
%
% Calculates Period histogram from PIC.
% Uses window_msec=[20,dur].
%
% Usage:PIC=calcPERhist(PIC)
% Input PIC: has PIC.x with stimulus parameters and spikes
% Output PIC: Stores calcs in PIC.PERhist, with paramaters used

PERhist_window_sec=[20 PIC.x.Hardware.Trigger.StmOn]/1000; % matches Wong et al 1998, Miller et al 1999a,b

%%%% PERIOD histogram usage for getting Synch/Phase %%%%%%%%%%%
% K=16 used for all frequencies: estimates 1st 7 harmonics, which is not enough for vowels.
% Johnson (1980) used 30-200 bins/cycle for the PERhist;
% Anderson et al (1971) used 10 bins/cycle for PERhist.
% With K bins/cycle for all frequencies, Synch/Phase from PERhist DFT gives values at harmonics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=256;  %% # bins/cycle
binWidth_sec=1/PIC.FundamentalFreq_Hz/K;
M=floor(diff(PERhist_window_sec)*PIC.FundamentalFreq_Hz);  % Integer number of cycles to include in driven-spike window
PERhist_window_sec(2)=PERhist_window_sec(1)+M/PIC.FundamentalFreq_Hz; % Reset EndTime to limit to integer number of cycles of F0

% Find driven spikes to use
spikeTimes = PIC.x.spikes{1};
drivenSpikeIndices = find( (spikeTimes(:,2) >= PERhist_window_sec(1)) & (spikeTimes(:,2) <= PERhist_window_sec(2)) );
drivenSpikes=spikeTimes(drivenSpikeIndices,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Period Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivenSpikes_BINS=rem(floor(drivenSpikes/binWidth_sec),K)+1; %Convert times to PERhist bins (1:K)
[PERhist,xxx]=hist(drivenSpikes_BINS,(1:K));  % Make Histogram from BINS (1:K)
% This is the actual number of recorded spikes used to create this PERhist
% (use this for stats, etc ...)
NumDrivenSpikes=sum(PERhist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADJUST for temporal scaling bias (if this is a simFF condition)!!!
if isfield(PIC,'simFF_PICshift')
	PERhist = PERhist * PIC.simFF_PICshift.TimeFact;
	SCALED_NumDrivenSpikes=sum(PERhist);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Convert PERhist to spikes per second
PERhist=PERhist/PIC.x.Stimuli.fully_presented_lines/M/binWidth_sec; % Convert to sp/sec

% Store calcs
PIC.PERhist.NumDrivenSpikes=NumDrivenSpikes;
if isfield(PIC,'simFF_PICshift')
	PIC.PERhist.SCALED_NumDrivenSpikes=SCALED_NumDrivenSpikes;
end
PIC.PERhist.PERhist=PERhist;
PIC.PERhist.PERhist_X_sec=(0.5:K)*binWidth_sec;

% Store parameters used for calcs
PIC.PERhist.params.binWidth_sec=binWidth_sec;
PIC.PERhist.params.PERhist_window_sec=PERhist_window_sec;

return;

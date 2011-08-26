function [PERhist_sps,binWidth_sec,NumDrivenSpikes] = PERhist(SpikeTrains,StimDur_msec,StimF0_Hz);
% M.Heinz 09Jun2007.  Taken from PIC=calcPERhist(PIC)
%
% Calculates Period histogram from PIC.
% Uses window_msec=[20,dur].
% 
% Usage: [PH,binWidth_sec] = PERhist(SpikeTrains,StimDur_msec,StimF0_Hz);
%
% Input: SpikeTrains: 1st col: rep #; 2nd col: spike times
% Output PH: PERhistogram

IGNOREonset_msec=20;
PERhist_window_sec=[IGNOREonset_msec StimDur_msec]/1000; % matches Wong et al 1998, Miller et al 1999a,b

%%%% PERIOD histogram usage for getting Synch/Phase %%%%%%%%%%%
% K=16 used for all frequencies: estimates 1st 7 harmonics, which is not enough for vowels.
% Johnson (1980) used 30-200 bins/cycle for the PERhist;
% Anderson et al (1971) used 10 bins/cycle for PERhist.
% With K bins/cycle for all frequencies, Synch/Phase from PERhist DFT gives values at harmonics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=256;  %% # bins/cycle
binWidth_sec=1/StimF0_Hz/K;
M=floor(diff(PERhist_window_sec)*StimF0_Hz);  % Integer number of cycles to include in driven-spike window
PERhist_window_sec(2)=PERhist_window_sec(1)+M/StimF0_Hz; % Reset EndTime to limit to integer number of cycles of F0

% Find driven spikes to use
drivenSpikeIndices = find( (SpikeTrains(:,2) >= PERhist_window_sec(1)) & (SpikeTrains(:,2) <= PERhist_window_sec(2)) );
drivenSpikes=SpikeTrains(drivenSpikeIndices,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Period Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivenSpikes_BINS=rem(floor(drivenSpikes/binWidth_sec),K)+1; %Convert times to PERhist bins (1:K)
[PERhist,xxx]=hist(drivenSpikes_BINS,(1:K));  % Make Histogram from BINS (1:K)
% This is the actual number of recorded spikes used to create this PERhist
% (use this for stats, etc ...)
NumDrivenSpikes=sum(PERhist);

%%% Convert PERhist to spikes per second
PERhist_sps=PERhist/max(SpikeTrains(:,1))/M/binWidth_sec; % Convert to sp/sec

return;

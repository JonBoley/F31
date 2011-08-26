function [SynchRate_sps,DFT_FreqBinwidth_Hz]=SynchRate_PERhist(PERhist,PH_binwidth_sec)
% M.Heinz 10Jun2007.  
% FROM: function PIC=calcSynchRate_PERhist(PIC)
% Calculates Synchronized Rate vs frequency from a Period histogram.
%
% Returns: - SynchRate function
%          - Frequency binwidth

SynchRate_sps=fft(PERhist)/length(PERhist);  % scaled in spikes-per-sec

FundPeriod_sec = PH_binwidth_sec*length(PERhist);  

DFT_FreqBinwidth_Hz=1/FundPeriod_sec;

return;

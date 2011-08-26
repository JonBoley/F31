function [ALSR,INDsUsed,ALSR_oct_range]=ALSRcalc(SynchRates_sps,chan_BFs,CenterFreq)
% M.Heinz 10Jun2007.
% FROM: function [rate_sps,synch,phase_cyc,Rayleigh_P]=RSP_SynchRate(SynchRate_sps,freq,DFT_FreqBinwidth_Hz,NumSpikes);
% Calculates ALSR.
%

ALSR_oct_range=0.28;  % Slightly higher than 0.25 to allow for slight mismatches in sampling rates

INDsUsed=find((chan_BFs>=CenterFreq*2^-ALSR_oct_range)&(chan_BFs<=CenterFreq*2^ALSR_oct_range));

ALSR=mean(SynchRates_sps(INDsUsed));

return;

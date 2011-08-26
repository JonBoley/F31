function [rate_sps,synch,phase_cyc,Rayleigh_P]=RSP_SynchRate(DFTs,freq,DFT_FreqBinwidth_Hz,NumSpikes);
% M.Heinz 10Jun2007.
% FROM: function PIC=calcSynchRate_PERhist(PIC)
% Calculates Rate, Synch, and Phase (at a given frequency) from DFTs.
%

FreqCRIT=0.01;  % fraction of FreqBinwidth allowance for a match to the frequency

%%%% Calculate DFT from PERIOD histogram - to get Synch and phase
DFTfreqs=(0:length(DFTs)-1)*DFT_FreqBinwidth_Hz;

FREQind=find(abs(DFTfreqs-freq)<(FreqCRIT*DFT_FreqBinwidth_Hz));
if isempty(FREQind)
   error(sprintf('in RSP_SynchRate - requested frequency is not with %.3f Hz of a component in the DFT',FreqCRIT*DFT_FreqBinwidth_Hz))
end

Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
RayleighCRIT=chi2inv(1-Rayleigh_P,2);
%%%% Compute Synchs and Phases for each feature

rate_sps=DFTs(1);  % DC component;
synch=abs(DFTs(FREQind))/DFTs(1);
phase_cyc=angle(DFTs(FREQind))/pi;
RayleighStat=2*NumSpikes*synch^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
if RayleighStat<RayleighCRIT
   synch=NaN;   % If not significant, mark as NaN
   phase_cyc=NaN;
end

return;

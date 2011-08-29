function simTCx = simTC(freq,simBF,simthr)
%
%% File simTC.m
%
% Returns the threshold at the passed frequency, for the passed BF and threshold
% Based on simple triangular TC with tail (e.g. Goldstein, 1974)

TAILslope=2;
LFslope=15;
HFslope=40;

simTCx=(log2(freq/simBF)<-1).*(simthr-10*LFslope*log10(0.5)-10*TAILslope*log10(2*freq/simBF)) + ...  % Tail response
   ((log2(freq/simBF)>=-1)&(freq<simBF)).*(simthr-10*LFslope*log10(freq/simBF)) + ...  %% LF side
   (freq>=simBF).*(simthr+10*HFslope*log10(freq/simBF));  %% HF side



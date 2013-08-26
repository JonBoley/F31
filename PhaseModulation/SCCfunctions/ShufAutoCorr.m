function [NSAC,NSACdelays,AVGrate,TOTALspikes] = ShufAutoCorr(SpikeTrain1,DELAYbinwidth,Duration)
% File: ShufAutoCorr
% 27Sep2004: M Heinz - updated from Recruitment - NoiseCorr Analysis
% Calls: InnerSACmex.dll  (MEX-C file: InnerSACmex.c)
%
% Created: 9/2/03 M. Heinz
%
% Computes Normalized Shuffled Auto-Coprrelogram (NSAC) from a set of Spike Trains and Duration

NUMspikeREPS=length(SpikeTrain1);
%%%%%%%%%%%% Compute AVGrate
Kmax=0;
for spikeREPind=1:NUMspikeREPS
   if length(SpikeTrain1{spikeREPind})>Kmax
      Kmax=length(SpikeTrain1{spikeREPind});
   end
end

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
SpikeMAT1=NaN*ones(NUMspikeREPS,Kmax);
for REPindREF=1:NUMspikeREPS
   SpikeMAT1(REPindREF,1:length(SpikeTrain1{REPindREF}))=SpikeTrain1{REPindREF};
end
NUMspikes=sum(~isnan(SpikeMAT1(:,:))');  % Count number of real spikes in each line
TOTALspikes=sum(NUMspikes);
AVGrate=TOTALspikes/NUMspikeREPS/Duration;

%%%%%%%%%%%%% Compute Shuffled Auto-Correlogram
% SpikeMAT1 is setup in Matlab as: rows hold each spike train; 
% C/MEX: indexing goes down 1st column 1st, so we need to pass the transpose of SpikeMAT1, to get easy indexing in MEXfile
% SLOW Mfile: ints = InnerShufAutoCorr1(SpikeMAT1,NUMspikes)
%
% 27Sep04: ??? Can we speed this up eventually by not loading a full matrix (Kmax), but use 1-D vector and NumSpikes?
%

[intsMEX,TOTALints] = InnerSACmex(SpikeMAT1',NUMspikes,TOTALspikes); 
ints=intsMEX(1:TOTALints);  % Remove extra ints due to NaN's ini SpikeMAT1 matrix
clear intsMEX
ints=ints(find(~isnan(ints)))/1e-6;  % Convert into micro-seconds

%%% Normalize SAC such that no temporal correlation = 1
NSACdelays=0:DELAYbinwidth:Duration;
%NSACdelays=0:DELAYbinwidth:Duration/4;
NSACdelays=[-fliplr(NSACdelays(2:end)) NSACdelays]/1e-6;  % Convert into micro-seconds
SAC=hist(ints,NSACdelays);
NSAC=SAC/(NUMspikeREPS*(NUMspikeREPS-1)*Duration*AVGrate^2*DELAYbinwidth);  % From Louage et al (2004: J. Neurophysiol)

return;

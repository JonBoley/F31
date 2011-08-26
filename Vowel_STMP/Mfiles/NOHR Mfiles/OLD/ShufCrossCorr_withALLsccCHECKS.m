function [NSCC,NSCCdelays,AVGrate,TOTALspikes] = ShufCrossCorr(SpikeTrains,DELAYbinwidth,Duration)
% File: ShufCrossCorr
% 27Sep2004: M Heinz - updated from Recruitment - NoiseCorr Analysis
%
% Created: 9/3/03 M. Heinz
% Calls: InnerSACmex.dll  (MEX-C file: InnerSACmex.c)
%
% Computes Normalized Shuffled Cross-Coprrelogram (NSCC) from a set of Spike Trains and Duration

for UNITind=1:2
   NUMspikeREPS{UNITind}=length(SpikeTrains{UNITind});
end

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
for UNITind=1:2
   Kmax{UNITind}=0;
   for spikeREPind=1:NUMspikeREPS{UNITind}
      if length(SpikeTrains{UNITind}{spikeREPind})>Kmax{UNITind}
         Kmax{UNITind}=length(SpikeTrains{UNITind}{spikeREPind});
      end
   end
end

%%%% Compute AVGrate
for UNITind=1:2
   SpikeMAT{UNITind}=NaN*ones(NUMspikeREPS{UNITind},Kmax{UNITind});
   for REPindREF=1:NUMspikeREPS{UNITind}
      SpikeMAT{UNITind}(REPindREF,1:length(SpikeTrains{UNITind}{REPindREF}))=SpikeTrains{UNITind}{REPindREF};
   end
   NUMspikes{UNITind}=sum(~isnan(SpikeMAT{UNITind}(:,:))');  % Count number of real spikes in each line
   TOTALspikes{UNITind}=sum(NUMspikes{UNITind});   
   AVGrate{UNITind}=TOTALspikes{UNITind}/NUMspikeREPS{UNITind}/Duration;
end



%%% SLOW M-file
tic
ints1 = InnerShufCrossCorr1(SpikeMAT,NUMspikes);
disp(sprintf('InnerShufCrossCorr1: elapsed time = %.3f sec',toc))



tic
[intsMEX,TOTALints] = InnerSCCmex(SpikeMAT{1}',NUMspikes{1},TOTALspikes{1},SpikeMAT{2}',NUMspikes{2},TOTALspikes{2}); 
ints=intsMEX(1:TOTALints);  % Remove extra ints due to NaN's in SpikeMAT1 matrix
clear intsMEX
disp(sprintf('InnerSCCmex: elapsed time = %.3f sec',toc))

%%%% CHECK ALL THE SAME
if std(ints-ints1)
   warning('ints, and ints1  NOT ALL THE SAME!!!!!!!')
else
   warning('BOTH methods the SAME **********')
end





ints=ints(find(~isnan(ints)))/1e-6;  % Convert into micro-seconds

%%% Normalize SAC such that no temporal correlation = 1
NSCCdelays=0:DELAYbinwidth:Duration;
NSCCdelays=[-fliplr(NSCCdelays(2:end)) NSCCdelays]/1e-6;  % Convert into micro-seconds
SCC=hist(ints,NSCCdelays);
NSCC=SCC/(NUMspikeREPS{1}*NUMspikeREPS{2}*Duration*AVGrate{1}*AVGrate{2}*DELAYbinwidth);  % From Louage et al (2004: J. Neurophysiol)

return;

function [NSAC,NSACdelays,AVGrate,TOTALspikes,ETIMEinner] = ShufAutoCorr(SpikeTrain1,DELAYbinwidth,Duration)
% File: ShufAutoCorr
% 27Sep2004: M Heinz - updated from Recruitment - NoiseCorr Analysis
% Calls: InnerSACmex.dll  (MEX-C file: InnerSACmex.c)
%
% Created: 9/2/03 M. Heinz
%
% Computes Normalized Shuffled Auto-Coprrelogram (NSAC) from a set of Spike Trains and Duration
% Assumes all reps have spikes here, and setup Calls to take out empty spikes


NUMspikeREPS=length(SpikeTrain1);

%%%%%%%%%%%% Compute AVGrate
TOTALspikes=0;
Kmax=0;
for spikeREPind=1:NUMspikeREPS
   TOTALspikes=TOTALspikes+length(SpikeTrain1{spikeREPind});
   if length(SpikeTrain1{spikeREPind})>Kmax
      Kmax=length(SpikeTrain1{spikeREPind});
   end
end
AVGrate=TOTALspikes/NUMspikeREPS/Duration;

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
SpikeMAT1=NaN*ones(NUMspikeREPS,Kmax);
for REPindREF=1:NUMspikeREPS
   SpikeMAT1(REPindREF,1:length(SpikeTrain1{REPindREF}))=SpikeTrain1{REPindREF};
end
NUMspikes=sum(~isnan(SpikeMAT1(:,:))');  % Count number of real spikes in each line
TOTALspikes=sum(NUMspikes);


%%%%%%%%%%%%% Compute Shuffled Auto-Correlogram

tic; ints1 = InnerShufAutoCorr1(SpikeMAT1,NUMspikes); ETIMEinner1=toc;
disp(sprintf('InnerShufAutoCorr1: elapsed time = %.3f sec',ETIMEinner1))
tic; ints2 = InnerShufAutoCorr2(SpikeMAT1,NUMspikes); ETIMEinner2=toc;
disp(sprintf('InnerShufAutoCorr2: elapsed time = %.3f sec',ETIMEinner2))

% 08Nov2003: 
%           checked several examples, and InnerShufAutoCorr1&2 were identical, and InnerSACmex.
%           InnerSACmex is much much faster than both 1&2!!!
% 27Sep2004:
%           re-checked, same result!!

% SpikeMAT1 is setup in Matlab as: rows hold each spike train; 
% C/MEX: indexing goes down 1st column 1st, so we need to pass the transpose of SpikeMAT1, to get easy indexing in MEXfile
% 27Sep04: Can we speed this up eventually by not loading a full matrix (Kmax), but use 1-D vector and NumSpikes?

tic; 
[intsMEX,TOTALints] = InnerSACmex(SpikeMAT1',NUMspikes,TOTALspikes); 
ETIMEinnerMEX=toc;
ints=intsMEX(1:TOTALints);  % Remove extra ints due to NaN's ini SpikeMAT1 matrix
clear intsMEX
disp(sprintf('InnerSACmex: elapsed time = %.3f sec',ETIMEinnerMEX))
ETIMEinner=ETIMEinnerMEX;


%%%% CHECK ALL THE SAME
if std(ints-ints2)|std(ints-ints2)
   warning('ints, ints1, ints 2 NOT ALL THE SAME!!!!!!!')
else
   warning('ALL 3 methods the SAME **********')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(sprintf('NUMspikeREPS = %d;   Kmax = %d',NUMspikeREPS,Kmax))
% disp(sprintf('(1,1): %.6f; (1,last): %.6f; (end,1): %.6f; (end,last): %.6f',SpikeMAT1(1,1),SpikeMAT1(1,NUMspikes(1)),SpikeMAT1(end,1),SpikeMAT1(end,NUMspikes(end))))
% disp(sprintf('NUMspikes[%dx%d]: %s',size(NUMspikes),mat2str(NUMspikes)))
% disp(sprintf('TOTALspikes = %d',TOTALspikes))
% disp(sprintf('ints(1) = %.6f',ints(1)))

% if (mean(ints-intsMEX)|std(ints-intsMEX))
%    error('MISMATCH BETWEEN ints and intsMEX')
% end

ints=ints(find(~isnan(ints)))/1e-6;  % Convert into micro-seconds

%%% Normalize SAC such that no temporal correlation = 1
NSACdelays=0:DELAYbinwidth:Duration;
NSACdelays=[-fliplr(NSACdelays(2:end)) NSACdelays]/1e-6;  % Convert into micro-seconds
SAC=hist(ints,NSACdelays);
NSAC=SAC/(NUMspikeREPS*(NUMspikeREPS-1)*Duration*AVGrate^2*DELAYbinwidth);  % From ISH Louage et al 



return;

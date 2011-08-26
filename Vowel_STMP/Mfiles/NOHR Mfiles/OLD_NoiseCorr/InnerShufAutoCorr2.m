function ints = InnerShufAutoCorr1(SpikeMAT1,NUMspikes)
% This Matlab code is the brute force way, and is the most slow

[NUMspikeREPS,Kmax]=size(SpikeMAT1);
ints=[];

for REPindREF=1:NUMspikeREPS
   for SpikeIND=1:NUMspikes(REPindREF)
      COMPAREinds=find([1:NUMspikeREPS]~=REPindREF);  % compare to all spike trains EXCEPT current spike train
      for REPindCOMP=COMPAREinds
         ints=[ints SpikeMAT1(REPindCOMP,~isnan(SpikeMAT1(REPindCOMP,:)))-SpikeMAT1(REPindREF,SpikeIND)];   % Avoid all NaNs!!
      end
   end
end


return;
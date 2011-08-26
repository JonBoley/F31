function ints = InnerShufAutoCorr1(SpikeMAT1,NUMspikes)
% This Matlab code is speeded up by using repmat, which avoids 1 for loop, but is still slow!

[NUMspikeREPS,Kmax]=size(SpikeMAT1);
ints=[];

for REPindREF=1:NUMspikeREPS
   COMPAREinds=find([1:NUMspikeREPS]~=REPindREF);  % compare to all spike trains EXCEPT current spike train
   CompMAT=SpikeMAT1(COMPAREinds,:);
   for SpikeIND=1:NUMspikes(REPindREF)
      SpikeTimeMAT=repmat(SpikeMAT1(REPindREF,SpikeIND),size(CompMAT));
      INTS=CompMAT-SpikeTimeMAT;
      
      ints=[ints reshape(INTS',1,(NUMspikeREPS-1)*Kmax)];
   end
end

ints=ints(find(~isnan(ints)));  % Remove all NaNs!!

return;
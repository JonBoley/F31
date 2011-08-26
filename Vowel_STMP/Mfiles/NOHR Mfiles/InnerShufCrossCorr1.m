function ints = InnerShufCrossCorr1(SpikeMAT,NUMspikes)
% This Matlab code is speeded up by using repmat, which avoids 1 for loop, but is still slow!

for i=1:2
   [NUMspikeREPS{i},Kmax{i}]=size(SpikeMAT{i});
end
ints=[];

%% Use UNITind=1 as REF, UNITind=2 as COMP
for REPindREF=1:NUMspikeREPS{1}
   for SpikeIND=1:NUMspikes{1}(REPindREF)
      SpikeTimeMAT=repmat(SpikeMAT{1}(REPindREF,SpikeIND),size(SpikeMAT{2}));
      INTS=SpikeMAT{2}-SpikeTimeMAT;      
      ints=[ints reshape(INTS',1,NUMspikeREPS{2}*Kmax{2})];
   end
end

ints=ints(find(~isnan(ints)));  % Remove all NaNs!!

return;


% File: TestABR_Anal
% M Heinz 11Aug2004

ABRdata=xp.AD_Data.AD_All_V;
[Nreps,Npts]=size(ABRdata);
   
figure(4);
subplot(311)
plot(mean(ABRdata),'b')
subplot(312)
plot(mean(ABRdata(1:Nreps/4,:)),'b')
subplot(313)
plot(mean(ABRdata(Nreps/4+1:Nreps/2,:)),'b')


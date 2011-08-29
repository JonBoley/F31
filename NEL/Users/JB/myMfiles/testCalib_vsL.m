% File: testCalib_vsL.m
%
% M.Heinz (10Nov2004)
%

cdd

x20=loadpic(1);
x40=loadpic(2);
x0=loadpic(3);
x50=loadpic(4);
x10=loadpic(5);
x30=loadpic(6);

figure(1); clf
set(gcf,'pos',[459     4   560   695])
subplot(211)
semilogx(x50.CalibData(:,1),x50.CalibData(:,2),'b-')
hold on
semilogx(x40.CalibData(:,1),x40.CalibData(:,2),'r-')
semilogx(x30.CalibData(:,1),x30.CalibData(:,2),'g-')
semilogx(x20.CalibData(:,1),x20.CalibData(:,2),'k-')
semilogx(x10.CalibData(:,1),x10.CalibData(:,2),'y-')
semilogx(x0.CalibData(:,1),x0.CalibData(:,2),'c-')
hleg1=legend('50 dB','40','30','20','10','0',3)


xlim([.04 20])
subplot(212)
semilogx(x50.CalibData(:,1),x50.CalibData(:,3),'b-')
hold on
semilogx(x40.CalibData(:,1),x40.CalibData(:,3),'r-')
semilogx(x30.CalibData(:,1),x30.CalibData(:,3),'g-')
semilogx(x20.CalibData(:,1),x20.CalibData(:,3),'k-')
semilogx(x10.CalibData(:,1),x10.CalibData(:,3),'y-')
semilogx(x0.CalibData(:,1),x0.CalibData(:,3),'c-')
xlim([.04 20])
xlabel('Frequency (kHz)')
ylabel('Phase (cycles)')
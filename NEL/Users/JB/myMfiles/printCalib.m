function printCalib(picNUMs)

COLORS = {'r','b','g','k','c','m','y'};
cdd
figure(picNUMs(1))
legendTEXT='';
for i =1:length(picNUMs)
    x=loadpic(picNUMs(i))
    semilogx(x.CalibData(:,1),x.CalibData(:,2),'Color',COLORS{i})
    hold on
    legendTEXT{i}=sprintf('P%d',picNUMs(i));
end
hold off
if length(picNUMs)==1
    title(sprintf('P%d: Calib with Equalizer for %s',picNUMs,date))
else
    title(sprintf('P%d-%d: Calib with Equalizer for %s',picNUMs(1),picNUMs(end),date))
end
xlabel('Freq (kHz)')
ylabel('dB SPL')
grid
xlim([.200 20])
ylim([75 120])
set(gca,'XTickLabels',[1 10])
legend(legendTEXT,3)


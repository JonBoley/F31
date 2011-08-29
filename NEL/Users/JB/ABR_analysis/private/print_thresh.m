function [error] = print_ABR

global abr_FIG abr_Stimuli ABRmag

figure('NumberTitle','off','Name','Pretty Print Figure','Units','inches','position',[2.5 1 5 5]);
axes('Position',[.2 .2 .6 .6]);

plot(ABRmag(:,1),ABRmag(:,2),'b-','LineWidth',2);
hold on;
plot(ABRmag(:,1),ABRmag(:,3),'ro');
plot(ABRmag(:,1),ABRmag(:,4),'r--','LineWidth',2);

xlabel(get(abr_FIG.ax1.xlab,'string'),'FontSize',14);
ylabel(get(abr_FIG.ax1.ylab,'string'),'FontSize',14,'Interpreter','tex');
title(get(abr_FIG.ax1.title,'string'),'FontSize',14);
drawnow;

orient portrait;
if ispc
    print('-dwinc','-r300');
else
    print('-dpsc','-r300');
end

close('Pretty Print Figure'); 
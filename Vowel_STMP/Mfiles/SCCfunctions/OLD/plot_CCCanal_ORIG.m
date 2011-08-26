function plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsIN)
% File: plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsIN)
%
% M. Heinz May 23, 2008
% Plot 15-panel plot for SAC/SCC_CCC analysis of 4 spike trains.  Assumes
% all SAC/SCC/CCC analyiss has been done already by CCCanal.m, and the
% functions and metrics are just passed.

% User-specified plot LIMITS - can be specified from outside by including
% in paramsIN
if isfield(paramsIN,'XLIMIT_delay'), XLIMIT_delay=paramsIN.XLIMIT_delay, else XLIMIT_delay=3;   end
if isfield(paramsIN,'XLIMIT_PSDhigh'), XLIMIT_PSDhigh=paramsIN.XLIMIT_PSDhigh, else XLIMIT_PSDhigh=500;   end
if isfield(paramsIN,'YLIMIT_SClow'), YLIMIT_SClow=paramsIN.YLIMIT_SClow, else YLIMIT_SClow=0.9;   end
% Set X/Y LIMITS to be consistent across panels
YLIMIT_SAC=max([max(SACSCCfunctions.SAC_A_avg) max(SACSCCfunctions.XpAC_A_avg) max(SACSCCfunctions.SAC_B_avg) ...
	max(SACSCCfunctions.XpAC_B_avg) max(SACSCCfunctions.SCC_AB_avg) max(SACSCCfunctions.XpCC_AB_avg)]);
YLIMIT_DC=max([max(abs(SACSCCfunctions.DIFCOR_A)) max(abs(SACSCCfunctions.DIFCOR_B)) max(abs(SACSCCfunctions.DIFCOR_AB))]);
YLIMIT_SC=max([max(SACSCCfunctions.SUMCOR_A) max(SACSCCfunctions.SUMCOR_B) max(SACSCCfunctions.SUMCOR_AB)]);
YLIMIT_PSD=max([max(SACSCCfunctions.PSDsc_A) max(SACSCCfunctions.PSDsc_B) max(SACSCCfunctions.CSDsc_AB)]);

figure; clf
%% Shuffled Correlograms
subplot(531)
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SAC_A_avg,'k','LineWidth', 2.5); hold on
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.XpAC_A_avg,'k','LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
title(sprintf('CONDITION A\n'),'FontSize',14);
text(0.9,1.1,'SAC and XpAC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)

subplot(532)
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SAC_B_avg,'k','LineWidth', 2.5); hold on
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.XpAC_B_avg,'k','LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CONDITION B\n'),'FontSize',14);

subplot(533)
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SCC_AB_avg,'k','LineWidth', 2.5); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.XpCC_AB_avg,'k','LineWidth', 1);
plot(SACSCCmetrics.CDscc_usec/1000,SACSCCfunctions.SCC_AB_avg(find(SACSCCfunctions.delays_usec==SACSCCmetrics.CDscc_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CROSS-CORR (A,B)\n'),'FontSize',14);
text(0.5,1.1,'SCC and XpCC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12)
text(1,0,sprintf('CDscc=%.f usec',SACSCCmetrics.CDscc_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%%DIFCORs
subplot(5,3,4);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.DIFCOR_A,'k'); hold on;
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_A (peak=%.2f)',SACSCCmetrics.DCpeak_A),'Interpreter','none');
ylabel(sprintf('NORMALIZED\n# COINCIDENCES'),'FontSize',12)

subplot(5,3,5);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.DIFCOR_B,'k');hold on;
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_B (peak=%.2f)',SACSCCmetrics.DCpeak_B),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(5,3,6);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.DIFCOR_AB,'k'); hold on
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(SACSCCmetrics.CDtfs_usec/1000,SACSCCfunctions.DIFCOR_AB(find(SACSCCfunctions.delays_usec==SACSCCmetrics.CDtfs_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('DIFCOR_AB (peak=%.2f)',SACSCCmetrics.DCpeak_AB),'Interpreter','none');
text(1,1,sprintf('CCCtfs=%.2f',SACSCCmetrics.CCCtfs),'units','norm','VerticalAlignment','top','HorizontalAlignment','right','FontSize',12,'Color','red')
text(1,0,sprintf('CDtfs=%.f usec',SACSCCmetrics.CDtfs_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%%SUMCORs
subplot(5,3,7);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_A,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_A (adj. peak=%.2f)',SACSCCmetrics.SCpeak_A),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(5,3,8);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_B,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_B (adj. peak=%.2f)',SACSCCmetrics.SCpeak_B),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('DELAY (ms)','FontSize',12)

subplot(5,3,9);
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_AB,'k');hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(SACSCCmetrics.CDenv_usec/1000,SACSCCfunctions.SUMCOR_AB(find(SACSCCfunctions.delays_usec==SACSCCmetrics.CDenv_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k');hold off
xlim(XLIMIT_delay*[-1 1]);  ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('SUMCOR_AB (adj. peak=%.2f)',SACSCCmetrics.SCpeak_AB),'Interpreter','none');
text(1,0,sprintf('CDenv=%.f usec',SACSCCmetrics.CDenv_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%% PSDs
subplot(5,3,10);
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.PSDsc_A,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('SPECTRAL DENSITY\nAMPLITUDE'),'FontSize',12,'HorizontalAlignment','center')
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumPSDsc_A_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,11);
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.PSDsc_B,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumPSDsc_B_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,12);
%figure;
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.CSDsc_AB,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1,1,sprintf('CCCenv=%.2f',SACSCCmetrics.CCCenv),'units','norm','VerticalAlignment','top','HorizontalAlignment','right','FontSize',12,'Color','red')
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumCSDsc_AB_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

%% PSDs for random spikes (noise-bias removal)
subplot(5,3,13);
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.rand.PSDsc_A,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumPSDsc_Arand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,14);
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.rand.PSDsc_B,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('FREQUENCY (Hz)','FontSize',12)
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumPSDsc_Brand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,15);
plot(SACSCCfunctions.freqVEC,SACSCCfunctions.rand.CSDsc_AB,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([ 0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics.sums.sumCSDsc_ABrand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

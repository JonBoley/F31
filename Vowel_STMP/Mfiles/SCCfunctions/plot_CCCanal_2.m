function plot_CCCanal_2(SACSCCfunctions,SACSCCmetrics,paramsIN)
% File: plot_CCCanal(SACSCCfunctions,SACSCCmetrics,paramsIN)
%
% plot_CCCanal_0: used with CCCanal_0 - plots ANY individual rep
% plot_CCCanal_1: used with CCCanal_1 - plots 1 rep, but using RS.avgPSD/CSD
% plot_CCCanal_2: used with CCCanal_2 - plots rep 1, but also shows AVG/STD vals 
% plot_CCCanal_3: for CCCanal_3 - plots all BOOTSTRAP reps, along with all AVGs 
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
YLIMIT_SAC=max([max(SACSCCfunctions{1}.SAC_A_avg) max(SACSCCfunctions{1}.XpAC_A_avg) max(SACSCCfunctions{1}.SAC_B_avg) ...
	max(SACSCCfunctions{1}.XpAC_B_avg) max(SACSCCfunctions{1}.SCC_AB_avg) max(SACSCCfunctions{1}.XpCC_AB_avg)]);
YLIMIT_DC=max([max(abs(SACSCCfunctions{1}.DIFCOR_A)) max(abs(SACSCCfunctions{1}.DIFCOR_B)) max(abs(SACSCCfunctions{1}.DIFCOR_AB))]);
YLIMIT_SC=max([max(SACSCCfunctions{1}.SUMCOR_A) max(SACSCCfunctions{1}.SUMCOR_B) max(SACSCCfunctions{1}.SUMCOR_AB)]);
YLIMIT_PSD=max([max(SACSCCfunctions{1}.PSDsc_A) max(SACSCCfunctions{1}.PSDsc_B) max(SACSCCfunctions{1}.CSDsc_AB)]);

figure; clf
%% Shuffled Correlograms
subplot(531)
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SAC_A_avg,'k','LineWidth', 2.5); hold on
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.XpAC_A_avg,'k','LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
title(sprintf('CONDITION A\n'),'FontSize',14);
text(0.9,1.1,'SAC and XpAC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)

subplot(532)
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SAC_B_avg,'k','LineWidth', 2.5); hold on
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.XpAC_B_avg,'k','LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CONDITION B\n'),'FontSize',14);

subplot(533)
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SCC_AB_avg,'k','LineWidth', 2.5); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.XpCC_AB_avg,'k','LineWidth', 1);
plot(SACSCCmetrics{1}.CDscc_usec/1000,SACSCCfunctions{1}.SCC_AB_avg(find(SACSCCfunctions{1}.delays_usec==SACSCCmetrics{1}.CDscc_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CROSS-CORR (A,B)\n'),'FontSize',14);
text(0.5,1.1,'SCC and XpCC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12)
text(1,0,sprintf('CDscc=%.f usec',SACSCCmetrics{1}.CDscc_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%%DIFCORs
subplot(5,3,4);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.DIFCOR_A,'k'); hold on;
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_A (peak=%.2f)',SACSCCmetrics{1}.DCpeak_A),'Interpreter','none');
ylabel(sprintf('NORMALIZED\n# COINCIDENCES'),'FontSize',12)

subplot(5,3,5);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.DIFCOR_B,'k');hold on;
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_B (peak=%.2f)',SACSCCmetrics{1}.DCpeak_B),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(5,3,6);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.DIFCOR_AB,'k'); hold on
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(SACSCCmetrics{1}.CDtfs_usec/1000,SACSCCfunctions{1}.DIFCOR_AB(find(SACSCCfunctions{1}.delays_usec==SACSCCmetrics{1}.CDtfs_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('DIFCOR_AB (peak=%.2f)',SACSCCmetrics{1}.DCpeak_AB),'Interpreter','none');
text(1,1,sprintf('CCCtfs=%.2f(%.2f)',SACSCCmetrics{1}.CCCtfsAVG,SACSCCmetrics{1}.CCCtfsSTD),'units','norm','VerticalAlignment','top','HorizontalAlignment','right','FontSize',12,'Color','red')
text(1,0,sprintf('CDtfs=%.f usec',SACSCCmetrics{1}.CDtfs_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%%SUMCORs
subplot(5,3,7);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SUMCOR_A,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_A (adj. peak=%.2f)',SACSCCmetrics{1}.SCpeak_A),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(5,3,8);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SUMCOR_B,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_B (adj. peak=%.2f)',SACSCCmetrics{1}.SCpeak_B),'Interpreter','none');
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('DELAY (ms)','FontSize',12)

subplot(5,3,9);
plot(SACSCCfunctions{1}.delays_usec/1000,SACSCCfunctions{1}.SUMCOR_AB,'k');hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(SACSCCmetrics{1}.CDenv_usec/1000,SACSCCfunctions{1}.SUMCOR_AB(find(SACSCCfunctions{1}.delays_usec==SACSCCmetrics{1}.CDenv_usec)),'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k');hold off
xlim(XLIMIT_delay*[-1 1]);  ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('SUMCOR_AB (adj. peak=%.2f)',SACSCCmetrics{1}.SCpeak_AB),'Interpreter','none');
text(1,0,sprintf('CDenv=%.f usec',SACSCCmetrics{1}.CDenv_usec),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12,'Color','blue')

%% PSDs
subplot(5,3,10);
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.PSDsc_A,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('SPECTRAL DENSITY\nAMPLITUDE'),'FontSize',12,'HorizontalAlignment','center')
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumPSDsc_A_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,11);
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.PSDsc_B,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumPSDsc_B_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,12);
%figure;
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.CSDsc_AB,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1,1,sprintf('CCCenv=%.2f(%.2f)',SACSCCmetrics{1}.CCCenvAVG,SACSCCmetrics{1}.CCCenvSTD),'units','norm','VerticalAlignment','top','HorizontalAlignment','right','FontSize',12,'Color','red')
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumCSDsc_AB_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

%% PSDs for random spikes (noise-bias removal)
subplot(5,3,13);
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.rand.PSDsc_A,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumPSDsc_Arand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,14);
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.rand.PSDsc_B,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlabel('FREQUENCY (Hz)','FontSize',12)
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumPSDsc_Brand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

subplot(5,3,15);
plot(SACSCCfunctions{1}.freqVEC,SACSCCfunctions{1}.rand.CSDsc_AB,'k');  xlim([0 XLIMIT_PSDhigh]); ylim([ 0 YLIMIT_PSD]); hold on
plot(ones(1,2)*paramsIN.CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
plot(ones(1,2)*paramsIN.CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1.2*paramsIN.CCCenv_LOWmod_Hz,YLIMIT_PSD,sprintf('%.f',SACSCCmetrics{1}.sums.sumCSDsc_ABrand_CCCenv),'units','data','VerticalAlignment','top','HorizontalAlignment','left','FontSize',10)

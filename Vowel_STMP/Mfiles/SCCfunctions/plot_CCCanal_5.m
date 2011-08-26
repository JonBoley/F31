function plot_CCCanal_3(SACSCCfunctions,SACSCCmetrics,paramsIN)
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


% Find indices and Freq Ranges for SCpeak and CCCenvs TOUSE
DCpeak_TOUSEindex=find(strcmp(paramsIN.DCpeak_TOUSE,SACSCCmetrics{end}.DCpeaks_legend));
CCCtfs_TOUSEindex=find(strcmp(paramsIN.CCCtfs_TOUSE,SACSCCmetrics{end}.CCCtfs_legend));
SCpeak_TOUSEindex=find(strcmp(paramsIN.SCpeak_TOUSE,SACSCCmetrics{end}.SCpeaks_legend));
CCCenv_TOUSEindex=find(strcmp(paramsIN.CCCenv_TOUSE,SACSCCmetrics{end}.CCCenvs_legend));
DC0peak_TOUSEindex=find(strcmp('DCpeak_0',SACSCCmetrics{end}.DCpeaks_legend));
SC0peak_TOUSEindex=find(strcmp('IFFTraw_0',SACSCCmetrics{end}.SCpeaks_legend));
if ceil(CCCenv_TOUSEindex/2)>length(SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz)
	doSUMS=0;  % These are "add-ons" at the end of CCCenvs (e.g., "raw", "adjSC"), where sums from LOWmodHz to HIGHmod_Hz are undefined
else
	doSUMS=1;
end
if doSUMS
	CCCenv_LOWmod_Hz=SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz(ceil(CCCenv_TOUSEindex/2),1);
	CCCenv_HIGHmod_Hz=SACSCCmetrics{end}.sums.PSD_LHfreqs_Hz(ceil(CCCenv_TOUSEindex/2),2);
else
	CCCenv_LOWmod_Hz=NaN;  
	CCCenv_HIGHmod_Hz=NaN;
end

% Params for all reps and AVG
REPcolors={'b','r','g','c','y','m','b','r','g','c','y','m',};
AVGcolor='k';
Nreps=length(SACSCCfunctions)-1;
AVGind=length(SACSCCfunctions);

% User-specified plot LIMITS - can be specified from outside by including
% in paramsIN
if isfield(paramsIN,'XLIMIT_delay'), XLIMIT_delay=paramsIN.XLIMIT_delay;, else XLIMIT_delay=3;   end
if isfield(paramsIN,'XLIMIT_PSDhigh'), XLIMIT_PSDhigh=paramsIN.XLIMIT_PSDhigh;, else XLIMIT_PSDhigh=500;   end
if isfield(paramsIN,'YLIMIT_SClow'), YLIMIT_SClow=paramsIN.YLIMIT_SClow;, else YLIMIT_SClow=0.9;   end
% Set X/Y LIMITS to be consistent across panels
YLIMIT_SAC=-999;
for repIND=1:Nreps+1
	YLIMIT_SAC=max([YLIMIT_SAC max(SACSCCfunctions{repIND}.SAC_A_avg) max(SACSCCfunctions{repIND}.XpAC_A_avg) max(SACSCCfunctions{repIND}.SAC_B_avg) ...
		max(SACSCCfunctions{repIND}.XpAC_B_avg) max(SACSCCfunctions{repIND}.SCC_AB_avg) max(SACSCCfunctions{repIND}.XpCC_AB_avg)]);
end
YLIMIT_DC=-999;
for repIND=1:Nreps+1
	YLIMIT_DC=max([YLIMIT_DC max(abs(SACSCCfunctions{repIND}.DIFCOR_A)) max(abs(SACSCCfunctions{repIND}.DIFCOR_B)) max(abs(SACSCCfunctions{repIND}.DIFCOR_AB))]);
end
YLIMIT_SC=-999;
for repIND=1:Nreps+1
	YLIMIT_SC=max([YLIMIT_SC max(SACSCCfunctions{repIND}.SUMCORadj_A) max(SACSCCfunctions{repIND}.SUMCORadj_B) max(SACSCCfunctions{repIND}.SUMCORadj_AB)]);
end
YLIMIT_PSD=-999;
% for repIND=1:Nreps+1
% 	YLIMIT_PSD=max([YLIMIT_PSD max(SACSCCfunctions{repIND}.PSDsc_A) max(SACSCCfunctions{repIND}.PSDsc_B) max(SACSCCfunctions{repIND}.CSDsc_AB)]);
% end


figure; clf
%% Shuffled Correlograms
subplot(3,3,1)
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SAC_A_avg,'Color',REPcolors{repIND},'LineWidth', 2.5); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.XpAC_A_avg,'Color',REPcolors{repIND},'LineWidth', 1);
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SAC_A_avg,'Color',AVGcolor,'LineWidth', 2.5);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.XpAC_A_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
title(sprintf('CONDITION A\n'),'FontSize',14);
text(0.9,1.1,'SAC and XpAC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',12)
text(0.02,0.95,sprintf('%d/%d',round(SACSCCmetrics{AVGind}.NumDrivenSpikes(1,1)),round(SACSCCmetrics{AVGind}.NumDrivenSpikes(1,2))), ...
	'units','norm','FontSize',8)
text(.5,1,sprintf('CF = %.3f kHz',paramsIN.CF_A_Hz/1000),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8)

subplot(3,3,2)
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SAC_B_avg,'Color',REPcolors{repIND},'LineWidth', 2.5); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.XpAC_B_avg,'Color',REPcolors{repIND},'LineWidth', 1);
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SAC_B_avg,'Color',AVGcolor,'LineWidth', 2.5);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.XpAC_B_avg,'Color',AVGcolor,'LineWidth', 1);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1])
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CONDITION B\n'),'FontSize',14);
text(0.02,0.95,sprintf('%d/%d',round(SACSCCmetrics{AVGind}.NumDrivenSpikes(2,1)),round(SACSCCmetrics{AVGind}.NumDrivenSpikes(2,2))), ...
	'units','norm','FontSize',8)
text(.5,1,sprintf('CF = %.3f kHz',paramsIN.CF_B_Hz/1000),'units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',8)

subplot(3,3,3)
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SCC_AB_avg,'Color',REPcolors{repIND},'LineWidth', 2.5); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.XpCC_AB_avg,'Color',REPcolors{repIND},'LineWidth', 1);
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SCC_AB_avg,'Color',AVGcolor,'LineWidth', 2.5);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.XpCC_AB_avg,'Color',AVGcolor,'LineWidth', 1);
plot(SACSCCmetrics{AVGind}.CDscc_usec/1000,SACSCCfunctions{AVGind}.SCC_AB_avg(find(SACSCCfunctions{AVGind}.delays_usec==SACSCCmetrics{AVGind}.CDscc_usec)), ...
	'bx','MarkerSize',10,'LineWidth',2)
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SAC*[0 1],'--k');
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_SAC*[0 1]); hold off
set(gca, 'Box', 'off', 'TickDir', 'out');
title(sprintf('CROSS-CORR (A,B)\n'),'FontSize',14);
text(0.5,1.1,'SCC and XpCC','units','norm','VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',12)
text(1,0,sprintf('CDscc=%.f usec',SACSCCmetrics{AVGind}.CDscc_usec),'units','norm','VerticalAlignment','bottom', ...
	'HorizontalAlignment','right','FontSize',10,'Color','blue')

%%DIFCORs
subplot(3,3,4);
DCpeak_A_TEXT='[';
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.DIFCOR_A,'Color',REPcolors{repIND}); hold on
	DCpeak_A_TEXT=sprintf('%s%.2f, ',DCpeak_A_TEXT,SACSCCmetrics{repIND}.DCpeak_A(DC0peak_TOUSEindex));   % ACF always at ZERO
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.DIFCOR_A,'Color',AVGcolor,'Linewidth',2);
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_A (peak=%.2f %s] ["%s"] )',SACSCCmetrics{AVGind}.DCpeak_A(DC0peak_TOUSEindex),DCpeak_A_TEXT(1:end-2),paramsIN.DCpeak_TOUSE), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');
ylabel(sprintf('NORMALIZED\n# COINCIDENCES'),'FontSize',12)

subplot(3,3,5);
DCpeak_B_TEXT='[';
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.DIFCOR_B,'Color',REPcolors{repIND}); hold on
	DCpeak_B_TEXT=sprintf('%s%.2f, ',DCpeak_B_TEXT,SACSCCmetrics{repIND}.DCpeak_B(DC0peak_TOUSEindex));  % ACF always at ZERO
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.DIFCOR_B,'Color',AVGcolor,'Linewidth',2);
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
title(sprintf('DIFCOR_B (peak=%.2f %s])',SACSCCmetrics{AVGind}.DCpeak_B(DC0peak_TOUSEindex),DCpeak_B_TEXT(1:end-2)),'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(3,3,6);
DCpeak_AB_TEXT='[';
CCCtfsTEXT='[';
for repIND=1:Nreps
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.DIFCOR_AB,'Color',REPcolors{repIND}); hold on
	DCpeak_AB_TEXT=sprintf('%s%.2f, ',DCpeak_AB_TEXT,SACSCCmetrics{repIND}.DCpeak_AB(DCpeak_TOUSEindex));
	CCCtfsTEXT=sprintf('%s%.2f, ',CCCtfsTEXT,SACSCCmetrics{repIND}.CCCtfs(CCCtfs_TOUSEindex));
end
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.DIFCOR_AB,'Color',AVGcolor,'Linewidth',2);
plot(SACSCCmetrics{AVGind}.CDtfs_usec/1000,SACSCCfunctions{AVGind}.DIFCOR_AB(find(SACSCCfunctions{AVGind}.delays_usec==SACSCCmetrics{AVGind}.CDtfs_usec)), ...
	'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_DC*[-1 1],'--k');
plot(XLIMIT_delay*[-1 1],zeros(1,2),'-k');hold off
xlim(XLIMIT_delay*[-1 1]); ylim(YLIMIT_DC*[-1 1]);
set(gca, 'Box', 'off', 'TickDir', 'out');
% title(sprintf('DIFCOR_AB (peak=%.2f %s])',SACSCCmetrics{AVGind}.DCpeak_AB,DCpeak_AB_TEXT(1:end-2)),'Interpreter','none','FontSize',8);
% text(0.05,1,sprintf('%s=%.2f %s; %.2f(%.2f)]',texlabel('rho_{TFS}'),SACSCCmetrics{AVGind}.CCCtfs,CCCtfsTEXT(1:end-2), ...
% 	SACSCCmetrics{AVGind}.CCCtfsAVG,SACSCCmetrics{AVGind}.CCCtfsSTD),'units','norm', ...
% 	'VerticalAlignment','top','HorizontalAlignment','left','FontSize',8,'Color','red')
title(sprintf('%s=%.2f %s; %.2f(%.2f)]',texlabel('rho_{TFS}'),SACSCCmetrics{AVGind}.CCCtfs(CCCtfs_TOUSEindex),CCCtfsTEXT(1:end-2), ...
	SACSCCmetrics{AVGind}.CCCtfsAVGs(CCCtfs_TOUSEindex),SACSCCmetrics{AVGind}.CCCtfsSTDs(CCCtfs_TOUSEindex)),'FontSize',10,'Color','red');
text(0.02,1,sprintf('DIFCOR_AB (peak=%.2f %s])',SACSCCmetrics{AVGind}.DCpeak_AB(DCpeak_TOUSEindex),DCpeak_AB_TEXT(1:end-2)),'units','norm', ...
	'VerticalAlignment','top','HorizontalAlignment','left','FontSize',8,'Interpreter','none')
text(1,0,sprintf('CDtfs=%.f usec',SACSCCmetrics{AVGind}.CDtfs_usec),'units','norm','VerticalAlignment','bottom', ...
	'HorizontalAlignment','right','FontSize',10,'Color','blue')
text(.02,.02,sprintf('["%s"]',paramsIN.CCCtfs_TOUSE),'Interpreter','none','HorizontalAlignment','left', ...
	'VerticalAlignment','bottom','FontSize',8,'Color','red','Units','norm')

%%SUMCORs
subplot(3,3,7);
SCpeak_A_TEXT='[';
for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCOR_A,'Color',REPcolors{repIND}); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCORadj_A,'Color',REPcolors{repIND}); hold on
	SCpeak_A_TEXT=sprintf('%s%.2f, ',SCpeak_A_TEXT,SACSCCmetrics{repIND}.SCpeaks_A(SC0peak_TOUSEindex));  % ACF always at ZERO
end
% plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCOR_A,'Color',AVGcolor,'Linewidth',2);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCORadj_A,'Color',AVGcolor,'Linewidth',2);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_A (peak=%.2f %s] ["%s"] )',SACSCCmetrics{AVGind}.SCpeaks_A(SC0peak_TOUSEindex),SCpeak_A_TEXT(1:end-2),paramsIN.SCpeak_TOUSE), ...
	'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(3,3,8);
SCpeak_B_TEXT='[';
for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCOR_B,'Color',REPcolors{repIND}); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCORadj_B,'Color',REPcolors{repIND}); hold on
	SCpeak_B_TEXT=sprintf('%s%.2f, ',SCpeak_B_TEXT,SACSCCmetrics{repIND}.SCpeaks_B(SC0peak_TOUSEindex));  % ACF always at ZERO
end
% plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCOR_B,'Color',AVGcolor,'Linewidth',2);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCORadj_B,'Color',AVGcolor,'Linewidth',2);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
title(sprintf('SUMCOR_B (peak=%.2f %s])',SACSCCmetrics{AVGind}.SCpeaks_B(SC0peak_TOUSEindex),SCpeak_B_TEXT(1:end-2)),'Interpreter','none','FontSize',8);
set(gca, 'Box', 'off', 'TickDir', 'out');

subplot(3,3,9);
SCpeak_AB_TEXT='[';
CCCenvTEXT='[';
for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCOR_AB,'Color',REPcolors{repIND}); hold on
	plot(SACSCCfunctions{repIND}.delays_usec/1000,SACSCCfunctions{repIND}.SUMCORadj_AB,'Color',REPcolors{repIND}); hold on
	SCpeak_AB_TEXT=sprintf('%s%.2f, ',SCpeak_AB_TEXT,SACSCCmetrics{repIND}.SCpeaks_AB(SCpeak_TOUSEindex));
	CCCenvTEXT=sprintf('%s%.2f, ',CCCenvTEXT,SACSCCmetrics{repIND}.CCCenvs(CCCenv_TOUSEindex));
end
% plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCOR_AB,'Color',AVGcolor,'Linewidth',2);
plot(SACSCCfunctions{AVGind}.delays_usec/1000,SACSCCfunctions{AVGind}.SUMCORadj_AB,'Color',AVGcolor,'Linewidth',2);
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(SACSCCmetrics{AVGind}.CDenv_usec/1000,SACSCCfunctions{AVGind}.SUMCORadj_AB(find(SACSCCfunctions{AVGind}.delays_usec==SACSCCmetrics{AVGind}.CDenv_usec)), ...
	'bx','MarkerSize',10,'LineWidth',2)
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k');hold off
xlim(XLIMIT_delay*[-1 1]);  ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
text(1,0,sprintf('CDenv=%.f usec',SACSCCmetrics{AVGind}.CDenv_usec), ...
	'units','norm','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',10,'Color','blue')

title(sprintf('%s=%.2f %s; %.2f(%.2f)]',texlabel('rho_{ENV}'),SACSCCmetrics{AVGind}.CCCenvs(CCCenv_TOUSEindex),CCCenvTEXT(1:end-2), ...
	SACSCCmetrics{AVGind}.CCCenvAVGs(CCCenv_TOUSEindex),SACSCCmetrics{AVGind}.CCCenvSTDs(CCCenv_TOUSEindex)),'FontSize',10,'Color','red');
text(0.02,1,sprintf('SUMCOR_AB (peak=%.2f %s])',SACSCCmetrics{AVGind}.SCpeaks_AB(SCpeak_TOUSEindex),SCpeak_AB_TEXT(1:end-2)),'units','norm', ...
	'VerticalAlignment','top','HorizontalAlignment','left','FontSize',8,'Interpreter','none')
text(.02,.02,sprintf('["%s"]',paramsIN.CCCenv_TOUSE),'Interpreter','none','HorizontalAlignment','left', ...
	'VerticalAlignment','bottom','FontSize',8,'Color','red','Units','norm')


% title(sprintf('["%s"]',paramsIN.CCCenv_TOUSE),'Interpreter','none','HorizontalAlignment','center','FontSize',8,'Color','red')
% title(sprintf('%s=%.2f %s; %.2f(%.2f)]',texlabel('rho_{ENV}'),SACSCCmetrics{AVGind}.CCCenvs(CCCenv_TOUSEindex),CCCenvTEXT(1:end-2), ...
% 	SACSCCmetrics{AVGind}.CCCenvAVGs(CCCenv_TOUSEindex),SACSCCmetrics{AVGind}.CCCenvSTDs(CCCenv_TOUSEindex)),'FontSize',10,'Color','red')



% %% PSDs
% subplot(5,3,10);
% sumPSDsc_A_CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.PSDsc_A,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumPSDsc_A_CCCenvTEXT=sprintf('%s%.f, ',sumPSDsc_A_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumPSD_A(ceil(CCCenv_TOUSEindex/2)));
% 	end
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.PSDsc_A,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% ylabel(sprintf('SPECTRAL DENSITY\nAMPLITUDE'),'FontSize',12,'HorizontalAlignment','center')
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumPSD_A(ceil(CCCenv_TOUSEindex/2)),sumPSDsc_A_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left'
% 		,'FontSize',8)
% end
% 
% subplot(5,3,11);
% sumPSDsc_B_CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.PSDsc_B,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumPSDsc_B_CCCenvTEXT=sprintf('%s%.f, ',sumPSDsc_B_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumPSD_B(ceil(CCCenv_TOUSEindex/2)));
% 	end
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.PSDsc_B,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumPSD_B(ceil(CCCenv_TOUSEindex/2)),sumPSDsc_B_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left','FontSize',8)
% end
% title(sprintf('["%s"]',paramsIN.CCCenv_TOUSE),'Interpreter','none','HorizontalAlignment','center','FontSize',8,'Color','red')
% 
% subplot(5,3,12);
% sumCSDsc_AB_CCCenvTEXT='[';
% CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.CSDsc_AB,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumCSDsc_AB_CCCenvTEXT=sprintf('%s%.f, ',sumCSDsc_AB_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumCSD_AB(ceil(CCCenv_TOUSEindex/2)));
% 	end
% 	CCCenvTEXT=sprintf('%s%.2f, ',CCCenvTEXT,SACSCCmetrics{repIND}.CCCenvs(CCCenv_TOUSEindex));
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.CSDsc_AB,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumCSD_AB(ceil(CCCenv_TOUSEindex/2)),sumCSDsc_AB_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left','FontSize',8)
% end
% title(sprintf('%s=%.2f %s; %.2f(%.2f)]',texlabel('rho_{ENV}'),SACSCCmetrics{AVGind}.CCCenvs(CCCenv_TOUSEindex),CCCenvTEXT(1:end-2), ...
% 	SACSCCmetrics{AVGind}.CCCenvAVGs(CCCenv_TOUSEindex),SACSCCmetrics{AVGind}.CCCenvSTDs(CCCenv_TOUSEindex)),'FontSize',10,'Color','red')
% 
% 
% %% PSDs for random spikes (noise-bias removal)
% subplot(5,3,13);
% sumPSDsc_Arand_CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.rand.PSDsc_A,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumPSDsc_Arand_CCCenvTEXT=sprintf('%s%.f, ',sumPSDsc_Arand_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumPSDrand_A(ceil(CCCenv_TOUSEindex/2)));
% 	end
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.rand.PSDsc_A,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumPSDrand_A(ceil(CCCenv_TOUSEindex/2)),sumPSDsc_Arand_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left','FontSize',8)
% end
% 
% subplot(5,3,14);
% sumPSDsc_Brand_CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.rand.PSDsc_B,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumPSDsc_Brand_CCCenvTEXT=sprintf('%s%.f, ',sumPSDsc_Brand_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumPSDrand_B(ceil(CCCenv_TOUSEindex/2)));
% 	end
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.rand.PSDsc_B,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% xlabel('FREQUENCY (Hz)','FontSize',12)
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumPSDrand_B(ceil(CCCenv_TOUSEindex/2)),sumPSDsc_Brand_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left','FontSize',8)
% end
% 
% subplot(5,3,15);
% sumCSDsc_ABrand_CCCenvTEXT='[';
% for repIND=1:Nreps
% 	plot(SACSCCfunctions{repIND}.freqVEC,SACSCCfunctions{repIND}.rand.CSDsc_AB,'Color',REPcolors{repIND});  hold on
% 	if doSUMS
% 		sumCSDsc_ABrand_CCCenvTEXT=sprintf('%s%.f, ',sumCSDsc_ABrand_CCCenvTEXT,SACSCCmetrics{repIND}.sums.sumCSDrand_AB(ceil(CCCenv_TOUSEindex/2)));
% 	end
% end
% plot(SACSCCfunctions{AVGind}.freqVEC,SACSCCfunctions{AVGind}.rand.CSDsc_AB,'Color',AVGcolor,'Linewidth',2);  xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]);
% plot(ones(1,2)*CCCenv_LOWmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2);
% plot(ones(1,2)*CCCenv_HIGHmod_Hz,[0 YLIMIT_PSD],'--k','linewidth',2); hold off
% xlim([0 XLIMIT_PSDhigh]); ylim([0 YLIMIT_PSD]); 
% set(gca, 'Box', 'off', 'TickDir', 'out');
% if doSUMS
% 	text(0.025,1,sprintf('%.f %s]',SACSCCmetrics{AVGind}.sums.sumCSDrand_AB(ceil(CCCenv_TOUSEindex/2)),sumCSDsc_ABrand_CCCenvTEXT(1:end-2)), ...
% 		'units','norm','VerticalAlignment','top','HorizontalAlignment','left','FontSize',8)
% end

orient landscape

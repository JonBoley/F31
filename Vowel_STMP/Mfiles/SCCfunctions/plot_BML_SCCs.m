function plot_BML_SCCs(DATAfilename)
% File: plot_BML_SCCs

%% SAVE BML Data
disp(sprintf('LOADING BML Data: %s_BML.mat',DATAfilename))
eval(['load ' DATAfilename sprintf('_BML.mat')])


%% Plot figures
figure
LEGfont=6;
subplot(511)
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.rate_A,'b-'); hold on
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.rate_B,'r-'); hold off
ylabel('Rate (sps)')
h1=legend('A','B','Location','best');
set(h1,'FontSize',LEGfont)
title(DATAfilename,'Interpreter','none')
subplot(512)
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.DCpeak_A,'b-'); hold on
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.DCpeak_B,'r-'); 
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.DCpeak_AB,'k-'); hold off
ylabel('DC peak height')
h2=legend('A','B','AB','Location','best');
set(h2,'FontSize',LEGfont)
subplot(513)
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.SCpeak_A,'b-'); hold on
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.SCpeak_B,'r-'); 
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.SCpeak_AB,'k-'); hold off
ylabel('SC peak height (Adj)')
h3=legend('A','B','AB','Location','best');
set(h3,'FontSize',LEGfont)
subplot(514)
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.CCCtfs,'k-'); 
ylabel(texlabel('rho_{TFS}'))
ylim([0 1])
set(gca,'YTick',[0:.2:1])
grid on 
subplot(515)
plot(BMLdata.OALevel_dBSPL_VEC,BMLdata.CCCenv,'k-');
ylabel(texlabel('rho_{ENV}'))
xlabel('Overall Level (dB SPL)')
ylim([0 1])
set(gca,'YTick',[0:.2:1])
grid on 
orient tall


eval(['saveas(gcf,''' DATAfilename ''',''fig'')'])
eval(['saveas(gcf,''' DATAfilename ''',''pdf'')'])

%%
XLIMIT_delay=3;
YLIMIT_SClow=0.9;

load ihc_test_pos20.mat
YLIMIT_SC=max([max(SACSCCfunctions.SUMCOR_A) max(SACSCCfunctions.SUMCOR_B) max(SACSCCfunctions.SUMCOR_AB)]);

subplot(2,3,1), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_B,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Normal)');

subplot(2,3,2), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_A,'b'); hold on
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Impaired)');

subplot(2,3,3), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_AB,'b'); hold on
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Across Condition)');

load mixed_test_pos20.mat

subplot(2,3,2), hold on;
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_A,'r'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Impaired)');
legend('IHC','Mixed');
hold off;

subplot(2,3,3), hold on;
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_AB,'r'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Across Condition)');
legend('IHC','Mixed');
hold off;


load ihc_test_neg20.mat

subplot(2,3,4), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_B,'k'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Normal)');

subplot(2,3,5), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_A,'b'); hold on
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Impaired)');

subplot(2,3,6), plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_AB,'b'); hold on
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Across Condition)');

load mixed_test_neg20.mat

subplot(2,3,5), hold on;
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_A,'r'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Impaired)');
legend('IHC','Mixed');
hold off;

subplot(2,3,6), hold on;
plot(SACSCCfunctions.delays_usec/1000,SACSCCfunctions.SUMCOR_AB,'r'); hold on
plot(XLIMIT_delay*[-1 1],ones(1,2),'-k');
plot(zeros(1,2),YLIMIT_SC*[0 1],'--k'); hold off
xlim(XLIMIT_delay*[-1 1]); ylim([YLIMIT_SClow YLIMIT_SC]);
set(gca, 'Box', 'off', 'TickDir', 'out');
title('SumCor (Across Condition)');
legend('IHC','Mixed');
hold off;

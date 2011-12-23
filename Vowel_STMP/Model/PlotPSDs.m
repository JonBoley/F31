figure;
PSD_A = abs(fft(SACSCCfunctions.SAC_A_avg));
PSD_freqs=1./(min(setdiff(abs(SACSCCfunctions.delays_usec*1e-6),0)):...
    min(setdiff(abs(SACSCCfunctions.delays_usec*1e-6),0)):...
    max(SACSCCfunctions.delays_usec*1e-6));
PSD_freqs=fliplr(PSD_freqs);
PSD_freqs = min(PSD_freqs)/2:min(PSD_freqs)/2:max(PSD_freqs)/2;
semilogx(PSD_freqs,20*log10(trifilt(PSD_A(1:length(PSD_freqs)),5)),'b:');

hold on;
PSD_B = abs(fft(SACSCCfunctions.SAC_B_avg));
semilogx(PSD_freqs,20*log10(trifilt(PSD_B(1:length(PSD_freqs)),5)),'b--');

CSD = abs(fft(SACSCCfunctions.SCC_AB_avg));
semilogx(PSD_freqs,20*log10(trifilt(CSD(1:length(PSD_freqs)),5)),'k');

loglog([500 500],[1 1e2],'b:');
loglog([1700 1700],[1 1e2],'b:');
loglog([2500 2500],[1 1e2],'b:');
% axis([200 5000 1 1e3]);
legend('PSD (F3)','PSD (2.06kHz)','CSD')
title('Normal hearing, -6dB SNR, 2.06kHz CF re F3')


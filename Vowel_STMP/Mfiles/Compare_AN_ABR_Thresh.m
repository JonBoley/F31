%% compare ABR and AN threshold shifts

switch lower(getenv('computername'))
    case {'m4500'}
        addpath('C:\Research\MATLAB\Vowel_STMP\Mfiles');
        addpath('C:\Research\MATLAB\NEL\Users\JB\ABR_analysis');
    otherwise %e.g., {'1353-neurophys'}
        addpath('C:\JonBoley\Vowel_STMP');
        addpath('C:\NEL\Users\JB\ABR_analysis');
end

PlotAllExpThresholds;
plot_abr_threshes;
% AN results: [f;-audiogram] and [f;-audiogram2]
% ABR results: [freqs_all'/1e3;nanmean(thresh_shift_all,2)']

abr_thresh_shift_avg = nanmean(thresh_shift_all');
[f_khz,ia,ib] = intersect(f,freqs_all'/1e3);

figure, plot(0:100,'k:'); hold on;
scatter(-audiogram(ia),abr_thresh_shift_avg(ib),50,f_khz);
scatter(-audiogram2(ia),abr_thresh_shift_avg(ib),25,f_khz);
hold off;
MaxShift = max(max(-audiogram(ia)),max(abr_thresh_shift_avg(ib)));
MaxShift = max(MaxShift,max(-audiogram(ia)));
axis([0 MaxShift+10 0 MaxShift+10]);
xlabel('AN Threshold'); ylabel('ABR Threshold');
colorbar; colormap(jet(2*length(f_khz)));
title('AN/ABR threshold comparison. Frequency(kHz) indicated by color')

figure, hold on;
plot(freqs_all,-thresh_shift_all,'k*-'); 
h_abr=plot(freqs_all,-nanmean(thresh_shift_all')','k-','LineWidth',3); 
h_an=plot(f*1e3,audiogram,'go-','LineWidth',3);
h_an2=plot(f*1e3,audiogram2,'go:','LineWidth',3);
h_an3=plot(f*1e3,(audiogram+audiogram2)/2,'bo-','LineWidth',3);
hold off;
legend([h_abr;h_an;h_an2;h_an3],['ABR    ';'AN(avg)';'AN(min)';'AN-both']);
title('ABR Threshold Shifts');
xlabel('frequency (Hz)');
ylabel('threshold shift (dB)');
set(gca,'XScale','log'); xlim([min(freqs_all) max(freqs_all)]);

% dock all windows
set(findobj('WindowStyle','normal'),'WindowStyle','docked')


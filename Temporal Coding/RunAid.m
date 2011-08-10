function real_aid=RunAid(tgt,Fs,FREQUENCIES,NAL_IG)

orig_len = length(tgt);

DIAGNOSTIC=0;
FIR_ORDER = 64;
edges_1k = [ 1.1890, 1.6818, 2.3781, 3.3636, 4.7561, 6.7272, 8.0000]*1e3;  %% edges for compression channels
chan_crs_1k = [ ... %% half-octave, overlapped, designed using Camfit 12-01-07
    %% slope = 8dB/oct %% *** NEW *** higher slope 22-01-07
    1.23 1.41 1.53 1.42 1.61 1.67; %% slope,gains = 8,12

    1.23 1.41 1.53 1.55 1.81 1.90; %% slope,gains = 8,18
    1.23 1.41 1.53 1.55 1.91 2.03; %% slope,gains = 8,24 **** FULL GAIN ONLY 20 dB @6k ****
    %% slope = 13dB/oct %% *** NEW *** higher slope 22-01-07
    1.34 1.55 1.55 1.42 1.61 1.67; %% slope,gains = 13,12
    1.34 1.55 1.68 1.55 1.81 1.90; %% slope,gains = 13,18
    1.34 1.55 1.68 1.80 2.09 2.19; %% slope,gains = 13,24
    %% slope = 18dB/oct %% UNCHANGED
    1.42 1.55 1.55 1.42 1.61 1.67; %% slope,gains = 18,12
    1.42 1.70 1.70 1.55 1.81 1.90; %% slope,gains = 18,18
    1.42 1.70 1.86 1.80 2.22 2.38; %% slope,gains = 18,24
    ];
edges_2k = [1.6818, 2.3781, 3.3636, 4.7561, 6.7272, 8.0000]*1e3; %%  %% edges for compression channels
chan_crs_2k = [ ... %% half-octave spacing, designed using Camfit 12-01-07
    %% slope = 8dB/oct
    1.04 1.14 1.31 1.58 1.67; %% slope,gains = 8,12
    1.04 1.24 1.36 1.59 1.67; %% slope,gains = 8,18 **** FULL GAIN NOT ACHIEVABLE ****
    1.04 1.24 1.36 1.59 1.67; %% slope,gains = 8,24 **** FULL GAIN NOT ACHIEVABLE ****
    %% slope = 13dB/oct (no real change from above either, despite slope change)
    1.04 1.34 1.42 1.61 1.67; %% slope,gains = 13,12
    1.04 1.34 1.42 1.77 1.90; %% slope,gains = 13,18
    1.04 1.34 1.42 1.87 2.03; %% slope,gains = 13,24 *** FULL GAIN NOT ACHIEVABLE (only @ 6k) ***
    %% slope = 18dB/oct UNCHANGED
    1.04 1.38 1.42 1.61 1.67; %% slope,gains = 18,12
    1.04 1.38 1.55 1.81 1.90; %% slope,gains = 18,18
    1.04 1.38 1.63 2.17 2.38; %% slope,gains = 18,24
    ];

relative_level = 24; % max prescription dB
slope = 18;% dB/octave
REF_RMSdB = -22.5; %% RMS of signal file to be processed
deltaCTdB = -15; %% how far compression threshold is below channel RMS. Ccompression thresholds always the same relative to RMS channel powers
nse_burst = gen_sii_noise(2, Fs); %% 2 seconds
nse_burst = nse_burst * (10^(.05*REF_RMSdB)) / (sqrt(sum(nse_burst.^2)./length(nse_burst)));
% set up compression thresholds, in linear domain, for each edge frequency, no_eq
[dB_band_lvls_1k, band_cfs] = band_levels(nse_burst, Fs, edges_1k);
lin_thresh_1k_pre = 10.^(.05*(REF_RMSdB + dB_band_lvls_1k + deltaCTdB));
[dB_band_lvls_2k, band_cfs] = band_levels(nse_burst, Fs, edges_2k);
lin_thresh_2k_pre = 10.^(.05*(REF_RMSdB + dB_band_lvls_2k + deltaCTdB));
% so also pad around noise burst with 50 msec of silence
silence = zeros(1,(50*Fs/1000));
% insert impulse between silences: should be decayed by noise burst
tgt = [silence .3162 silence nse_burst silence tgt']; %% nse_burst is 2 secs, silence 50msec

eval(sprintf('cmpr_tgt = NchanFbankAid(tgt, edges_%dk, Fs, chan_crs_%dk(%d,:), lin_thresh_%dk_pre);', 1, 1, 1, 1));
%function aperture = prescription_design_function(fref, hf_dB, slope, filt_order, fs, diagnostic);
aperture = prescription_design_function(1e3, relative_level, slope, FIR_ORDER, Fs, DIAGNOSTIC,FREQUENCIES,NAL_IG);  %% MAStone designed filter for DEMO only
real_aid = filter(aperture, 1, cmpr_tgt); %% apply insertion gain
clear cmpr_tgt

real_aid = real_aid(end-orig_len:end);
% subplot(2,1,1),plot(tgt(end-orig_len:end));
% subplot(2,1,2),plot(real_aid);

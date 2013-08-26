% impairment initialization
% mild hearing loss (Bruce, ISAAR 2007)
Audiogram_freq = [250 500 1000 2000 3000 4000 6000];
Audiogram = [20 20 30 40 45 50 50]; %dBHL
Dsd_OHC_Loss = 2/3*Audiogram; % assume 2/3 OHC damage for now
[Cohc_impaired,Cihc_impaired,OHC_Loss]=...
    fitaudiogram(Audiogram_freq,Audiogram,Dsd_OHC_Loss);

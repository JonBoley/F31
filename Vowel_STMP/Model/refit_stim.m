% refit the stimulus for the model
dur_sec=length(signal)/Fs;
dBSPL_before=20*log10(sqrt(mean(signal.^2))/(20e-6));
sfreq=Fs;
sfreqNEW=ANmodel_Fs_Hz;
P=round(sfreqNEW/10);
Q=round(sfreq/10);
if (P/Q*sfreq~=sfreqNEW), disp('Integer sfreq conversion not exact'); end
Nfir=30;
signal_model=resample(signal,P,Q,Nfir);
dBSPL_after=20*log10(sqrt(mean(signal_model.^2))/(20e-6));
if abs(dBSPL_before-dBSPL_after)>2
    error('RESAMPLING CHANGED input by %f dB',dBSPL_before-dBSPL_after);
end

signal_model = signal_model*10^((OALevel_dBSPL-dBSPL_after)/20);

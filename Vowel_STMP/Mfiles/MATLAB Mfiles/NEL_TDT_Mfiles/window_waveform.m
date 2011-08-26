%##########################################################################################
function out_wave = window_waveform(in_wave,ur_Hz,dur_ms)
%##########################################################################################
stimOn_sec = dur_ms / 1000;
rf_time_sec = default_rise_time(stimOn_sec * 1000)/1000;
ramp_nSamples = ceil(ur_Hz * rf_time_sec);
end_sample = size(in_wave, 1);
windowing_vector = zeros(size(in_wave));
windowing_vector(1:ramp_nSamples) = [0:(ramp_nSamples-1)]/ramp_nSamples;
windowing_vector((ramp_nSamples+1):(end_sample-ramp_nSamples)) = 1;
windowing_vector(((end_sample-ramp_nSamples)+1):end) = [(ramp_nSamples-1):-1:0]/ramp_nSamples;
out_wave = windowing_vector .* in_wave;
return;

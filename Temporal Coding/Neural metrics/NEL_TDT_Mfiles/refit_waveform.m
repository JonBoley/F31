%##########################################################################################
function out_wave = refit_waveform(in_wave,ur_Hz,dur_ms)
%##########################################################################################
stimOn_sec = dur_ms / 1000;
required_nSamples = floor(ur_Hz*stimOn_sec);
current_nSamples = size(in_wave, 1);
if (current_nSamples > required_nSamples) % Waveform needs to be truncated.
   out_wave = in_wave(1:required_nSamples);
%    disp(['In function ''refit_waveform'': Input waveform has been truncated ' ...
%       'to fit requested duration.']);
elseif (current_nSamples < required_nSamples) % Waveform needs to be "repeated".
   nRepeats = ceil(required_nSamples / current_nSamples);
   out_wave = repmat(in_wave, nRepeats, 1);
   out_wave = out_wave(1:required_nSamples); % truncate any extra of the final repeat.
%    disp(['In function ''refit_waveform'': Input waveform has been repeated ' ...
%       'to fill requested duration.']);
else
   out_wave = in_wave;
end
return;

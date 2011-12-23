%##########################################################################################
function out_wave = rescale_waveform(in_wave,Level_dBSPL,dBreTONE)
% x.Stimuli.BASELINE.dBreTONE
% x.Stimuli.Condition.Level_dBSPL
%##########################################################################################
global x
% This assumes WAV file input, e.g. peak =1, and uses dBreTONE for this
% vowel to scale to Pascals 
VowelWAV_dBSPL = 20*log10(sqrt(2)/2/20e-6) + dBreTONE;
scale_factor = (10^(Level_dBSPL/20))/(10^(VowelWAV_dBSPL/20));
out_wave = scale_factor * in_wave;
return;
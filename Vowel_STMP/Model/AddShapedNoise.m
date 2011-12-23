function [output,delta_dB]=AddShapedNoise(input,snr)
%
%%% This function takes an input signal 
% and adds shaped noise at specified SNR.
%
%  function [output,delta_dB]=AddShapedNoise(input,snr)
% input - This is is the original signal to create speech shaped noise from
% snr -  Required SNR
% 
% output - The input signal with shaped noise added
% delta_dB - The level difference in dB (positive = higher output)
%
% Ex: [output,delta_dB]=AddShapedNoise(input,-2)
% creates shaped noise with the long term spectrum of noise matching
% that of the input. This noise is added to the wavform at a SNR of -2 dB.
%
% Programmed by Jon Boley on 05/11/10
% (based on SNR_ASA_Compute by Jayaganesh Swaminathan)

dB_before=20*log10(sqrt(mean(input.^2))/(20e-6));

% randomize the phase to create shaped noise
SSN = real(ifft(abs(fft(input))...
    .* exp(1i*2*pi*rand(size(input))))); %steady state noise

% compute the RMS of the original signal
rms_signal=sqrt(mean(input.^2));

% scale the noise based on SNR
rms_scale=10^(snr/20);
rms_noise=sqrt(mean(SSN.^2));
SSN=((1/rms_scale).*rms_signal.*SSN./rms_noise).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output=input+SSN(1:length(input)).';
delta_dB = 20*log10(sqrt(mean(output.^2))/(20e-6))-dB_before;

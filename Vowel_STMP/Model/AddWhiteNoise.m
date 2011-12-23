function [output,delta_dB]=AddWhiteNoise(input,snr)
%
%%% This function takes an input signal
% and adds white noise at specified SNR.
%
%  function [output,delta_dB]=AddWhiteNoise(input,snr)
% input - This is is the original signal to create speech shaped noise from
% snr -  Required SNR
%
% output - The input signal with white noise added
% delta_dB - The level difference in dB (positive = higher output)
%
% Ex: [output,delta_dB]=AddWhiteNoise(input,-2)
% creates white noise to the input. This noise is added to the wavform
% at a SNR of -2 dB.
%
% Programmed by Jon Boley on 05/11/10
% (based on SNR_ASA_Compute by Jayaganesh Swaminathan)

if (snr<Inf)
    dB_before=20*log10(sqrt(mean(input.^2))/(20e-6));
    
    % create white noise
    BBN = rand(size(input));
    
    % compute the RMS of the original signal
    rms_signal=sqrt(mean(input.^2));
    
    % scale the noise based on SNR
    rms_scale=10^(snr/20);
    rms_noise=sqrt(mean(BBN.^2));
    BBN=((1/rms_scale).*rms_signal.*BBN./rms_noise).';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    output=input+BBN(1:length(input)).';
    delta_dB = 20*log10(sqrt(mean(output.^2))/(20e-6))-dB_before;
else
    output=input;
    delta_dB=0;
end

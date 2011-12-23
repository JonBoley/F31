% Set up hearing aid parameters
% NAL equations
% frequency shaping (250,500,1000,2000,3000,4000,6000Hz)
k_NAL = [-18 -8 1 -1 -2 -2 -2]; % dB 
H_3FA = sum(Audiogram(2:4)) / 3; % sum up loss at 500Hz,1kHz,2kHz
X = 0.15*H_3FA;
R = 0.31; % NAL-RP formula  (not quite half-gain rule)
NAL_IG = X + R.*Audiogram + k_NAL; % insertion gain
NAL_IG = max(0,NAL_IG); % no negative gain

% This assumes contant gain over the area of interest
Gain = interp1(Audiogram_freq,NAL_IG,midCF_kHz*1000);

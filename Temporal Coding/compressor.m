function out=compressor(x, threshold, ratio, gain, release, attack, Fs, bitdepth)
% Compressor Function
% Music Signal Processor - EE4896
% Group H
% Jason, Jess, Max, Cody, Jennica
%
% x                       [mono input signal]
% threshold               [dBFS]
% ratio                   [float]
% gain                    [dB]
% release                 [seconds]
% attack                  [seconds]
% Fs                      [Hz]
% bitdepth                [int]

samples = length(x);
T = 1 / Fs;

peak_attack  = 0.005;  % 5ms   from DAFX pp. 100
peak_release = 0.090;  % 130ms  "    "

hysterisis_attack_thresh  = 1 / 2^bitdepth;
hysterisis_release_thresh = 1 / 2^bitdepth;

% Compression Slope
CS = 1 - 1 / ratio;

absx      = abs(x);
y         = zeros(samples+2,1);
g         = ones(samples+2,1);
f         = zeros(samples+2,1);
xpeak     = zeros(samples+2,1);
xpeak_db  = zeros(samples+2,1);
x_peak(1) = abs(x(1));
x_peak(2) = abs(x(2));
in_attack = 0;
gain_lin  = (10^(gain/20));

% disp(sprintf('Thr:%f rat:%f gain:%f rel:%f att:%f FS:%d bit:%d CS:%f',threshold,ratio,gain,release,attack,Fs,bitdepth,CS));

for n=3:(samples)
	% Peak Detection
	xpeak(n)    = max(absx(n)-xpeak(n-1),0)*(1-exp(-2.2*T/peak_attack))+xpeak(n-1)*exp(-2.2*T/peak_release);
	xpeak_db(n) = 20*log10(xpeak(n));

	% Calculate Compression
	overshoot         = abs(min(xpeak_db(n) - threshold, 0));
	compression_in_db = overshoot * CS;
	compression       = 10^(compression_in_db/20);
    f(n) = compression;

    %disp(sprintf('[%d] xpeak:%f (%f dB) over:%f  comp :%f (%f dB)',n,xpeak(n),xpeak_db(n),overshoot,compression,compression_in_db));

	% hysteresis loop
	if (abs(f(n)) - abs(f(n-1))) > hysterisis_attack_thresh
		in_attack = 1;
    elseif (abs(f(n-1)) - abs(f(n))) > hysterisis_release_thresh
		in_attack = 0;
    end;

	% Attack and Release times:
	if in_attack == 1
        % AT
		multiplier = 1 - exp(-2.2*T/attack);
    else
        % RT
		multiplier = 1 - exp(-2.2*T/release);
    end;

    % Gain calculation`
    g(n) = (( f(n) - g(n-1) ) * multiplier ) + g(n-1);
    y(n) = x(n-2) * g(n) * gain_lin;

    %disp(sprintf('[%d] peak: %f (%f dB) mult : %f %-4f -> %-4f  %d%d',n,xpeak(n),xpeak_db(n),multiplier,x(n-2), y(n), in_attack, in_release ));
end;

% Transpose answer to 1-index
out = y(3:samples+2);
function out = TwoBandGain(in,fc,fs,gains)

% design 8th-order Linkwitz-Riley crossover
[b,a] = butter(4,fc/(fs/2),'low');  % low-pass
out_lo=filter(b,a,in);
out_lo=filter(b,a,out_lo);

[b,a] = butter(4,fc/(fs/2),'high');  % high-pass
out_hi=filter(b,a,in);
out_hi=filter(b,a,out_hi);

out = 10^(gains(1)/20)*out_lo + ...
    10^(gains(2)/20)*out_hi;



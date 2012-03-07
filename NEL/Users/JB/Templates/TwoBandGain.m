function out = TwoBandGain(in,fc,fs,gains)

% design 8th-order Linkwitz-Riley crossover
[b,a] = butter(4,fc/(fs/2),'low');  % low-pass
Hd_lp = dfilt.df2(b,a);
Hd_lp = dfilt.cascade(Hd_lp,Hd_lp);

[b,a] = butter(4,fc/(fs/2),'high');  % high-pass
Hd_hp = dfilt.df2(b,a);
Hd_hp = dfilt.cascade(Hd_hp,Hd_hp);

out = 10^(gains(1)/20)*filter(Hd_lp,in) + ...
    10^(gains(2)/20)*filter(Hd_hp,in);
end %function


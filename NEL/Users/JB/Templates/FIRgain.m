function out = FIRgain(in,gain_db,f_hz,Fs)
f = f_hz/(Fs/2);
a = 10.^(gain_db/20);

% make sure we're less than Nyquist
f=f(f<1);
a=a(f<1);

if mod(length(f),2) 
    %create freq pairs
    f=[f(1:end-1) f(end-1)+eps f(end) f(end)+eps 1];
    a=[a(1:end-1) a(end-1) a(end) a(end) 0];
else
    f=[f f(end)+eps 1];
    a=[a a(end) 0];
end

% design 32-order filter
b = firpm(32,f,a);

out=filter(b,1,in);


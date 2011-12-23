% model fiber parameters
CF    = 1e3;  % CF in Hz;   
cohc  = 1.0;  % normal ohc function
cihc  = 1.0;  % normal ihc function
spont = 50; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects
% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds
stimdb = 50; % stimulus intensity in dB SPL
% PSTH parameters
nrep = 50;            % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3; % binwidth in seconds;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

pin = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts; 
pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

[timeout,meout,c1filterout,c2filterout,c1vihc,c2vihc,vihc,synout,psth500k] ...
    = zbcatmodel(pin,CF,nrep,1/Fs,T*1.5,cohc,cihc,spont);
    
psthbins = round(psthbinwidth*Fs);  % number of psth500k bins per psth bin
psthtime = timeout(1:psthbins:end); % time vector for psth
pr = sum(reshape(psth500k,psthbins,length(psth500k)/psthbins))/nrep; % pr of spike in each bin
psth = pr/psthbinwidth; % psth in units of spikes/s
 
figure
subplot(2,1,1)
plot(timeout,[pin zeros(1,length(timeout)-length(pin))])
title('pin')
yl1 = ylim;
subplot(2,1,2)
plot(timeout,meout)
title('meout')
xlabel('Time (s)')
yl2 = ylim;
yl = [min(yl1(1),yl2(1)) max(yl1(2),yl2(2))];
subplot(2,1,1)
ylim(yl)
subplot(2,1,2)
ylim(yl)

figure
subplot(2,1,1)
plot(timeout,c1filterout)
title('c1filterout')
yl1 = ylim;
subplot(2,1,2)
plot(timeout,c2filterout)
title('c2filterout')
xlabel('Time (s)')
yl2 = ylim;
yl = [min(yl1(1),yl2(1)) max(yl1(2),yl2(2))];
subplot(2,1,1)
ylim(yl)
subplot(2,1,2)
ylim(yl)

figure
subplot(3,1,1)
plot(timeout,c1vihc)
title('c1vihc')
yl1 = ylim;
subplot(3,1,2)
plot(timeout,c2vihc)
title('c2vihc')
xlabel('Time (s)')
yl2 = ylim;
subplot(3,1,3)
plot(timeout,vihc)
title('vihc')
yl3 = ylim;
yl = [min([yl1(1) yl2(1) yl3(1)]) max([yl1(2) yl2(2) yl3(2)])];
subplot(3,1,1)
ylim(yl)
subplot(3,1,2)
ylim(yl)
subplot(3,1,3)
ylim(yl)

figure
subplot(2,1,1)
plot(timeout,synout)
title('synout')
xlabel('Time (s)')
subplot(2,1,2)
plot(psthtime,psth)
title('psth')
xlabel('Time (s)')

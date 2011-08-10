clear all;
tdres = 10e-6;

F1=5000;
F2=F1*2+5;

cf = mean([F1 F2]);

index = 1;
ifspike = 0; % No spikes, just sout (= synapse output)
%for spl = -10:10:60,
dur=0.5;
Npts=round(0.5/tdres);
for spl = 30,
  sig(1:Npts) = sin((1:Npts)*tdres*2*pi*F1) + sin((1:Npts)*tdres*2*pi*F2);

  sig(Npts+1:1000) = 0.;
  sig(1:200) = sig(1:200).*(0:199)/199; % add triangular ramps
  sig((Npts-199):Npts) = sig((Npts-199):Npts).*(199:-1:0)/199;
  sig = spl2a(spl)*sig; % scale the sinusoid
  sout(:,index) = an_arlo([tdres,cf,50,1,0,ifspike],sig'); % human
%   sout2(:,index) = an_arlo([tdres,cf,50,1,9,ifspike],sig'); % cat
  rate(index) = mean(sout(:,index)); % find mean rate for rate-level function
%   rate2(index) = mean(sout2(:,index));
  index = index+1;
end;
% figure;
% plot(-10:10:60,rate);
% title('human - Sout (spikes not generated)');
% xlabel('Tone level(dB SPL)')
% ylabel('Rate(sp/sec)')

figure(1);
subplot(211)
plot(1e-3*tdres*(1:length(sout)),sig);
xlabel('Time(msec)')
ylabel('Sig')
subplot(212)
plot(1e-3*tdres*(1:length(sout)),sout);
xlabel('Time(msec)')
ylabel('Rate(sp/sec)')
title('human - Sout (spikes not generated)');

% figure;
% plot(1e-3*tdres*(1:length(sout)),sout2);
% xlabel('Time(msec)')
% ylabel('Rate(sp/sec)')
% title('species 9 - Sout (spikes not generated)');

% figure;
% plot(-10:10:60,rate2);
% title('species 9 - Sout (spikes not generated)');
% xlabel('Tone level(dB SPL)')
% ylabel('Rate(sp/sec)')
% 
% ifspike = 1; %Now, Generate spikes (use last stimulus from above)
% nrep = 300;
% sout3 = an_arlo([tdres,cf,50,1,9,1],sig'); % cat, with spikes
% [sptime,nspikes] = sgmodel([tdres, nrep],sout3);
% mysptime = sptime(1:nspikes);
% 
% if(1)
% len = length(sout3);
% 
% binsize = 0.0001;  % binsize for PST histogram, in sec
% nbins = floor(length(sout3) * tdres/binsize);
% last = 0;
% for i = 1:nbins
%   time = binsize*i;
%   x = find(mysptime<time);
%   y0 = length(x); % total # of spikes in PST up to this bin.
%   y(i) = y0 - last; % subtract out spike in previous bins to get # for this bin.
%   last = y0; 
% end;
% end
% tdres
% nspikes
% 
% figure;
% plot(1e-3* tdres * (1:length(y)),y);
% xlabel('Time(msec)')
% ylabel('Spikes/bin (0.1 ms bin)')



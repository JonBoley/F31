%% PSTH
level = num2str(70);
gain1 = num2str(20);
gain2 = num2str(-20);
CFnum = 20;
load(['archive\phone12\' level 'dBSPL\' 'ihc' '\' gain1 'dBgain_' 'Dec_08_08']);
dur_sec = .400;
SRnum = 3;
binsize = 200e-6;
CFs=(250*2.^(0:log2(4000/250)/29:log2(4000/250)));
cf = num2str(CFs(CFnum));
Period = 1/CFs(CFnum);
PerHist = zeros(1,ceil(Period/binsize)+1);
PSTH = zeros(1,ceil(dur_sec/binsize));
for i=1:length(SpikesA_plus{CFnum,SRnum})
    NumSpikes = length(SpikesA_plus{CFnum,SRnum}{i});
    for j=1:NumSpikes
        PSTH(floor(SpikesA_plus{CFnum,SRnum}{i}(j)/binsize+1)) = ...
            PSTH(floor(SpikesA_plus{CFnum,SRnum}{i}(j)/binsize+1)) + 1;
        PerHist(floor(mod(SpikesA_plus{CFnum,SRnum}{i}(j),Period)/binsize)+1) = ...
            PerHist(floor(mod(SpikesA_plus{CFnum,SRnum}{i}(j),Period)/binsize)+1) + 1;
    end
end

dur_sec = find(PSTH>0,1,'last')*binsize;
PSTH = PSTH/binsize;  % Convert to spikes per second
PerHist = PerHist/max(PerHist);

figure(77), plot(binsize:binsize:dur_sec,PSTH(1:floor(dur_sec/binsize)),'b');
title(['PSTH (CF=' cf 'Hz, ' level 'dBSPL)']); 
xlabel('Time (sec)'); ylabel('Spikes/sec');

figure(78), plot(binsize:binsize:Period,PerHist(1:floor(Period/binsize)),'b');
title(['Period Histogram (CF=' cf 'Hz, ' level 'dBSPL)']); 
xlabel('Time (sec)'); ylabel('Normalized Rate');

% PSTH (+20dB)
% clear all;
dur_sec = .400;
load(['archive\phone12\' level 'dBSPL\' 'ihc' '\' gain2 'dBgain_' 'Dec_08_08']);
PerHist = zeros(1,ceil(Period/binsize)+1);
PSTH = zeros(1,ceil(dur_sec/binsize)+1);
for i=1:length(SpikesA_plus{CFnum,SRnum})
    NumSpikes = length(SpikesA_plus{CFnum,SRnum}{i});
    for j=1:NumSpikes
        PSTH(floor(SpikesA_plus{CFnum,SRnum}{i}(j)/binsize+1)) = ...
            PSTH(floor(SpikesA_plus{CFnum,SRnum}{i}(j)/binsize+1)) + 1;
        PerHist(floor(mod(SpikesA_plus{CFnum,SRnum}{i}(j),Period)/binsize)+1) = ...
            PerHist(floor(mod(SpikesA_plus{CFnum,SRnum}{i}(j),Period)/binsize)+1) + 1;
    end
end

dur_sec = find(PSTH>0,1,'last')*binsize;
PSTH = PSTH/binsize;  % Convert to spikes per second
PerHist = PerHist/max(PerHist);

%
figure(77), hold on; plot(binsize:binsize:dur_sec,PSTH(1:floor(dur_sec/binsize)),'r');
legend([gain1 'dB Gain'],[gain2 'dB Gain']); hold off;

figure(78), hold on; plot(binsize:binsize:Period,PerHist(1:floor(Period/binsize)),'r');
legend([gain1 'dB Gain'],[gain2 'dB Gain']);  hold off;


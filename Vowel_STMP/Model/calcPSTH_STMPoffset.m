%%
% My question:
%     * How far off is PSTH for 0ms neural delay, 
%        as a function of STMP shift

% model 100 fibers, with 0ms neural delay
% calculate correlation delay between CF & STMP versions
% plot offset vs STMPfactor_time

clear all; %close all; home;
% CFfilename = 'STMPvsCF_CF_2012-09-22_171513'; % +/- 0.5oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-22_181815'; % +/- 0.5oct
% CFfilename = 'STMPvsCF_CF_2012-09-22_202543'; % +/- 1.0oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-22_212910'; % +/- 1.0oct
% CFfilename = 'STMPvsCF_CF_2012-09-23_150543'; % 500Hz +/- 1.0oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-23_151324'; % 500Hz +/- 1.0oct
% STMPfilename = 'STMPvsCF_STMP_2012-10-13_205706';
% CFfilename = 'STMPvsCF_CF_2012-10-13_204723';

CFfiles=dir('STMPvsCF_CF_2012-11-24*.mat');
STMPfiles=dir('STMPvsCF_STMP_2012-11-24*.mat');

% fileNum=5;
CFfilename = CFfiles(end).name;
STMPfilename = STMPfiles(end).name;

load(CFfilename,'-regexp','[^CFfilename^STMPfilename]');
PSTH_CF = PSTH;
load(STMPfilename,'-regexp','[^CFfilename^STMPfilename]');

maxlag = round(0.010/PSTHbinWidth_sec); % 10ms (# of bins)
offset = NaN*ones(length(PSTH_CF),1);
for i=1:length(PSTH_CF)
   [c,lag] = xcorr(PSTH_CF{i},PSTH_STMP{i});
   offset_ms(i) = 1e3*PSTHbinWidth_sec*lag(find(c==max(c)));
end

% figure, plot(deltaCF,offset_ms); xlabel('STMP Shift (octaves)'); ylabel('STMP offset (ms)');
% figure, plot(2.^deltaCF,offset_ms); xlabel('STMP freq factor'); ylabel('STMP offset (ms)');
figure, plot(2.^-deltaCF,offset_ms,'.'); xlabel('STMP time factor'); ylabel('STMP offset (ms)'); 
title(sprintf('%1.2fkHz',CF_kHz(1)));

% hmm, that's an interesting figure...
% maybe fit a line to it and see if it makes sense
% if you can't figure out where it comes from, maybe see if you can change it (via stimulus or impairment)



% Greenwood Function for Chinchilla:
A = 163.5;
k = 0.85;
a = 2.1;
x = 0:0.001:1; %proportion of cochlear length
F_Hz = A * (10.^(a*x) - k);
figure, plot(F_Hz,x);

freqs = (CF_kHz(1)*1e3*2.^-deltaCF);
dist = interp1(F_Hz,x,freqs);
figure, scatter(offset_ms,dist);

CF = [1 1.25 1.5 2 2.5 3 4];
m = [3.25 2.95 2.75 2.45 2.2 2.05 1.8]; %offset_ms
figure, plot(log10(CF*1e3/1e3),m,'o');

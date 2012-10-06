%%
% My question:
%     * How far off is PSTH for 0ms neural delay, 
%        as a function of STMP shift

% model 100 fibers, with 0ms neural delay
% calculate correlation delay between CF & STMP versions
% plot offset vs STMPfactor_time

clear all; close all; home;
% CFfilename = 'STMPvsCF_CF_2012-09-22_171513'; % +/- 0.5oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-22_181815'; % +/- 0.5oct
% CFfilename = 'STMPvsCF_CF_2012-09-22_202543'; % +/- 1.0oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-22_212910'; % +/- 1.0oct
% CFfilename = 'STMPvsCF_CF_2012-09-23_150543'; % 500Hz +/- 1.0oct
% STMPfilename = 'STMPvsCF_STMP_2012-09-23_151324'; % 500Hz +/- 1.0oct
STMPfilename = 'STMPvsCF_STMP_2012-10-06_163854';
CFfilename = 'STMPvsCF_CF_2012-10-06_162912';

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
figure, plot(2.^-deltaCF,offset_ms); xlabel('STMP time factor'); ylabel('STMP offset (ms)');

% hmm, that's an interesting figure...
% maybe fit a line to it and see if it makes sense
% if you can't figure out where it comes from, maybe see if you can change it (via stimulus or impairment)


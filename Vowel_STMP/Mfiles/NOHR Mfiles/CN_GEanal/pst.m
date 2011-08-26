function pst(picNum, figNum)
% GE 18Feb2004.  Plots pst of specified picture in current directory.

clear FIG PIC
global FIG PIC



PIC.num = picNum;
PIC.spikeChan = 1; % hard-coded for now.
FIG.binWidth_sec = 5e-4; % hard-coded for now.
FIG.isiMax_sec = 0.05;
FIG.fontSize = 8; % used for axis labelling.

if ~loadPic
   break;
end

if(~exist('figNum', 'var'))
   FIG.num = 100;   % default, if not specified.
else
   FIG.num = figNum;
end

layoutFigure;
do_raster;
do_PST_histo;
do_firstSpikeLatency_histo;
do_CV;
do_ISI_histogram;

%%################################################################################################
function bSuccessful = loadPic()
global PIC

bSuccessful = 1;
picSearchString = sprintf('p%04d*.m', PIC.num);
picMFile = dir(picSearchString);
if (~isempty(picMFile))
   eval( strcat('x = ',picMFile.name(1:length(picMFile.name)-2),';') );
else
   error = sprintf('Picture file p%04d*.m not found.', PIC.num)
   bSuccessful = 0;
   return;
end

PIC.x = x;

%%################################################################################################
function layoutFigure()
global FIG PIC

x = PIC.x;

FIG.handles.main = figure(FIG.num); clf;
set(gcf, 'Name', sprintf ('picture p%04d ''pst'' analysis', PIC.num));
FIG.handles.raster = subplot('Position',[0.1 0.1 0.85 0.1]);
FIG.handles.psth = subplot('Position',[0.1 0.25 0.85 0.15]);
FIG.handles.inst_rate = subplot('Position',[0.1 0.42 0.85 0.1]);
FIG.handles.CV= subplot('Position',[0.1 0.55 0.85 0.15]);
FIG.handles.CV_rect = rectangle;
FIG.handles.isih= subplot('Position',[0.7 0.75 0.25 0.15]);
FIG.handles.isih_rect = rectangle;
FIG.fslh_truncationFactor = 1.1*x.Hardware.Trigger.StmOn/(x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff);
FIG.handles.fslh = subplot('Position',[0.1 0.75 0.85*FIG.fslh_truncationFactor 0.15]);
FIG.handles.fslh_rect = rectangle;

textString = sprintf('picture %03d recorded on %s\n%s', ...
               x.General.picture_number,...
               x.General.date, ...
               x.Stimuli.description{1});
textString = strrep(textString, '_', '\_');
FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 1, textString, 'Units', 'normalized', 'FontSize', FIG.fontSize-1,...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');


%%################################################################################################
function do_raster()
global FIG PIC

x = PIC.x;

subplot(FIG.handles.raster);

plot(x.spikes{1}(:,2),x.spikes{1}(:,1), 'k.', 'MarkerSize', 4);

FIG.xmax = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
xlim([0 FIG.xmax]);
set(gca, 'FontSize', FIG.fontSize);
XLabel('time (sec)');
YLabel('raster');
set(gca, 'YTickLabel', '');
set(gca, 'TickDir', 'out');

%%################################################################################################
function do_PST_histo()
global FIG PIC

x = PIC.x;

subplot(FIG.handles.psth);

lastBin_sec = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
pst_X = [0:FIG.binWidth_sec:lastBin_sec];

plot(pst_X, hist(x.spikes{PIC.spikeChan}(:,2), pst_X), 'k');

xlim([0 FIG.xmax]);
set(gca, 'FontSize', FIG.fontSize);
YLabel('PSTH');
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
set(gca, 'TickDir', 'out');

%%################################################################################################
function do_firstSpikeLatency_histo()
global FIG PIC

x = PIC.x;

subplot(FIG.handles.fslh); hold on;

for (i=1:x.Stimuli.fully_presented_lines)
   lineSpikeTimes = x.spikes{PIC.spikeChan}(find(x.spikes{PIC.spikeChan}(:,1) == i), 2);
   if ~isempty(lineSpikeTimes)
      firstSpikes(i) = lineSpikeTimes(1);
   else
      firstSpikes(i) = NaN;
   end
end
firstSpikes = firstSpikes(find(~isnan(firstSpikes)));
lastBin_sec = max([x.Hardware.Trigger.StmOn/1000 max(firstSpikes)]);
fslh_X = [0:FIG.binWidth_sec:lastBin_sec];

fslh_Y = hist(firstSpikes, fslh_X);
plot(fslh_X, fslh_Y, 'k');

if (max(firstSpikes)*1000 > x.Hardware.Trigger.StmOn)
   FIG.fslh_truncationFactor = 1.1*(max(firstSpikes)*1000)/(x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff);
   FIG.fslh_truncationFactor = min([0.5 FIG.fslh_truncationFactor]);
   set(gca, 'Position', [0.1 0.75 0.85*FIG.fslh_truncationFactor 0.15]);
end



xlim([0 FIG.fslh_truncationFactor*FIG.xmax]);
meanFSL_sec = mean(firstSpikes);
stdFSL_sec = std(firstSpikes);
filt_FSL = trifilt(fslh_Y,3);
peakFSL_sec = fslh_X(find(filt_FSL==max(filt_FSL)));
text_string = sprintf('mean FSL = %.1f +/- %.2f msec   peak FSL @ %.1fmsec',...
   meanFSL_sec*1000, stdFSL_sec*1000, peakFSL_sec*1000);
set(gca, 'FontSize', FIG.fontSize);
YLabel('first-spike');
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
set(gca, 'XTick', get(FIG.handles.raster, 'XTick'));
set(gca, 'TickDir', 'out');

currYLim = YLim;
plot([meanFSL_sec meanFSL_sec], currYLim, 'b');
plot([peakFSL_sec(1) peakFSL_sec(1)], currYLim, 'r');

a = axis;
set(FIG.handles.fslh_rect, 'Position', [meanFSL_sec-(stdFSL_sec/2) a(3) stdFSL_sec a(4)-a(3)]);
shad = 0.9;
set(FIG.handles.fslh_rect, 'EraseMode', 'xor', 'FaceColor', [shad shad shad], 'EdgeColor', [shad shad shad]);
YLim([currYLim(1) 1.33*currYLim(2)]);

text_string_1 = sprintf('FSLH mean = %.1f +/- %.2f msec', meanFSL_sec*1000, stdFSL_sec*1000);
text_string_2 = sprintf('FSLH peak @ %.1f msec', peakFSL_sec(1)*1000);

text('String', text_string_1, 'FontSize', FIG.fontSize-1, 'Units', 'normalized', 'Pos', [0.02 1], ...
      'Color', 'b', 'VerticalAlignment', 'bottom');
text('String', text_string_2, 'FontSize', FIG.fontSize-1, 'Units', 'normalized', 'Pos', [0.02 1], ...
      'Color', 'r', 'VerticalAlignment', 'top');

subplot(FIG.handles.raster);
hold on;
currYLim = YLim;
plot([peakFSL_sec(1) peakFSL_sec(1)], currYLim, 'r');
subplot(FIG.handles.psth);
hold on;
currYLim = YLim;
plot([peakFSL_sec(1) peakFSL_sec(1)], currYLim, 'r');


%%################################################################################################
function do_CV()
global FIG PIC

x = PIC.x;

subplot(FIG.handles.CV); hold on;

lastBin_sec = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
pst_X = [FIG.binWidth_sec/2:FIG.binWidth_sec:lastBin_sec+(FIG.binWidth_sec/2)];
isi_sec = diff(x.spikes{PIC.spikeChan}(:,2));
isi_sec(find(isi_sec<0)) = isi_sec(find(isi_sec<0)) + lastBin_sec;  % "wrap" across lines
orig_isi_sec = isi_sec;
isi_sec(find(isi_sec>FIG.isiMax_sec)) = nan;
isi_spikeTimes_sec = x.spikes{PIC.spikeChan}(find(~isnan(isi_sec)),2);
isi_sec = isi_sec(find(~isnan(isi_sec)));
isi_mean = pst_X; isi_mean(:) = nan;
isi_var = pst_X; isi_var(:) = nan;
inst_rate = pst_X; inst_rate(:) = nan;
for curr_spikeTime_bin = 1:size(pst_X, 2)
   curr_centerTime_sec = pst_X(curr_spikeTime_bin);
   currBin_indices = find(abs(isi_spikeTimes_sec-curr_centerTime_sec)<FIG.binWidth_sec/2);
   if (length(currBin_indices) > 2)
      isi_mean(curr_spikeTime_bin) = mean(isi_sec(currBin_indices));
      isi_var(curr_spikeTime_bin) = std(isi_sec(currBin_indices));
   end
   if (length(currBin_indices) > 0)
      inst_rate(curr_spikeTime_bin) = 1/mean(isi_sec(currBin_indices));
   end
end

set(gca, 'FontSize', FIG.fontSize);
scatter(pst_X, isi_var./isi_mean, 1, 'b');
xlim([0 FIG.xmax]);
ylabel('isi CV');
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');

a = axis;
set(FIG.handles.CV_rect, 'Position', [a(1) 0.35 a(2)-a(1) 0.15]);
shad = 0.6;
set(FIG.handles.CV_rect, 'EraseMode', 'xor', 'FaceColor', [shad 1 shad], 'EdgeColor', [shad 1 shad]);

PIC.isi_sec = orig_isi_sec;

subplot(FIG.handles.inst_rate);
set(gca, 'FontSize', FIG.fontSize);
plot(pst_X, inst_rate, 'k');
% plot(pst_X, ones(1,length(pst_X))./isi_mean, 'k');
xlim([0 FIG.xmax]);
ylabel('inst. rate');
set(gca, 'XTickLabel', '');

%%################################################################################################
function do_ISI_histogram()
global FIG PIC

x = PIC.x;
isi_sec = PIC.isi_sec;

subplot(FIG.handles.isih);
hold on;

isi_X = [0:0.001:max(isi_sec)];
isi_Y = hist(isi_sec, isi_X);
plot(isi_X, isi_Y, 'k');
xlim([-0.02*isi_X(end) isi_X(end)]);
set(gca, 'FontSize', FIG.fontSize-1);


YLabel('isi histo');
isi_X_tickInc = 5*10^(floor(log10(isi_X(end))-1));
while (length([0:isi_X_tickInc:isi_X(end)]) > 5)
   isi_X_tickInc = isi_X_tickInc * 2;
end
set(gca, 'YTickLabel', '', 'XTick', [0:isi_X_tickInc:isi_X(end)]);
set(gca, 'TickDir', 'out');

currYLim = YLim;
a = axis;
set(FIG.handles.isih_rect, 'Position', [0 a(3) FIG.isiMax_sec a(4)-a(3)]);
shad = 0.9;
set(FIG.handles.isih_rect, 'EraseMode', 'xor', 'FaceColor', [shad shad 1], 'EdgeColor', [shad shad 1]);
YLim([currYLim(1) 1.3*currYLim(2)]);

peak_isi_sec = isi_X(find(isi_Y==max(isi_Y)));
text_string = sprintf('ISIH peak @ %.1f msec', peak_isi_sec(1)*1000);
text('String', text_string, 'FontSize', FIG.fontSize-1, 'Units', 'normalized', 'Pos', [0.02 1], ...
      'Color', 'k', 'VerticalAlignment', 'top');

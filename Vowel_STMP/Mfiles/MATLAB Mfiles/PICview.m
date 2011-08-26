function PICview(picNum,excludeLines, figNum)
% File: PICview.m
% M.Heinz 04Jan2005.  Modified from PSTview.m 
%
% Plots raster,and rate vs. rep for a generic picture, just to look at it
%
% Usage: PICview(picNum,excludeLines, figNum)
% picNums: picture number
% excludeLines: vector of any lines to be excluded

temp = excludeLines;
clear FIG PIC excludeLines
global FIG PIC excludeLines
excludeLines = temp;

if ~exist('excludeLines','var')
   excludeLines=[];
end

if(~exist('figNum', 'var'))
   FIG.num = 100;   % default, if not specified.
else
   FIG.num = figNum;
end
FIG.fontSize = 8; % used for axis labelling.

PIC.x=loadPic(picNum);
PIC.num=picNum;

layoutFigure;  % *Keep here
do_raster;     % *Keep here
do_rate;       % *Plot here, calc from Mfile

return;

%%################################################################################################
function layoutFigure()
global FIG PIC

figure_prop_name = {'PaperPositionMode','units','Position'};
% figure_prop_val =  { 'auto'            ,'inches', [8.7083    2.20    5.8333    4]};
if FIG.num==100
	figure_prop_val =  { 'auto'            ,'inches', [0.05    1.0    5.8333    9]};
else
	figure_prop_val =  { 'auto'            ,'inches', [6    1.0    5.8333    9]};
end
FIG.handles.main = figure(FIG.num); clf;
set(gcf,figure_prop_name,figure_prop_val);

NameText=sprintf('picture: %04d; filename: %s', PIC.num,getfileName(PIC.num));
set(gcf, 'Name', NameText);

% Yshift=0.05;
rasterXcorner=0.35; rasterYcorner=0.1; rasterWidth=0.6; rasterHeight=0.8;
rateXcorner=0.07; rateYcorner=rasterYcorner; rateWidth=0.20; rateHeight=rasterHeight;
FIG.handles.rate   = subplot('Position',[rateXcorner rateYcorner rateWidth rateHeight]);
FIG.handles.raster = subplot('Position',[rasterXcorner rasterYcorner rasterWidth rasterHeight]);

titleString = sprintf('picture %s recorded on %s\n%s', ...
               mat2str(PIC.num),...
               PIC.x.General.date, ...
               PIC.x.Stimuli.description{1});
titleString = strrep(titleString, '_', '\_');
FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 1, titleString, 'Units', 'normalized', 'FontSize', FIG.fontSize-1,...
   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

return;

%%################################################################################################
function do_raster()
global FIG PIC excludeLines

subplot(FIG.handles.raster);

SpikeINDs=find(~ismember(PIC.x.spikes{1}(1:end,1)',excludeLines));
plot(PIC.x.spikes{1}(SpikeINDs,2),PIC.x.spikes{1}(SpikeINDs,1), 'k.', 'MarkerSize', 4);

FIG.raster.xmax = (PIC.x.Hardware.Trigger.StmOn + PIC.x.Hardware.Trigger.StmOff) / 1000;
FIG.raster.ymax = ceil(PIC.x.Stimuli.fully_presented_lines/10)*10;
% FIG.raster.ymax = ceil(max(PIC.x.spikes{1}(SpikeINDs,1))/10)*10;
xlim([0 FIG.raster.xmax]);
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
FIG.raster.hx=xlabel('time (sec)');
% set(FIG.raster.hx,'units','norm','pos',[0.4926   -0.01         0])
set(gca,'YTick',0:10:FIG.raster.ymax)
set(gca, 'TickDir', 'out');
return;

%%################################################################################################
function do_rate()
global FIG PIC excludeLines

PIC=calcRatePerLine(PIC);  % Uses driv=[10,410],spont=[600,1000]

% NaN excluded lines
dRateExcludeINDs = intersect(excludeLines,1:length(PIC.RatePerLine.driv));
PIC.RatePerLine.driv(dRateExcludeINDs)=NaN;
sRateExcludeINDs = intersect(excludeLines,1:length(PIC.RatePerLine.spont));
PIC.RatePerLine.spont(sRateExcludeINDs)=NaN;

disp(sprintf('MEAN DRIV RATE = %.f sp/sec',mean(PIC.RatePerLine.driv)))

subplot(FIG.handles.rate);
plot(PIC.RatePerLine.driv,1:length(PIC.RatePerLine.driv),'r','LineWidth',2)
hold on
plot(PIC.RatePerLine.spont,1:length(PIC.RatePerLine.spont),'g','LineWidth',2)
hold off
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse')
set(gca,'YTick',0:10:FIG.raster.ymax)
ylabel('Rep Number');
FIG.rate.hx=xlabel('Rate (sp/s)');
% set(FIG.rate.hx,'units','norm','pos',[-0.2   -0.018         0])

return;



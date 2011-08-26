function SACview(picNums,excludeLines, figNum)
% M.Heinz 29Jul2004.  Modified from pst.m by GE
%
% updated 5/26/05 for SAC analysis ... to be expanded, now just to show PSTview
% updated 4/11/05 to allow TrBFi and new EHrBFi
%
% Plots raster,pst,rate vs. rep, and PerHist for a given condition
% Allows multiple pst-style pics to be concatenated.  Also allows lines to be excluded for each picture.
% Assumes current directory is the data directory
% All calcs are done from external MFile calls, while all plotting is done here
%
% Usage: PSTview(picNums,excludeLines)
% picNums: vector of picture numbers
% excludeLines: cell array with vectors for each picture with any lines to be excluded

% General function to do various calcs
%    *- concatPICS_NOHR.m
%    *- calcRatePerLine.m
%    *- calcPSThist
%    *- calcPERhist
%    *- calcSynchRate_PSThist
%    *- calcSynchRate_PERhist

clear FIG PIC
global FIG PIC

NUMpics=length(picNums);
if ~exist('excludeLines','var')
   excludeLines=cell(1,NUMpics);
end

if(~exist('figNum', 'var'))
   FIG.num = 100;   % default, if not specified.
else
   FIG.num = figNum;
end
FIG.fontSize = 8; % used for axis labelling.
FIG.FeatureColors={'g','r'};  % Troughs,Formants
FIG.FeatureLines={':',':'};  % Troughs,Formants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PIC=concatPICS_NOHR(picNums,excludeLines);

layoutFigure;  % *Keep here
do_raster;     % *Keep here
do_rate;       % *Plot here, calc from Mfile
do_PST_histo;  % *Plot here, calc from Mfile

return;



%%################################################################################################
function layoutFigure()
global FIG PIC

x = PIC.x;

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'inches', [8.7083    2.20    5.8333    8]};
FIG.handles.main = figure(FIG.num); clf;
set(gcf,figure_prop_name,figure_prop_val);

NameText=sprintf('''pst'' analysis for picture(s): p%04d, ', PIC.nums(1));
for PICind=2:length(PIC.nums)
   NameText=[NameText sprintf('p%04d ', PIC.nums(PICind))];
end
NameText=NameText(1:end-1);
set(gcf, 'Name', NameText);

Yshift=0.05;
rasterXcorner=0.35; rasterYcorner=0.65; rasterWidth=0.6; rasterHeight=0.3;
rateXcorner=0.07; rateYcorner=rasterYcorner; rateWidth=0.20; rateHeight=rasterHeight;
psthXcorner=rasterXcorner; psthWidth=rasterWidth; psthHeight=0.15; psthYcorner=rasterYcorner-psthHeight-Yshift; 
perhXcorner=rateXcorner; perhYcorner=psthYcorner; perhWidth=rateWidth; perhHeight=psthHeight;
synchrXcorner=rasterXcorner; synchrWidth=rasterWidth; synchrHeight=psthHeight; synchrYcorner=psthYcorner-synchrHeight-Yshift; ; 
DFTXcorner=rasterXcorner; DFTWidth=rasterWidth; DFTHeight=psthHeight; DFTYcorner=synchrYcorner-DFTHeight-Yshift; 
FIG.handles.rate   = subplot('Position',[rateXcorner rateYcorner rateWidth rateHeight]);
FIG.handles.raster = subplot('Position',[rasterXcorner rasterYcorner rasterWidth rasterHeight]);
FIG.handles.psth   = subplot('Position',[psthXcorner psthYcorner psthWidth psthHeight]);
FIG.handles.perh   = subplot('Position',[perhXcorner perhYcorner perhWidth perhHeight]);
FIG.handles.SynchR = subplot('Position',[synchrXcorner synchrYcorner synchrWidth synchrHeight]);
FIG.handles.SynchR_right = subplot('Position',[synchrXcorner+synchrWidth synchrYcorner .001 synchrHeight]);
FIG.handles.DFT    = subplot('Position',[DFTXcorner DFTYcorner DFTWidth DFTHeight]);
FIG.handles.DFT_right    = subplot('Position',[DFTXcorner+DFTWidth DFTYcorner 0.001 DFTHeight]);

titleString = sprintf('picture(s) %s recorded on %s\n%s', ...
               mat2str(PIC.nums),...
               x.General.date, ...
               x.Stimuli.description{1});
titleString = strrep(titleString, '_', '\_');
FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 1, titleString, 'Units', 'normalized', 'FontSize', FIG.fontSize-1,...
   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

OffsetText={'+','-'};
if strcmp(PIC.TEMPLATE(end-1:end),'BF')
   OffsetIND=strcmp(x.Stimuli.Condition.Offset_Direction,'below ')+1;
end

if strcmp(PIC.TEMPLATE,'SACrlv');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nPolarity Inverted: %s\nNoise Atten= %.f dB\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      PIC.x.Stimuli.Condition.InvertPolarity, ...
      PIC.x.Stimuli.attens(PIC.CONDind));
end
text(0.03, synchrYcorner+synchrHeight-Yshift/3, ConditionString, 'Units', 'normalized', 'FontSize', FIG.fontSize,...
   'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','b');

return;

%%################################################################################################
function do_raster()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.raster);

plot(x.spikes{1}(:,2),x.spikes{1}(:,1), 'k.', 'MarkerSize', 4);

FIG.raster.xmax = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
FIG.raster.ymax = ceil(x.Stimuli.fully_presented_lines/10)*10;
xlim([0 FIG.raster.xmax]);
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
FIG.raster.hx=xlabel('time (sec)');
set(FIG.raster.hx,'units','norm','pos',[0.4926   -0.01         0])
set(gca,'YTick',0:10:FIG.raster.ymax)
set(gca, 'TickDir', 'out');
return;

%%################################################################################################
function do_rate()
global FIG PIC

PIC=calcRatePerLine(PIC);  % Uses driv=[10,410],spont=[600,1000]

subplot(FIG.handles.rate);
plot(PIC.RatePerLine.driv,1:length(PIC.RatePerLine.driv),'r','LineWidth',2)
hold on
plot(PIC.RatePerLine.spont,1:length(PIC.RatePerLine.spont),'g','LineWidth',2)
hold off
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse')
set(gca,'YTick',0:10:FIG.raster.ymax)
YLabel('Rep Number');
FIG.rate.hx=xlabel('(sp/s)');
set(FIG.rate.hx,'units','norm','pos',[-0.2   -0.018         0])

return;


%%################################################################################################
function do_PST_histo()
global FIG PIC FeaturesText

% PIC=calcSynchRate_PST(PIC);  % Calculates PST as well
PIC=calcPST(PIC);

%%%% Plot PST Histogram
subplot(FIG.handles.psth);
plot(PIC.PST.pst_X_sec, PIC.PST.pst_Y_sps, 'k');
hold on
% plot(PIC.PST.pst_X_sec(PIC.PST.drivenPSTinds), PIC.PST.pst_Y_sps(PIC.PST.drivenPSTinds), 'b');
xlim([0 FIG.raster.xmax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'TickDir', 'out');
FIG.psth.hy=ylabel('(sp/s)');
set(FIG.psth.hy,'units','norm','pos',[-0.108    0.4874         0])
FIG.psth.hx=xlabel('time (sec)');
set(FIG.psth.hx,'units','norm','pos',[0.4926   -0.01         0])


return;





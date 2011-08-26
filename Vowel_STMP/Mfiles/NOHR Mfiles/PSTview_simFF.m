function PSTview_simFF(picNums,excludeLines,figNum)
% M.Heinz 28Sep2004.  Modified from PSTview.m 
% Plots original PSTview of 'reBF' data, and plots simulated 'reFF' data.
% - reBF data is a stimulus that has been shifted above (or below) BF.
% - simulated FF data is translated backwards, such that it models a neuron with BF shifted below (or above)
%   the nominal BF, responding to a stimulus at the nominal BF.
%
% ORIGINAL PSTview.m
% Plots raster,pst,rate vs. rep, and PerHist for a given condition
% Allows multiple pst-style pics to be concatenated.  Also allows lines to be excluded for each picture.
% Assumes current directory is the data directory
% All calcs are done from external MFile calls, while all plotting is done here
%
% Usage: PSTview(picNums,excludeLines)
% picNums: vector of picture numbers
% excludeLines: cell array with vectors for each picture with any lines to be excluded


NUMpics=length(picNums);
if ~exist('excludeLines','var')
   excludeLines=cell(1,NUMpics);
end

if(~exist('figNum', 'var'))
   figNum = 100;   % default, if not specified.
end

%% Run PSTview on original data
PSTview(picNums,excludeLines,figNum-1)
clear FIG PIC  % PSTview uses same
global FIG PIC

FIG.num=figNum;
FIG.fontSize = 8; % used for axis labelling.
FIG.FeatureColors={'g','r'};  % Troughs,Formants
FIG.FeatureLines={':',':'};  % Troughs,Formants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PIC=concatPICS_NOHR(picNums,excludeLines);

% Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
PIC=simFF_PICshift(PIC);


layoutFigure;  % *Keep here
do_raster;     % *Keep here
do_rate;       % *Plot here, calc from Mfile
do_PST_histo;  % *Plot here, calc from Mfile
do_PerHist(excludeLines);    % *Plot here, calc from Mfile

return;



%%################################################################################################
function layoutFigure()
global FIG PIC

x = PIC.x;

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'inches', [2.8    2.20    5.8333    8]};
FIG.handles.main = figure(FIG.num); clf;
set(gcf,figure_prop_name,figure_prop_val);

NameText=sprintf('''pst_simFF'' analysis for picture(s): p%04d, ', PIC.nums(1));
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

if strcmp(PIC.TEMPLATE,'EHrBF');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %s%.2f octaves\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      OffsetText{OffsetIND}, ...
      x.Stimuli.Condition.FreqOffset_octs, ...
      x.Stimuli.Condition.Feature, ...
      x.Stimuli.Condition.FormsAtHarmonics, ...
      x.Stimuli.Condition.InvertPolarity, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'TrBF')
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %s%.2f octaves\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      OffsetText{OffsetIND}, ...
      x.Stimuli.Condition.FreqOffset_octs, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'TrBFi')
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %.2f octaves\nFeature: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      PIC.x.Stimuli.Used.OctShifts_List(PIC.CONDind), ...
      PIC.x.Stimuli.Used.Features_List{PIC.CONDind}, ...
      PIC.x.Stimuli.Used.Levels_dBSPL_List(PIC.CONDind));
elseif strcmp(PIC.TEMPLATE,'TrFF')
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nFF = %.3f kHz\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      x.Stimuli.Condition.BaseFrequency_kHz, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'EHrFF');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nFF = %.3f kHz\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      x.Stimuli.Condition.BaseFrequency_kHz, ...
      x.Stimuli.Condition.Feature, ...
      x.Stimuli.Condition.FormsAtHarmonics, ...
      x.Stimuli.Condition.InvertPolarity, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'EHrBFi');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %.2f octaves\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      PIC.x.Stimuli.Used.OctShifts_List(PIC.CONDind), ...
      PIC.x.Stimuli.Used.Features_List{PIC.CONDind}, ...
      PIC.x.Stimuli.Condition.FormsAtHarmonics, ...
      PIC.x.Stimuli.Condition.InvertPolarity, ...
      PIC.x.Stimuli.Used.Levels_dBSPL_List(PIC.CONDind));
elseif strcmp(PIC.TEMPLATE,'EHvNrBFi');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %.2f octaves\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\nNoise Atten= %.f dB\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      PIC.x.Stimuli.Used.OctShifts_List(PIC.CONDind), ...
      PIC.x.Stimuli.Used.Features_List{PIC.CONDind}, ...
      PIC.x.Stimuli.Condition.FormsAtHarmonics, ...
      PIC.x.Stimuli.Condition.InvertPolarity, ...
      PIC.x.Stimuli.Used.Levels_dBSPL_List(PIC.CONDind), ...
      PIC.x.Stimuli.Used.NoiseAttens_dB_List(PIC.CONDind));
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

PIC=calcSynchRate_PST(PIC);  % Calculates PST as well

%%%% Plot PST Histogram
subplot(FIG.handles.psth);
plot(PIC.PST.pst_X_sec, PIC.PST.pst_Y_sps, 'k');
hold on
plot(PIC.PST.pst_X_sec(PIC.PST.drivenPSTinds), PIC.PST.pst_Y_sps(PIC.PST.drivenPSTinds), 'b');
xlim([0 FIG.raster.xmax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'TickDir', 'out');
FIG.psth.hy=ylabel('(sp/s)');
set(FIG.psth.hy,'units','norm','pos',[-0.108    0.4874         0])
FIG.psth.hx=xlabel('time (sec)');
set(FIG.psth.hx,'units','norm','pos',[0.4926   -0.01         0])

%%%% Plot Synchronized Rate from PST Histogram
subplot(FIG.handles.SynchR);
plot(PIC.SynchRate_PST.FFTfreqs,abs(PIC.SynchRate_PST.SynchRate_PST))
hold on
% ADD Synch Coef Labels to right Axis
set(FIG.handles.SynchR_right,'YAxisLocation','right','YTick',[0:.2:1],'FontSize',FIG.fontSize,'TickDir','out','TickLength',[.03 0.025])

% Plot BF 
plot(PIC.BF_Hz*ones(1,2),[0 1e4],'k--')
FIG.SynchR.ymax=max(abs(PIC.SynchRate_PST.SynchRate_PST));
text(PIC.BF_Hz,FIG.SynchR.ymax*1.06,'BF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')

% Plot FF (Fixed Frequency, if used)
if strcmp(PIC.TEMPLATE(end-1:end),'FF')
   plot(PIC.FixedFreq_Hz*ones(1,2),[0 1e4],'k-.')
   text(PIC.FixedFreq_Hz,FIG.SynchR.ymax*1.06,'FF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
end

% Plot Features (e.g., Tone Freq, Formants and Troughs)
for i=1:length(PIC.FeatureFreqs_Hz)
   plot(PIC.FeatureFreqs_Hz(i)*ones(1,2),[0 1e4],':','Color',FIG.FeatureColors{mod(i,2)+1})
   text(PIC.FeatureFreqs_Hz(i),FIG.SynchR.ymax*.98,FeaturesText{i},'FontSize',FIG.fontSize, ...
      'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center','Color',FIG.FeatureColors{mod(i,2)+1})
end

FIG.SynchR.xmax=min([ceil(max([PIC.FeatureFreqs_Hz PIC.BF_Hz])/5000)*5000 20000]);

xlim([0 FIG.SynchR.xmax])
ylim([0 FIG.SynchR.ymax])
set(gca, 'FontSize', FIG.fontSize);
% xlabel('Frequency (Hz)');
ylabel('Synchronized Rate (sp/sec)')
hold off

return;


%%################################################################################################
function do_PerHist(excludeLines)
global FIG PIC FeaturesText

PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well

%%%% Plot PERIOD histogram
subplot(FIG.handles.perh);
bar(PIC.PERhist.PERhist_X_sec*1000,PIC.PERhist.PERhist,1)
xlim([0 1/PIC.FundamentalFreq_Hz*1000])
set(gca, 'FontSize', FIG.fontSize);
text(.05,.92,sprintf('F0 = %.f (Hz)',PIC.FundamentalFreq_Hz),'units','norm', 'FontSize', FIG.fontSize)
hy=ylabel('(sp/s)');
set(hy,'units','norm')
set(hy,'pos',[-0.27 0.4914 0])
set(gca,'XTick',[0 .5 1]*1/PIC.FundamentalFreq_Hz*1000)
set(gca, 'TickDir', 'out');
xlabel('time (msec)')

%%%% Plot Synchronized Rate from PERIOD histogram
subplot(FIG.handles.DFT);
%%%%%%% Get both Synhchronized Rate and SynchCoef labels
plot(PIC.SynchRate_PERhist.FFTfreqs,abs(PIC.SynchRate_PERhist.SynchRate_PERhist),'b-x')
hold on
% Plot BF (Target Frequency)
FIG.DFT.ymax=max(abs(PIC.SynchRate_PERhist.SynchRate_PERhist));
text(PIC.BF_Hz,FIG.DFT.ymax*1.06,'BF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
plot(PIC.BF_Hz*ones(1,2),[0 1e4],'k--')

% Plot FF (Fixed Frequency, if used)
if strcmp(PIC.TEMPLATE(end-1:end),'FF')
   plot(PIC.FixedFreq_Hz*ones(1,2),[0 1e4],'k-.')
   text(PIC.FixedFreq_Hz,FIG.DFT.ymax*1.06,'FF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
end

% Plot Features (e.g., Tone Freq, Formants and Troughs)
for i=1:length(PIC.FeatureFreqs_Hz)
   plot(PIC.FeatureFreqs_Hz(i)*ones(1,2),[0 1e4],':','Color',FIG.FeatureColors{mod(i,2)+1})
   text(PIC.FeatureFreqs_Hz(i),FIG.SynchR.ymax*.98,FeaturesText{i},'FontSize',FIG.fontSize, ...
      'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center','Color',FIG.FeatureColors{mod(i,2)+1})
end

xlim([0 FIG.SynchR.xmax])
ylim([0 FIG.DFT.ymax])
set(gca, 'FontSize', FIG.fontSize);
xlabel('Frequency (Hz)');
ylabel('Synchronized Rate (sp/sec)')

% ADD Synch Coef Labels to right Axis
set(FIG.handles.DFT_right,'YAxisLocation','right','YTick',[0:.2:1],'FontSize',FIG.fontSize,'TickDir','out','TickLength',[.03 0.025])
subplot(FIG.handles.DFT_right)
ht=title('Synch','FontSize',FIG.fontSize);
set(ht,'units','norm','position',[7    1.0427 0])

SigText={' ','*'};
%%% List: Numspikes,PIC.FundamentalFreq_Hz,Synch,Phase,Signif
if sum(strcmp(PIC.TEMPLATE,{'EHrBF','EHrFF'}))
   FeatIND=find(strcmp(FeaturesText,PIC.x.Stimuli.Condition.Feature(1:2)));
   SynchText=sprintf('Nsps=%d, %s =%.2f Hz, Synch=%.2f[%s], Ph= %.2f rad',PIC.PERhist.NumDrivenSpikes, ...
      FeaturesText{FeatIND},PIC.FeatureFreqs_Hz(FeatIND),PIC.SynchRate_PERhist.FeatureSynchs(FeatIND), ...
      SigText{PIC.SynchRate_PERhist.FeatureRaySig(FeatIND)+1},PIC.SynchRate_PERhist.FeaturePhases(FeatIND));
elseif sum(strcmp(PIC.TEMPLATE,{'TrBF','TrFF','TrBFi'}))
   SynchText=sprintf('Nsps=%d, Tone =%.2f Hz, Synch=%.2f[%s], Ph= %.2f rad',PIC.PERhist.NumDrivenSpikes, ...
      PIC.FundamentalFreq_Hz,PIC.SynchRate_PERhist.FeatureSynchs(1), ...
      SigText{PIC.SynchRate_PERhist.FeatureRaySig(1)+1},PIC.SynchRate_PERhist.FeaturePhases(1));
   FeatIND=1;
elseif sum(strcmp(PIC.TEMPLATE,{'EHrBFi','EHvNrBFi'}))
   FeatIND=find(strcmp(FeaturesText,deblank(PIC.x.Stimuli.Used.Features_List{PIC.CONDind})));
   SynchText=sprintf('Nsps=%d, %s =%.2f Hz, Synch=%.2f[%s], Ph= %.2f rad',PIC.PERhist.NumDrivenSpikes, ...
      FeaturesText{FeatIND},PIC.FeatureFreqs_Hz(FeatIND),PIC.SynchRate_PERhist.FeatureSynchs(FeatIND), ...
      SigText{PIC.SynchRate_PERhist.FeatureRaySig(FeatIND)+1},PIC.SynchRate_PERhist.FeaturePhases(FeatIND));
end
subplot(FIG.handles.psth);
FIG.psth.htext=text(0.99,0.96,SynchText,'FontSize',FIG.fontSize,'HorizontalAlign','right', ...
   'Units','norm','VerticalAlign','top');

subplot(FIG.handles.figBox)
clear FullSynchText
FullSynchText{1}=sprintf('Feature   Synch   Phase (rad)');
for i=1:length(PIC.FeatureFreqs_Hz)
   FullSynchText{i+1} = sprintf('    %s%s        %.2f%s       %.2f',FeaturesText{i}, ...
      SigText{(i==FeatIND)+1},PIC.SynchRate_PERhist.FeatureSynchs(i), ...
      SigText{PIC.SynchRate_PERhist.FeatureRaySig(i)+1},PIC.SynchRate_PERhist.FeaturePhases(i));
end
text(0.01, .03,strvcat(FullSynchText), 'Units', 'normalized', 'FontSize', FIG.fontSize,...
   'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color','k');

return;



function rc = url(trackNum, unitNum, data_folder)

% Written by GE 19Apr2004.

global uinfo FIG;

if ~exist('data_folder')
   uinfo.dataFolder = [cd '\'];
else
   uinfo.dataFolder = data_folder;
end
uinfo.TNum = trackNum;
uinfo.UNum = unitNum;
uinfo.unitName = sprintf('%d.%02d', trackNum, unitNum);
uinfo.trigFlag = 0;

uinfo.stimTypes = {'T2' 'T1' 'TB' 'FN' 'SP' 'TX' 'B2' 'B1'};
offset = 0;
uinfo.stimPanels = [1+offset:2:7+offset 2+offset:2:8+offset];
uinfo.stimAttenField = {'Tone' 'Tone' 'Tone' 'list' 'list' 'Tone' 'Tone' 'Tone' };


error = get_unit_info;
if (error)
   break;
end
layout_figure;

FIG.comments.count = 0;
FIG.comments.text = '';
FIG.triggerWarns.count = 0;
FIG.triggerWarns.text = sprintf('Trigger warnings:\n\n');

for i = 1:length(uinfo.stimTypes)
   plot_RLs(i);
end

if (FIG.comments.count > 0)
   warndlg(FIG.comments.text);
end

if (FIG.triggerWarns.count > 0)
   warndlg(FIG.triggerWarns.text);
end


%%##################################################################################################
function error = get_unit_info()
global uinfo Data;

error = 0;

% Temporarily switch to experiment data directory and load "info file":
origDir = cd;
eval (['cd ''' uinfo.dataFolder '''']);
load DataInfoFile;
eval (['cd ''' origDir '''']);

% Setup appropriate parameters based on loaded unit info data:

     % General info:
if (uinfo.TNum > size(Data.Info, 1) |  uinfo.UNum > size(Data.Info, 2))
	error_report = 'Invalid track or unit number.'
	error = -1;
	return
end
infoField = Data.Info{uinfo.TNum, uinfo.UNum};
if (isempty(infoField))
   error_report = 'Invalid track or unit number.'
   error = -1;
   return
end
uinfo.text = sprintf('Experiment %s     Unit %s     %s %s unit     BF = %.1fkHz     threshold = -%ddB ', ...
               strrep(Data.General.date,'_', '\_'), ...
               uinfo.unitName, ...
               infoField.Location, ...
               infoField.UnitType, ...
               infoField.BF_kHz, ...
               infoField.Threshold_dBatten ...
             );
uinfo.BF_kHz = infoField.BF_kHz;
% uinfo.spikeChan = infoField.spikeChan;  % ge debug 02Apr2004, enable line when available.
% uinfo.ignoreComments = infoField.ignoreComments;  % ge debug 02Apr2004, enable line when available.
% uinfo.ignoreTriggers = infoField.ignoreTriggers;  % ge debug 02Apr2004, enable line when available.
uinfo.spikeChan = 1;
uinfo.ignoreComments = 0;
uinfo.ignoreTriggers = 0;  

for i = 1:length(uinfo.stimTypes)
   eval(sprintf('uinfo.pics.%s = get_pics(''%s'');', uinfo.stimTypes{i}, uinfo.stimTypes{i}));
end

%%##################################################################################################
function layout_figure()
global uinfo FIG;

FIG.fontSize = 8;

FIG.handles.main = figure(5); clf; hold on;  % ge debug 19Mar2004, figure # is hard-coded, for now
set(gcf, 'Color', 'w');
set(gcf, 'Name', 'general rate-levels (url.m)');


FIG.M = 4;
FIG.N = 2;
for i = 1:length(uinfo.stimTypes)
   eval(sprintf('[FIG.handles.subplot.%s] = setupSubplot(''%s'', %d);', ...
      uinfo.stimTypes{i}, uinfo.stimTypes{i}, uinfo.stimPanels(i)));
end

FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 0.99, uinfo.text, ...
   'Units', 'normalized', 'FontName', 'Times', 'FontSize', FIG.fontSize,...
   'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');


%%##################################################################################################
function plot_RLs(stimIndex)
global uinfo FIG Data;

stimType = uinfo.stimTypes{stimIndex};

eval(['subplot(FIG.handles.subplot.' stimType '.panel);']);
pics = eval(['uinfo.pics.' stimType ';']);
textHandle = eval(['FIG.handles.subplot.' stimType '.text;']);
rectHandle = eval(['FIG.handles.subplot.' stimType '.rect;']);
set(gca, 'FontSize', FIG.fontSize-1);

uinfo.allSpont = [];
uinfo.allSpont_att = [];

legendLabels = cell(1, length(pics));
picCounter = 0;
for currPic = pics
   [error x] = loadPic(currPic);
   if (~error)
      picCounter = picCounter + 1;
      checkPic(x);
      set(textHandle, 'Visible', 'Off')
      [driv spont] = compute_driv_spont_rates(x);
      attenField = eval(['x.Line.attens.' uinfo.stimAttenField{stimIndex}]);
      if (~isnan(attenField(1,1)))
         att_dB = attenField(:,1);
      elseif (~isnan(attenField(1,2)))
         att_dB = attenField(:,2);
      end
      h_driv(picCounter) = plot(att_dB(1:size(driv,2)), trifilt(driv, 3), 'LineWidth', 2);
      colorOrd = get(gca,'ColorOrder');
      colorOrd = [colorOrd(2,:);colorOrd(3:end,:);colorOrd(1,:)];
      set(gca, 'ColorOrder', colorOrd);
      if ~strcmp(stimType, 'TX')
         legendLabels{picCounter} = sprintf('%s (p%03d)', stimType, currPic);
      else
         legendLabels{picCounter} = sprintf('T%.1f (p%03d)', x.Stimuli.main.tone.freq/1000, currPic);
      end
      plot(att_dB(1:size(spont,2)), trifilt(spont, 3), ':k');
      uinfo.allSpont = [uinfo.allSpont spont];
      uinfo.allSpont_att = [uinfo.allSpont_att att_dB(1:size(spont,2))'];
   end
end

if (picCounter > 0)
	legend(h_driv, legendLabels, 2);
	legend('boxoff');
end

if ~isempty(uinfo.allSpont)
   meanSpont = mean(uinfo.allSpont);
   stdSpont = std(uinfo.allSpont);
   set(textHandle, 'Visible', 'On', 'Position', [1.0 1.0], 'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'top', 'FontSize', FIG.fontSize-1);
   set(textHandle, 'String', sprintf('SR = %.1f +/- %.1f spk/s', meanSpont, std(uinfo.allSpont)));
   xMin = min(uinfo.allSpont_att)+1;
   xMax = max(uinfo.allSpont_att)-1;
   plot([xMin xMax], [meanSpont meanSpont], 'r');
	if (xMax-xMin > 0) & (stdSpont>0)
       set(rectHandle, 'Position', [xMin meanSpont-(stdSpont/2) xMax-xMin stdSpont], 'Visible', 'on');
       shad = 0.85;
       set(rectHandle, 'EraseMode', 'xor', 'FaceColor', [shad shad shad], 'EdgeColor', [shad shad shad]);
   end
end


%%#######################################################################
function [error, x] = loadPic(picNum)     % Load picture
global uinfo;

if (~isempty(dir(sprintf('%sp%04d*.m', uinfo.dataFolder, picNum))))
   error = 0;
   origDir = cd;
   eval (['cd ''' uinfo.dataFolder '''']);
   picMFile = dir(sprintf('p%04d*.m', picNum));
   eval( strcat('x = ', picMFile.name(1:length(picMFile.name)-2),';') );
   eval (['cd ''' origDir '''']);
else
   error_report = sprintf('Picture file p%04d*.m not found.', picNum)
   error = -1;
   x = [];
   return;
end


%%#######################################################################
function [driv, spont] = compute_driv_spont_rates(x)
global uinfo;
     % compute driven and spontaneous rates for each line
spikeTimes = x.spikes{uinfo.spikeChan};
driven_dur_sec = x.Hardware.Trigger.StmOn / 1000;
spont_dur_sec = x.Hardware.Trigger.StmOff / 1000;
for line = 1:x.Stimuli.fully_presented_lines
   spikeIndices = find( (spikeTimes(:,1) == line ) );
   driv(line) = length (find(spikeTimes(spikeIndices,2) <= driven_dur_sec));
   spont(line) = length(spikeIndices) - driv(line);
end
driv = driv / driven_dur_sec;
spont = spont / spont_dur_sec;


%%#######################################################################
function checkPic(x)
global uinfo FIG;

if ~(uinfo.ignoreTriggers)
   if strcmp(x.General.trigger, 'Fair ') | strcmp(x.General.trigger, 'Poor ')
      FIG.triggerWarns.count = FIG.triggerWarns.count + 1;
      picNum = x.General.picture_number;
      shortDesc = x.Stimuli.short_description;
      triggerText = x.General.trigger;
      picText = sprintf('      p%03d ('' %s'') -- trigger = ''%s''\n', picNum, shortDesc, triggerText);
      FIG.triggerWarns.text = [FIG.triggerWarns.text picText];
   end
end

if ~(uinfo.ignoreComments) & ~isempty(x.General.comment)
	FIG.comments.count = FIG.comments.count + 1;
   picNum = x.General.picture_number;
   shortDesc = x.Stimuli.short_description;
   comment = x.General.comment;
	commentText = sprintf('Comment for picture p%03d ( %s):\n      ''''%s''''\n\n', picNum, shortDesc, comment);
	FIG.comments.text = [FIG.comments.text commentText];
end

%%#########################################################################
function pics = get_pics(stimType)
global uinfo Data;

pics = [];
XXfield = eval(sprintf('Data.RLFs{uinfo.TNum,uinfo.UNum}.%s;', stimType));
for i = 1:size(XXfield, 2)
   if (isempty(find(Data.IgnorePicNums == XXfield{i}.picNUM)));
      pics = [pics XXfield{i}.picNUM];
   end
end

%%########################################################################
function [subplotHandles] = setupSubplot(stimType, panelNumber)
global FIG

subplotHandles.panel = subplot(FIG.M,FIG.N,panelNumber);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse');
YLabel('rate (sp/s)');
if (panelNumber==FIG.M*FIG.N) | (panelNumber==FIG.M*FIG.N-1)
   XLabel('atten (dB)');
end
subplotHandles.text = text(0.5, 0.5, sprintf('No valid %s pics\n available for this unit.', stimType),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);
subplotHandles.rect = rectangle('Position', [0 0 100 100], 'Visible', 'off');

function rc = uchar(trackNum, unitNum, data_folder)

% Written by GE 19Mar2004.  Generates a layout of characterizations for a given unit, including
%    basic info, post-stimulus time histograms, response maps, rate-level functions.

global uinfo FIG;

if ~exist('data_folder')
   uinfo.dataFolder = [cd '\'];
else
   uinfo.dataFolder = data_folder;
end
uinfo.TNum = trackNum;
uinfo.UNum = unitNum;
uinfo.unitName = sprintf('%d.%02d', trackNum, unitNum);

error = get_unit_info;
if (error)
   break;
end
layout_figure;

FIG.comments.count = 0;
FIG.comments.text = '';
FIG.triggerWarns.count = 0;
FIG.triggerWarns.text = sprintf('Trigger warnings:\n\n');

%% ge debug 23Mar2004: incorporate appropriate CALIB file for each plot (esp. RMs?)
plot_RLs;
plot_PSTs;
plot_RMs;

if (FIG.comments.count > 0)
   warndlg(FIG.comments.text);
end

if (FIG.triggerWarns.count > 0)
   warndlg(FIG.triggerWarns.text);
end


%%##################################################################################################
function error = get_unit_info()
global uinfo;

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
uinfo.text = sprintf('UNIT %s (experiment date %s)  \n%s %s unit: BF = %.1fkHz, threshold = -%ddB ', ...
               uinfo.unitName, ...
               strrep(Data.General.date,'_', '\_'), ...
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

     % BF rate-level pics:
uinfo.pics.BFRL = [];
RLVfield = Data.characterization{uinfo.TNum,uinfo.UNum}.RLV;
for i = 1:size(RLVfield, 2)
   if (isempty(find(Data.IgnorePicNums == RLVfield{i}.picNUM)));
      uinfo.pics.BFRL = [uinfo.pics.BFRL RLVfield{i}.picNUM];
   end
end

     % noise rate-level pics:
uinfo.pics.NORL = [];
NOfield = Data.characterization{uinfo.TNum,uinfo.UNum}.NO;
for i = 1:size(NOfield, 2)
	if (isempty(find(Data.IgnorePicNums == NOfield{i}.picNUM)));
		uinfo.pics.NORL = [uinfo.pics.NORL NOfield{i}.picNUM];
   end
end

     % PST pics:
uinfo.pics.PST = [];
PSTfield = Data.characterization{uinfo.TNum,uinfo.UNum}.PST;
for i = 1:size(PSTfield, 2)
   if (isempty(find(Data.IgnorePicNums == PSTfield{i}.picNUM)));
      uinfo.pics.PST = [uinfo.pics.PST PSTfield{i}.picNUM];
   end
end
PSTafield = Data.characterization{uinfo.TNum,uinfo.UNum}.PSTa;
for i = 1:size(PSTafield, 2)
   if (isempty(find(Data.IgnorePicNums == PSTafield{i}.picNUM)));
      uinfo.pics.PST = [uinfo.pics.PST -1*PSTafield{i}.picNUM]; % "flag" async pics with a negative picture #
   end
end


     % response map pics:
uinfo.pics.RM = [];
RMfield = Data.characterization{uinfo.TNum,uinfo.UNum}.RM;
for i = 1:size(RMfield, 2)
   if (isempty(find(Data.IgnorePicNums == RMfield{i}.picNUM)));
      uinfo.pics.RM = [uinfo.pics.RM RMfield{i}.picNUM];
   end
end
RM2Tfield = Data.characterization{uinfo.TNum,uinfo.UNum}.RM_2T;
for i = 1:size(RM2Tfield, 2)
   if (isempty(find(Data.IgnorePicNums == RM2Tfield{i}.picNUM)));
      uinfo.pics.RM = [uinfo.pics.RM -1*RM2Tfield{i}.picNUM]; % "flag" 2T RM pics with a negative picture #
   end
end

%%##################################################################################################
function layout_figure()
global uinfo FIG;

FIG.fontSize = 8;

FIG.handles.main = figure(1); clf; hold on;  % ge debug 19Mar2004, figure # is hard-coded, for now
set(gcf, 'Color', 'w');

set(gcf, 'Name', 'general unit characterization (uchar.m)');

FIG.handles.RL = subplot('Position',[0.1 0.7 0.35 0.18]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse');
YLabel('rate (sp/s)');
XLabel('atten (dB)');
FIG.RLtext = text(0.5, 0.5, sprintf('No valid RL pics\n available for this unit.'),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);
FIG.handles.RLrect = rectangle('Position', [0 0 100 100]);


pstLeft = 0.1; pstBottom = 0.1; pstWidth = 0.34; pstHeight = 0.5;
FIG.handles.PSTright = subplot('Position',[pstLeft+pstWidth pstBottom 0.01 pstHeight]); % for labelling right axis
set(gca, 'FontSize', FIG.fontSize-1);
set(gca, 'XTickMode', 'manual');
set(gca, 'YAxisLocation', 'right');
set(gca, 'YTickMode', 'manual');
FIG.handles.PST = subplot('Position',[pstLeft pstBottom pstWidth pstHeight]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'YAxisLocation', 'left');
set(gca, 'XTickLabelMode', 'auto');
YLabel('db atten for PSTH');
XLabel('time (sec)');
set(gca, 'YTickLabel', '');
FIG.PSTtext = text(0.5, 0.5, sprintf('No valid PST pics\n available for this unit.'),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);
subplot(FIG.handles.PSTright);
    

rmLeft = 0.6; rmBottom = 0.1; rmWidth = 0.29; rmHeight = 0.78;
FIG.handles.RMright = subplot('Position',[rmLeft+rmWidth rmBottom 0.01 rmHeight]); % for labelling right axis
set(gca, 'FontSize', FIG.fontSize-1);
set(gca, 'XTickMode', 'manual');
set(gca, 'YAxisLocation', 'right');
set(gca, 'YTickMode', 'manual');
FIG.handles.RM = subplot('Position',[rmLeft rmBottom rmWidth rmHeight]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'YAxisLocation', 'left');
set(gca, 'XTickLabelMode', 'auto');
YLabel('dB atten for response map');
XLabel('freq (kHz)');
set(gca, 'XScale', 'log');
freqTicks = [.1 .2 .3 .5 1 2 3 5 10 20 30 50 100];
set(gca, 'XTick', 1000*freqTicks, 'XTickLabel', freqTicks, 'XMinorTick', 'off');
set(gca, 'YTickLabel', '');
FIG.RMtext = text(0.5, 0.5, sprintf('No valid RM pics\n available for this unit.'),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);


FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.02, 0.95, uinfo.text, 'Units', 'normalized', 'FontSize', FIG.fontSize-1);


%%##################################################################################################
function plot_RLs()
global uinfo FIG;

subplot(FIG.handles.RL);

uinfo.allSpont = [];
uinfo.allSpont_att = [];

legendLabels = cell(1, length(uinfo.pics.BFRL)+length(uinfo.pics.NORL));
picCounter = 0;
for currPic = uinfo.pics.BFRL
   [error x] = loadPic(currPic);
   if (~error)
      picCounter = picCounter + 1;
      checkPic(x);
      set(FIG.RLtext, 'Visible', 'Off')
      [driv spont] = compute_driv_spont_rates(x);
      if (~isnan(x.Line.attens.Tone(1,1)))
         att_dB = x.Line.attens.Tone(:,1);
      elseif (~isnan(x.Line.attens.Tone(1,2)))
         att_dB = x.Line.attens.Tone(:,2);
      end
      h_driv(picCounter) = plot(att_dB(1:size(driv,2)), trifilt(driv, 3), 'Color', 'b', 'LineWidth', 2);
      legendLabels{picCounter} = sprintf('BF (p%03d)', currPic);
      plot(att_dB(1:size(spont,2)), trifilt(spont, 3), ':b');
      uinfo.allSpont = [uinfo.allSpont spont];
      uinfo.allSpont_att = [uinfo.allSpont_att att_dB(1:size(spont,2))'];
   end
end

for currPic = uinfo.pics.NORL
   [error x] = loadPic(currPic);
   if (~error)
      picCounter = picCounter + 1;
      checkPic(x);
      set(FIG.RLtext, 'Visible', 'Off')
      [driv spont] = compute_driv_spont_rates(x);
      if (~isnan(x.Line.attens.Noise(1,1)))
         att_dB = x.Line.attens.Noise(:,1);
      elseif (~isnan(x.Line.attens.Noise(1,2)))
         att_dB = x.Line.attens.Noise(:,2);
      end
      h_driv(picCounter) = plot(att_dB(1:size(driv,2)), trifilt(driv, 3), 'Color', 'k', 'LineWidth', 2);
      legendLabels{picCounter} = sprintf('NO (p%03d)', currPic);
      plot(att_dB(1:size(spont,2)), trifilt(spont, 3), ':k');
      uinfo.allSpont = [uinfo.allSpont spont];
      uinfo.allSpont_att = [uinfo.allSpont_att att_dB(1:size(spont,2))'];
   end
end

if (picCounter > 0)
	set(gca, 'FontSize', FIG.fontSize-1);
	legend(h_driv, legendLabels, 2);
	legend('boxoff');
end

if ~isempty(uinfo.allSpont)
   meanSpont = mean(uinfo.allSpont);
   stdSpont = std(uinfo.allSpont);
   set(FIG.RLtext, 'Visible', 'On', 'Position', [1.0 1.0], 'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'bottom', 'FontSize', FIG.fontSize-1);
   set(FIG.RLtext, 'String', sprintf('SR = %.1f +/- %.1f spk/s', meanSpont, std(uinfo.allSpont)));
   plot(XLim, [meanSpont meanSpont], 'r');
   set(FIG.handles.RLrect, 'Position', [min(uinfo.allSpont_att)+1 meanSpont-(stdSpont/2) max(uinfo.allSpont_att)-2 stdSpont]);
   shad = 0.85;
   set(FIG.handles.RLrect, 'EraseMode', 'xor', 'FaceColor', [shad shad shad], 'EdgeColor', [shad shad shad]);
end


%%##################################################################################################
function plot_PSTs()
global uinfo FIG;

binWidth_sec = 1e-4; % hard-coded for now.
filt_size = 3;

pst = cell(1, size(uinfo.pics.PST,2));
i = 0;
offset = 0;
xmax = 0;
for currPic = uinfo.pics.PST
   [error x] = loadPic(abs(currPic));
   if (~error)
      checkPic(x);
      i = i + 1;
      lastBin_sec = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
      pst{i}.X = [0:binWidth_sec:lastBin_sec];
      pst{i}.Y = trifilt(hist(x.spikes{uinfo.spikeChan}(:,2), pst{i}.X), filt_size);
      if (currPic > 0)
         pst{i}.color = 'k';
      else
         pst{i}.color = 'r';  % asynchronous PST pic
      end   

      attens_dB(i) = x.Stimuli.main.attens;
      offset = max([offset max(pst{i}.Y)-min(pst{i}.Y)]);
      xmax = max([xmax lastBin_sec]);
   end
end 
offset = offset * 1.1;

if (i>0)
	subplot(FIG.handles.PST);
	set(gca, 'YGrid', 'on', 'GridLineStyle', ':');
	set(FIG.PSTtext, 'Visible', 'Off')
	% do a descending sort:
	attens_dB = -1*attens_dB;
	[attenTickLabels, sort_index] = sort(attens_dB);
	attenTickLabels = -1*attenTickLabels;
	
	YrightLabels_eval_string = ['YrightLabels = {'' '''];
	i = 0;
	for j = sort_index
		plot(pst{j}.X, pst{j}.Y + i*offset, pst{j}.color);
		YrightLabels_eval_string = [YrightLabels_eval_string ' ; ''' sprintf('p%03d',abs(uinfo.pics.PST(j))) ''''];
		i = i + 1;
	end
	YrightLabels_eval_string = [YrightLabels_eval_string '};'];
	eval(YrightLabels_eval_string);
	
	attenTicks = [-1*offset:offset:(i-1)*offset];
	set(gca, 'YTick', attenTicks, 'YTickLabel', [inf attenTickLabels], 'YMinorTick', 'off');
	YLim([-0.1*offset (i+0.1)*offset]);
   XLim([0 xmax]);
	
	subplot(FIG.handles.PSTright);
	set(gca, 'YTick', attenTicks, 'YTickLabel', YrightLabels, 'YMinorTick', 'off');
	YLim([-0.1*offset (i+0.1)*offset]);

end

%%##################################################################################################
function plot_RMs()
global uinfo FIG;

filt_size = 3;

rm = cell(1, size(uinfo.pics.RM,2));
i = 0;
offset = 0;
for currPic = uinfo.pics.RM
   [error x] = loadPic(abs(currPic));
   if (~error)
      checkPic(x);
      i = i + 1;
      [driv spont] = compute_driv_spont_rates(x);
      [rm{i}.freq_Hz sortIndex] = sort(x.Line.freq(1:x.Stimuli.fully_presented_lines));
      rm{i}.driv = trifilt(driv(sortIndex), filt_size);
      rm{i}.spont = trifilt(spont(sortIndex), filt_size);
      if (currPic > 0)
         rm{i}.color = 'b';
      else
         rm{i}.color = 'r';  % 2-tone RM pic
      end   
      attens_dB(i) = x.Stimuli.main.attens;
      offset = max([offset max(rm{i}.driv) max(rm{i}.driv)-min(rm{i}.driv)]);
      offset = max([offset max(rm{i}.spont) max(rm{i}.spont)-min(rm{i}.spont)]);
   end
end 
offset = offset * 1.1;

if (i>0)
   
	subplot(FIG.handles.RM);
   set(gca, 'YGrid', 'on', 'GridLineStyle', ':');
   set(FIG.RMtext, 'Visible', 'Off')

   % add a vertical line for the unit BF:
   if ~isempty(uinfo.BF_kHz)
      h_bfLine = plot([1000*uinfo.BF_kHz 1000*uinfo.BF_kHz], [0 100], 'Color', 0.6*[1 1 0]);
   end
   
	% do a descending sort:
	attens_dB = -1*attens_dB;
	[attenTickLabels, sort_index] = sort(attens_dB);
	attenTickLabels = -1*attenTickLabels;
	
	
   YrightLabels_eval_string = ['YrightLabels = {'' '];
	i=0;
	min_freq_Hz = 1e6;  % arbitrarily large...
	max_freq_Hz = 0;
   prev_attens_dB = -1;
   nRepeats = 0;
	for j = sort_index
       curr_attens_dB = attens_dB(j);
       rptFlag = 0;
       if (curr_attens_dB == prev_attens_dB)
          rptFlag = 1;
          nRepeats = nRepeats + 1;
       end
       prev_attens_dB = curr_attens_dB;
       plot(rm{j}.freq_Hz, rm{j}.driv + (i-nRepeats)*offset, rm{j}.color);
       plot(rm{j}.freq_Hz, rm{j}.spont + (i-nRepeats)*offset, 'k', 'LineWidth', 2);
       if ~(rptFlag)
          YrightLabels_eval_string = [YrightLabels_eval_string ''' ; ''' sprintf('p%03d',abs(uinfo.pics.RM(j))) ];
       else
          YrightLabels_eval_string = [YrightLabels_eval_string sprintf(',%d',abs(uinfo.pics.RM(j))) ];
       end
       i = i + 1;
       attenTickLabels(i-nRepeats) = attenTickLabels(i);
       min_freq_Hz = min([min_freq_Hz min(rm{j}.freq_Hz)]);
       max_freq_Hz = max([max_freq_Hz max(rm{j}.freq_Hz)]);
	end
	YrightLabels_eval_string = [YrightLabels_eval_string '''};'];
   eval(YrightLabels_eval_string);

	attenTicks = [-1*offset:offset:((i-nRepeats)-1)*offset];
	set(gca, 'YTick', attenTicks, 'YTickLabel', [inf attenTickLabels], 'YMinorTick', 'off');
   ymin = -0.1*offset;
   ymax = ((i-nRepeats)+0.1)*offset;
	YLim([ymin ymax]);
	XLim([0.9*min_freq_Hz 1.1*max_freq_Hz]);

   if ~isempty(uinfo.BF_kHz)
      set(h_bfLine, 'YData', [ymin ymax]);
   end

   % add a scale bar:
   barLength = 0.5*10^(floor(log10(offset)));  
	while (barLength < 0.35*offset)
       barLength = barLength * 2;
	end
   shad = 0.7;
   plot([max_freq_Hz max_freq_Hz], [(ymax-0.1*offset)-barLength ymax-0.1*offset],...
          'Color', [shad shad shad], 'LineWidth', 3);
   set(FIG.RMtext, 'String', sprintf('scale bar = %.1f spk/s', barLength));
   set(FIG.RMtext, 'Visible', 'On', 'Position', [1.0 1.0], 'HorizontalAlignment', 'right', ...
                     'VerticalAlignment', 'bottom', 'FontSize', FIG.fontSize-1);

   
   % add picture number labels on right side:
   subplot(FIG.handles.RMright);
 	set(gca, 'YTick', attenTicks, 'YTickLabel', YrightLabels, 'YMinorTick', 'off');
	YLim([ymin ymax]);

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


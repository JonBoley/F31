function rc = usmp(trackNum, unitNum, data_folder)

% Written by GE 19Mar2004.


global uinfo FIG Data;

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

plot_eh_smp;

%%##################################################################################################
function error = get_unit_info()
global uinfo Data;

error = 0;

uinfo.spikeChan = 1;

% Temporarily switch to experiment data directory and load "info file":
origDir = cd; eval (['cd ''' uinfo.dataFolder '''']); load DataInfoFile; eval (['cd ''' origDir '''']);

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
          
     % /eh/ vowel smp info:
uinfo.SMP.EH = Data.Vowels_SMP{uinfo.TNum,uinfo.UNum}.EH;

%%##################################################################################################
function layout_figure()
global uinfo FIG;

FIG.fontSize = 8;
FIG.handles.main = figure(3); clf; hold on;  % ge debug 19Mar2004, figure number hard-coded, for now
set(gcf, 'Name', 'vowel smp analysis (usmp.m)');

FIG.handles.smp_eh = subplot('Position',[0.1 0.1 0.8 0.8]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse');
YLabel('rate (spk/s)');
set(gca, 'YTickLabel', '', 'XTickLabel', '');
FIG.smp_eh_text = text(0.5, 0.5, sprintf('No valid vowel /eh/ pics\n available for this unit.'),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);

FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.02, 0.95, uinfo.text, 'Units', 'normalized', 'FontSize', FIG.fontSize-1);


%%##################################################################################################
function plot_eh_smp()
global uinfo FIG Data;

smp_eh = cell(1, size(uinfo.SMP.EH, 2));
picCounter = 0;
for i = 1:size(uinfo.SMP.EH, 2)
   currPic = uinfo.SMP.EH{i}.picNUM;
   if (abs(log2(uinfo.SMP.EH{i}.featureFreq_Hz) - log2(1000*uinfo.BF_kHz)) < 1/8)
       [error x] = loadPic(currPic);
       if (~error & isempty(x.Stimuli.bad_lines)) % ge debug 26Mar2004 -- temp fix -- should correct for bad lines.
          picCounter = picCounter + 1;
          smp_eh{picCounter}.picNum = currPic;
          smp_eh{picCounter}.feature  = uinfo.SMP.EH{i}.feature;
          [smp_eh{picCounter}.driv smp_eh{picCounter}.spont] = compute_driv_spont_rates(x);
          smp_eh{picCounter}.att_dB = x.Stimuli.attens;
       end
   end
end

if (picCounter > 0)
   subplot(FIG.handles.smp_eh);
   set(FIG.smp_eh_text, 'Visible', 'Off');
   set(gca, 'YTickLabelMode', 'auto', 'XTickLabelMode', 'auto');

   
	% Parameters specific for vowel /eh/:
	feature_names = {'T0' 'F1' 'T1' 'F2' 'T2' 'F3' 'T3'};
   feature_base_freq_Hz = [300 501 1202 1703 2203 2504 3005];
	feature_levels_dB = [-11.2 0 -30.6 -15.5 -33.6 -28.7 -42.7];
	feature_line_colors = ['k' 'r' 'r' 'g' 'g' 'b' 'b'];

   
   % Plot driven rates first:
   legendLabels = cell(1, picCounter);
	for i = 1:picCounter
      feature_ind = strcmp(smp_eh{i}.feature(1),'T') + 2*(str2num(smp_eh{i}.feature(2)));
      curr_feature_level_dB = feature_levels_dB(feature_ind);
      curr_feature_color = feature_line_colors(feature_ind);
      if strcmp(smp_eh{i}.feature(1),'F')  % i.e., is a formant.
         curr_feature_lineWidth = 2;
         curr_feature_style = '-';
      else
         curr_feature_lineWidth = 1;
         curr_feature_style = '-';
      end
      plot(smp_eh{i}.att_dB(1:size(smp_eh{i}.driv,2))-curr_feature_level_dB, trifilt(smp_eh{i}.driv, 3), ...
         'Color', curr_feature_color, ...
         'LineStyle', curr_feature_style, ...
         'LineWidth', curr_feature_lineWidth);
      legendLabels{i} = sprintf('%s (p%03d)', smp_eh{i}.feature, smp_eh{i}.picNum);
   end
   set(gca, 'FontSize', FIG.fontSize-1);
   legend(legendLabels, 2);
   legend('boxoff');

   % Plot spontaneous rates:
	for i = 1:picCounter
      feature_ind = strcmp(smp_eh{i}.feature(1),'T') + 2*(str2num(smp_eh{i}.feature(2)));
      curr_feature_level_dB = feature_levels_dB(feature_ind);
      plot(smp_eh{i}.att_dB(1:size(smp_eh{i}.spont,2))-curr_feature_level_dB, trifilt(smp_eh{i}.spont, 3), 'k:');
	end
   
   currYLim = get(gca, 'YLim');
   ylim([-3 currYLim(2)]);

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

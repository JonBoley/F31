function rc = urss(trackNum, unitNum, data_folder)

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

layout_figure(2, 'combined RSS analyses');
plot_RSS_combined_analysis;

layout_figure(4, 'separate RSS analyses', 0.45);
plot_RSS_separate_analysis;

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
          
     % RSS info:
uinfo.RSS = Data.RSS{uinfo.TNum,uinfo.UNum}.RSS;
if (size(uinfo.RSS,2) < 1)
   error_report = 'No RSS pictures are available.'
   error = -1;
   return
end   


%%##################################################################################################
function layout_figure(figNum, figTag, horizOffset_norm)
global uinfo FIG;

if ~exist('horizOffset_norm')
   horizOffset_norm = 0;
end

FIG.fontSize = 8;

if isempty(find(findobj==figNum))
	FIG.handles.main = figure(figNum); clf; hold on;  % ge debug 19Mar2004, figure number hard-coded, for now
	set(gcf, 'Units', 'normalized');
	set(gcf, 'Position', [0.05+horizOffset_norm 0.5 0.4 0.4]);
    set(gcf, 'Color', 'w');
else
    figure(figNum); clf; hold on;
end
set(gcf, 'Name', ['RSS analysis (urss.m) ' figTag]);

rss1Left = 0.1; rss1Bottom = 0.1; rss1Width = 0.7; rss1Height = 0.78;
FIG.handles.rss1right = subplot('Position',[rss1Left+rss1Width rss1Bottom 0.01 rss1Height]); % for labelling right axis
set(gca, 'FontSize', FIG.fontSize-1);
set(gca, 'XTickMode', 'manual');
set(gca, 'YAxisLocation', 'right');
set(gca, 'YTickMode', 'manual');
FIG.handles.RSS_first_order = subplot('Position',[rss1Left rss1Bottom rss1Width rss1Height]);
hold on;
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'YAxisLocation', 'left');
set(gca, 'XTickLabelMode', 'manual');
YLabel('dB atten for RSS stim');
XLabel('frequency (kHz)');
set(gca, 'XScale', 'log');
set(gca, 'YTickLabel', '');
YLim([-5 5]);
FIG.RSStext = text(0.5, 0.5, sprintf('No valid RSS pics\n available for this unit.'),...
       'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', FIG.fontSize-1);

FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.02, 0.95, uinfo.text, 'Units', 'normalized', 'FontSize', FIG.fontSize-1);


%%##################################################################################################
function plot_RSS_combined_analysis()
global uinfo FIG Data;

% Determine number of RSS computations to be done:
attenVect = uinfo.RSS{1}.dBatten;
urVect = uinfo.RSS{1}.updateRate_Hz;
nRSSgroups = 1;
for i = 2:size(uinfo.RSS, 2)
   currAtten = uinfo.RSS{i}.dBatten;
   curr_ur = uinfo.RSS{i}.updateRate_Hz;
   newAttenFlag = 0;
   attenIndices = find(attenVect == currAtten);
   if isempty(attenIndices)    % current attenuation has not yet been encountered.
      nRSSgroups = nRSSgroups + 1;
      attenVect(nRSSgroups) = currAtten;
      urVect(nRSSgroups) = curr_ur;
   elseif isempty(find(urVect(attenIndices)==curr_ur))% new update rate for a previously encountered attenuation.
      noteString = 'New update rate encountered.'
      nRSSgroups = nRSSgroups + 1;
      attenVect(nRSSgroups) = currAtten;
      urVect(nRSSgroups) = curr_ur;
   end
end
% Do a descending sort on the attenuations vector:
attenVect = -1*attenVect;
attenVect = sort(attenVect);
attenVect = -1*attenVect;


% Do the RSS analyses:
rss = cell(1, nRSSgroups);
picString = cell(1, nRSSgroups+1);
picString{1} = [''];
offset = 0;
ymin = 0;
nInvalidGroups = 0;
for j = 1:nRSSgroups
   currAtten = attenVect(j);
   curr_ur = urVect(j);
   full_spectra_matrix = [];
   full_driv_vector = [];
   c0_indices = [];
   prev_binFreqs_kHz = [];
   picString{j+1-nInvalidGroups} = ['p'];
   validPicsFound = 0;
   for i = 1:size(uinfo.RSS, 2)
      if (uinfo.RSS{i}.dBatten == currAtten) & (uinfo.RSS{i}.updateRate_Hz == curr_ur)
			currPic = uinfo.RSS{i}.picNUM;
         picString{j+1-nInvalidGroups} = [picString{j+1-nInvalidGroups} sprintf('%d',currPic) ','];
			stimClass = uinfo.RSS{i}.StimClass;
			stimSet = uinfo.RSS{i}.StimSet;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %    ge debug 23Mar2004: hard-coded for now -- load RSS synthesis info manually:
% MH26Apr2004               stimDir = sprintf('\\Signals\\GE_MH\\RSS_stim%s_set%02d\\', stimClass, stimSet);
              stimDir = sprintf('..\\..\\Signals\\GE_MH\\RSS_stim%s_set%02d\\', stimClass, stimSet);
              origDir = cd;  eval (['cd ''' stimDir '''']); load synthesis_info; eval (['cd ''' origDir '''']);
              nRSS = size(spectra_matrix, 2);
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if ~isempty(prev_binFreqs_kHz)
            if ~(binFreqs_kHz == prev_binFreqs_kHz)
               warning = 'Mismatch in bin frequencies of RSS sets.!!!!'
            end
            prev_binFreqs_kHz = binFreqs_kHz;
         end
			[error x] = loadPic(currPic);
			if (~error & ~(x.Stimuli.fully_presented_lines < nRSS))
            validPicsFound = 1;
            full_spectra_matrix = [full_spectra_matrix spectra_matrix];
            [driv spont] = compute_driv_spont_rates(x);
            c0_indices = [c0_indices [1 2 51 52]+length(full_driv_vector)];
            full_driv_vector = [full_driv_vector driv];
			end
      end
   end
	rss{j}.rescaled_freqVector = (curr_ur / 100e3) * binFreqs_kHz;
	picString{j+1-nInvalidGroups} = picString{j+1-nInvalidGroups}(1:end-1);
   if (validPicsFound)
       c0 = mean(full_driv_vector(c0_indices));
       driv_norm = full_driv_vector - c0;
       var_full_spectra_matrix = mean(diag(full_spectra_matrix*(full_spectra_matrix')));
       rss{j}.weight_vect = (full_spectra_matrix * driv_norm') / var_full_spectra_matrix;
       offset = max([offset max(rss{j}.weight_vect) max(rss{j}.weight_vect)-min(rss{j}.weight_vect)]);
       ymin = -1*max([-1*ymin -1*min(rss{j}.weight_vect)]);
   else
       nInvalidGroups = nInvalidGroups + 1;
       rss{j}.weight_vect = [];
   end
end
offset = offset * 1.2;
ymin = ymin * 1.1;
ymax = ((nRSSgroups-nInvalidGroups)+0.1)*offset;

% Do some preliminary plot modifications:
subplot(FIG.handles.RSS_first_order);
set(gca, 'YGrid', 'on', 'GridLineStyle', ':');
freqTicks = [.1 .2 .5 1 2 5 10 20 50];
set(gca, 'XTick', freqTicks, 'XTickLabel', freqTicks, 'XMinorTick', 'off');
set(FIG.RSStext, 'Visible', 'Off')
if ~isempty(uinfo.BF_kHz)
   h_bfLine = plot([uinfo.BF_kHz uinfo.BF_kHz], [0 100], 'Color', 0.6*[1 1 0]);
end
YLim([ymin ymax]);
attenTicks = [-10^6 0:offset:(nRSSgroups-1)*offset];
set(gca, 'YTick', attenTicks-ymin, 'YTickLabel', [inf attenVect], 'YMinorTick', 'off');
subplot(FIG.handles.rss1right);
set(gca, 'YTick', attenTicks-ymin, 'YTickLabel', picString, 'YMinorTick', 'off');
YLim([ymin ymax]);

% Plot the RSS weight vectors:
subplot(FIG.handles.RSS_first_order);
for j = 1:nRSSgroups
   if ~isempty(rss{j}.weight_vect);
      plot(rss{j}.rescaled_freqVector, rss{j}.weight_vect + (j-1)*offset - ymin, 'k', 'LineWidth', 2);
   end
end

% Add a scale bar:
barLength = 0.5*10^(floor(log10(offset)));  
while (barLength < 0.35*offset)
    barLength = barLength * 2;
end
shad = 0.7;
xlims = XLim;
plot([xlims(2) xlims(2)], [(ymax-0.1*offset)-barLength ymax-0.1*offset],...
       'Color', [shad shad shad], 'LineWidth', 3);
set(FIG.RSStext, 'String', sprintf('scale bar = %.1f (spk/s)/dB', barLength));
set(FIG.RSStext, 'Visible', 'On', 'Position', [1.0 1.0], 'HorizontalAlignment', 'right', ...
                  'VerticalAlignment', 'bottom', 'FontSize', FIG.fontSize-1);



%%##################################################################################################
function plot_RSS_separate_analysis()
global uinfo FIG Data;

rss = cell(1, size(uinfo.RSS, 2));
offset = 0;
ymin = 0;
picCounter = 0;
for i = 1:size(uinfo.RSS, 2)
   currPic = uinfo.RSS{i}.picNUM;
   stimClass = uinfo.RSS{i}.StimClass;
   stimSet = uinfo.RSS{i}.StimSet;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %    ge debug 23Mar2004: hard-coded for now -- load RSS synthesis info manually:
% MH26Apr2004       stimDir = sprintf('\\Signals\\GE_MH\\RSS_stim%s_set%02d\\', stimClass, stimSet);
      stimDir = sprintf('..\\..\\Signals\\GE_MH\\RSS_stim%s_set%02d\\', stimClass, stimSet);
    origDir = cd;  eval (['cd ''' stimDir '''']); load synthesis_info; eval (['cd ''' origDir '''']);
      nRSS = size(spectra_matrix, 2);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [error x] = loadPic(currPic);
   if (~error & ~(x.Stimuli.fully_presented_lines < nRSS))
      picCounter = picCounter + 1;
      [driv spont] = compute_driv_spont_rates(x);
      c0 = mean(driv([1 2 51 52]));
      driv_norm = driv - c0;
      var_spectra_matrix = mean(diag(spectra_matrix*(spectra_matrix')));
      rss{picCounter}.weight_vect = (spectra_matrix * driv_norm') / var_spectra_matrix;
      rss{picCounter}.rescaled_freqVector = (x.Stimuli.updateRate_Hz / 100e3) * binFreqs_kHz;
      attens_dB(picCounter) = x.Stimuli.attens;
      offset = max([offset max(rss{picCounter}.weight_vect)-min(rss{picCounter}.weight_vect)]);
      ymin = -1*max([-1*ymin -1*min(rss{picCounter}.weight_vect)]);
   else
      picCounter = picCounter + 1;
      attens_dB(picCounter) = -1;
   end
end
offset = offset * 1.1;
ymin = ymin * 1.1;

if (picCounter > 0) & (~isempty(find(attens_dB>0)))
   subplot(FIG.handles.RSS_first_order);
   set(gca, 'YGrid', 'on', 'GridLineStyle', ':');
   freqTicks = [.1 .2 .5 1 2 5 10 20 50];
   set(gca, 'XTick', freqTicks, 'XTickLabel', freqTicks, 'XMinorTick', 'off');
   set(FIG.RSStext, 'Visible', 'Off')
   
   % add a vertical line for the unit BF:
   if ~isempty(uinfo.BF_kHz)
      h_bfLine = plot([uinfo.BF_kHz uinfo.BF_kHz], [0 100], 'Color', 0.6*[1 1 0]);
   end
   
 	% do a descending sort:
	attens_dB = -1*attens_dB;
	[attenTickLabels, sort_index] = sort(attens_dB);
   attens_dB = -1* attens_dB;
	attenTickLabels = -1*attenTickLabels;

   YrightLabels_eval_string = ['YrightLabels = {'''];
	i = 0;
   prev_attens_dB = -1;
   nRepeats = 0;
	for j = sort_index
		curr_attens_dB = attens_dB(j);
      if ~(curr_attens_dB < 0)
			rptFlag = 0;
			if (curr_attens_dB == prev_attens_dB)
				rptFlag = 1;
				nRepeats = nRepeats + 1;
			end
			prev_attens_dB = curr_attens_dB;
         i = i + 1;
         plot(rss{j}.rescaled_freqVector, rss{j}.weight_vect + ((i-nRepeats)-1)*offset - ymin, 'k', 'LineWidth', 2);
			if ~(rptFlag)
				YrightLabels_eval_string = [YrightLabels_eval_string ''' ; ''' sprintf('p%d',abs(uinfo.RSS{j}.picNUM)) ];
			else
				YrightLabels_eval_string = [YrightLabels_eval_string sprintf(',%d',abs(uinfo.RSS{j}.picNUM)) ];
			end
         attenTickLabels(i-nRepeats) = attenTickLabels(i);
      end
	end
	YrightLabels_eval_string = [YrightLabels_eval_string '''};'];
   eval(YrightLabels_eval_string);

   ymax = ((picCounter-nRepeats)+0.1)*offset;

 	attenTicks = [-10^6 0:offset:((i-nRepeats)-1)*offset];
	set(gca, 'YTick', attenTicks-ymin, 'YTickLabel', [inf attenTickLabels], 'YMinorTick', 'off');
	YLim([ymin ymax]);

   if ~isempty(uinfo.BF_kHz)
      set(h_bfLine, 'YData', [ymin ymax]);
   end

   % Add a scale bar:
	barLength = 0.5*10^(floor(log10(offset)));  
	while (barLength < 0.35*offset)
        barLength = barLength * 2;
	end
	shad = 0.7;
	xlims = XLim;
	plot([xlims(2) xlims(2)], [(ymax-0.1*offset)-barLength ymax-0.1*offset],...
           'Color', [shad shad shad], 'LineWidth', 3);
	set(FIG.RSStext, 'String', sprintf('scale bar = %.1f (spk/s)/dB', barLength));
	set(FIG.RSStext, 'Visible', 'On', 'Position', [1.0 1.0], 'HorizontalAlignment', 'right', ...
                      'VerticalAlignment', 'bottom', 'FontSize', FIG.fontSize-1);

   subplot(FIG.handles.rss1right);
 	set(gca, 'YTick', attenTicks-ymin, 'YTickLabel', YrightLabels, 'YMinorTick', 'off');
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

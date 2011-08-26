function vowels_smp(picNum, figNum)

if ~exist('figNum', 'var')
   figNum = 200;
end
params.spikeChan = 1;

vowelLevels_dBatten = [90 70 50];
vowelLevels_colors = ['krb'];

%feature_names = {'T0' 'F1' 'T1' 'F2' 'T2' 'F3' 'T3'};
feature_levels_dB = [-11.2 0 -30.6 -15.5 -33.6 -28.7 -42.7];
line_colors = ['k' 'r' 'r' 'g' 'g' 'b' 'b'];


figure(figNum); clf
set(gcf, 'Name', sprintf ('pictures p%04d - p%04d; ''eh'' smp analysis', picNum(1), picNum(end)));
h_normLevels = subplot('Position',[0.1 0.6 0.85 0.35]);
hold on;
h_vowelLevels = subplot('Position',[0.1 0.1 0.85 0.4]);
hold on;

for PICind=1:length(picNum)
   [rc x] = loadPic(picNum(PICind));
   
   feature_ind = strcmp(x.Stimuli.feature(1),'T') + 2*(str2num(x.Stimuli.feature(2)));
   curr_feature_level_dB = feature_levels_dB(feature_ind);
   curr_feature_color = line_colors(feature_ind);
   if strcmp(x.Stimuli.feature(1),'F')
      curr_feature_style = '-';
   else
      curr_feature_style = ':';
   end
   
   [att_dB driv spont] = calc_RL(x, params);
   
   subplot(h_normLevels);
   plot(att_dB(1:size(driv,2))-curr_feature_level_dB, trifilt(driv, 3), ...
      'Color', curr_feature_color, 'LineStyle', curr_feature_style);
   plot(att_dB(1:size(spont,2))-curr_feature_level_dB, trifilt(spont, 3), ...
      'Color', 'y', 'LineStyle', '-');
       
   subplot(h_vowelLevels);
   color_index = 1;
   for (currLevel = vowelLevels_dBatten)
      feature_rate = driv(find(att_dB == currLevel));
      scatter(currLevel + curr_feature_level_dB, feature_rate, 4, vowelLevels_colors(color_index));
      color_index = color_index + 1;
   end
   
end

subplot(h_normLevels);
set(gca,'XDir','reverse');
ylabel('firing rate (spk/s)');
yChunks = 50;
currXLim = get(gca, 'XLim');
currYLim = get(gca, 'YLim');
ylim([0 yChunks*ceil(currYLim(2)/yChunks)]);

subplot(h_vowelLevels);
xlabel('feature level (dB vowel\_atten + atten\_re\_F1)');
xlim(currXLim);
set(gca,'XDir','reverse');
ylim([0 yChunks*ceil(currYLim(2)/yChunks)]);



%%#######################################################################
function [rc, x] = loadPic(picNum)     % Load picture
picSearchString = sprintf('p%04d*.m', picNum);
picMFile = dir(picSearchString);
if (~isempty(picMFile))
   rc = 0;
   eval( strcat('x = ',picMFile.name(1:length(picMFile.name)-2),';') );
else
   error = sprintf('Picture file p%04d*.m not found.', picNum)
   rc = -1;
   x = [];
   return;
end

%%#######################################################################
function [att_dB, driv, spont] = calc_RL(x, params)
[driv spont] = compute_driv_spont_rates(x, params);
att_dB = x.Stimuli.attens;

%%#######################################################################
function [driv, spont] = compute_driv_spont_rates(x, params)
     % compute driven and spontaneous rates for each line
spikeTimes = x.spikes{params.spikeChan};
driven_dur_sec = x.Hardware.Trigger.StmOn / 1000;
spont_dur_sec = x.Hardware.Trigger.StmOff / 1000;
for line = 1:x.Stimuli.fully_presented_lines
   spikeIndices = find( (spikeTimes(:,1) == line ) );
   driv(line) = length (find(spikeTimes(spikeIndices,2) <= driven_dur_sec));
   spont(line) = length(spikeIndices) - driv(line);
end
driv = driv / driven_dur_sec;
spont = spont / spont_dur_sec;

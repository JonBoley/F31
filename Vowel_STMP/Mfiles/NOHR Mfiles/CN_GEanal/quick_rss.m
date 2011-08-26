function quick_rss(freqVector, stimMatrix, weight_threshold, picNum, lineStyle)


params = [];
params.spikeChan = 1;
params.figNum = 100;

[rc x] = loadPic(picNum);
[driv spont] = compute_driv_spont_rates(x, params);

c0 = mean(driv([1 2 51 52]));  % average response to flat noise.
driv_norm = driv - c0;

weight_vect = (stimMatrix * driv_norm') / 9200;
figure(33);
nameStr = sprintf ('picture p%04d weighting function', picNum);
set(gcf, 'Name', nameStr);
rescaled_freqVector = (x.Stimuli.updateRate_Hz / 100e3) * freqVector;
 semilogx(rescaled_freqVector, weight_vect, lineStyle);
%semilogx(rescaled_freqVector, weight_vect);
freqTicks = [1 2 3 4 5 7 10 20 30 40 50];
set(gca, 'XTick', freqTicks, 'XTickLabel', freqTicks, 'XMinorTick', 'off');
YLim([-5 5]);
rescaled_freqVector(find(weight_vect > weight_threshold))';

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

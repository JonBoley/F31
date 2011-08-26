function params = picview(params)
% function params = picview(params)
%  Use zero input arguments to initialize params.

% Set up default input values:
if (nargin == 0)
   clear params;
   params = [];
   params.plotType = '';
   params.picArray = {};
   params.spikeChan = 1;
   params.figNum = 100;
   params.forceHold = 0;
   params.extern.func = '';
   params.checkPicType = 1;
   return;
end

% color definitions:
cSilver = [0.75 0.75 0.75]; cGray = [0.3 0.3 0.3]; cMaroon = [0.5 0 0]; cGreen = [0 0.5 0]; cNavy = [0 0 0.5];
cPurple = [0.5 0 0.5]; cOlive = [0.5 0.5 0]; cTeal = [0 0.5 0.5]; cBlack = [0 0 0]; cRed = [1 0 0];
cLime = [0 1 0]; cBlue = [0 0 1]; cMagenta = [1 0 1]; cYellow = [1 1 0]; cCyan = [0 1 1];
cMap = [cBlack; cRed; cBlue; cOlive; cMagenta; cLime; cPurple; cCyan; cTeal];
nColors = size(cMap, 1);

% Loop through pictures in picArray:
if (params.forceHold == 0)
    figure(params.figNum);
    hold off;
end
params.multi.offset = 0;
params.multi.newBlock = 0;
for iBlock = 1:size(params.picArray,1)
   params.multi.lineParams.color = cMap(mod(iBlock, nColors),:);
   for iCol = 1:size(params.picArray{iBlock},2)
      params.multi.picNum = params.picArray{iBlock}(1,iCol);
      params = picview_single(params);
      hold on;
      params.multi.newBlock = 0;
   end
   params.multi.newBlock = 1;
end


%%#######################################################################
function params = picview_single(params)
% function params = picview_single(plotType, picNum, params, figNum)

rc = 0;
[rc x] = loadPic(params.multi.picNum);
if (rc ~= 0)
   return
end

switch (params.plotType)
case 'raster' % raster plot of spike times vs line number
   figDefaultSettings(params.figNum);
   plot(x.spikes{1}(:,2),x.spikes{1}(:,1), 'k.', 'MarkerSize', 4);
   trig = x.Hardware.Trigger;
   XLim([0 (trig.StmOn+trig.StmOff)/1000]);
case 'rl'  % rate-level
   if (params.checkPicType==0)
      [att_dB driv spont] = calc_RL(x, params);
   elseif(checkPicType(params.multi.picNum, '', 'RLV', 1) == 0)   % tone rate-level
      [att_dB driv spont] = calc_RL(x, params);
   elseif (checkPicType(params.multi.picNum, 'rate-level', 'NO', 1) == 0) % noise rate-level
      [att_dB driv spont] = calc_RL(x, params);
   else
      return
   end
   figDefaultSettings(params.figNum);
   plot(att_dB(1:size(driv,2)), trifilt(driv, 3), 'Color', 'r', 'LineWidth', 2);
   hold on;
   plot(att_dB(1:size(spont,2)), trifilt(spont, 3), 'Color', 'k');
   hold off;
   set(gca,'XDir','reverse');
case 'rm' % response map
   if(checkPicType(params.multi.picNum, 'response map', 'RM', 1) == 0)
      [freq_kHz driv spont] = calc_RM(x, params);
   else
      return
   end
   figDefaultSettings(params.figNum);
   plot(freq_kHz, driv, freq_kHz, spont);
   set(gca, 'XScale', 'log');
case 'pst' % peri-stimulus time histogram
   if (params.checkPicType==0)
      [params, pst_X, pst_Y] = calc_PST(x, params);
   elseif(checkPicType(params.multi.picNum,'pst','PST*',1)==0)
      [params, pst_X, pst_Y] = calc_PST(x, params);
   else
      return
   end
   figDefaultSettings(params.figNum);
   plot(pst_X, pst_Y, params.multi.lineParams);
   if (params.PST.YLimManual ~= 0)
      set(gca, 'YLim', params.PST.YLim);
   end
case 'ep'   % acquired A/D signal ("evoked potential")
	if (~isfield(x.EP,'aveY'))
       error = sprintf('No EP data available in picture p%04d.', params.multi.picNum) 
       return;
	else
       [params ep_X ep_Y] = calc_EP(x, params);
       figDefaultSettings(params.figNum);
       plot(ep_X, ep_Y, params.multi.lineParams);
	end
case 'cap'   % acquired A/D signal ("compound action potential")
	if (~isfield(x.CAPData,'CAPAvg_V'))
       error = sprintf('No CAP data available in picture p%04d.', params.multi.picNum) 
       return;
	else
       [params cap_X cap_Y] = calc_CAP(x, params);
       figDefaultSettings(params.figNum);
       plot(cap_X, cap_Y, 'k');
       set(gca, params.CAP.axis);
       xlabel(params.CAP.XLabel);
       ylabel(params.CAP.YLabel);
	end
case 'extern'   % external analysis function
    try
        [params extern_X extern_Y] = eval([params.extern.func '(x, params);']);
        figDefaultSettings(params.figNum);
        plot(extern_X, extern_Y, params.extern.lineParams);
    catch
        error = sprintf('Could not properly evaluate external function.')
    end
otherwise
   error = sprintf ('%s is not a valid plot type.', params.plotType)
   return
end

updateFigName(params.multi.picNum);

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
function figDefaultSettings(figNum)
figure(figNum);
set(gca, 'XDir', 'normal');
set(gca, 'XScale', 'linear');

%%#######################################################################
function updateFigName(picNum)
nameStr = sprintf ('picture p%04d', picNum);
set(gcf, 'Name', nameStr);


%%#######################################################################
function rc = checkPicType(picNum, picType, fileSuffix, rprtErr)     % check filename
picSearchString = sprintf('p%04d*%s.m', picNum, fileSuffix);
picMFile = dir(picSearchString);
if (~isempty(picMFile))
   rc = 0;
else
   if (rprtErr ~= 0)
      warning = sprintf('Picture file p%04d*.m is not an official ''%s'' file.', picNum, picType)
   end
   rc = -1;
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

%%#######################################################################
function [att_dB, driv, spont] = calc_RL(x, params)
[driv spont] = compute_driv_spont_rates(x, params);
% Set an appropriate attenuations vector:
if (isfield(x.Line.attens, 'Tone'))
    if (~isnan(x.Line.attens.Tone(1,1)))
        att_dB = x.Line.attens.Tone(:,1);
    elseif (~isnan(x.Line.attens.Tone(1,2)))
        att_dB = x.Line.attens.Tone(:,2);
    end
elseif (isfield(x.Line.attens, 'Noise'))
    if (~isnan(x.Line.attens.Noise(1,1)))
        att_dB = x.Line.attens.Noise(:,1);
    elseif (~isnan(x.Line.attens.Noise(1,2)))
        att_dB = x.Line.attens.Noise(:,2);
    end
else
    att_dB = x.Stimuli.attens;
end

%%#######################################################################
function [freq_kHz, driv, spont] = calc_RM(x, params)
[driv spont] = compute_driv_spont_rates(x, params);
[freq_kHz sortIndex] = sort(x.Line.freq(1:x.Stimuli.fully_presented_lines));
driv = driv(sortIndex);
spont = spont(sortIndex);
[maxDriv maxIndex] = max(driv);
info = sprintf ('\tMaximum: %d/sec at %.1f kHz.', maxDriv, freq_kHz(maxIndex))

%%#######################################################################
function [params, pst_X, pst_Y] = calc_PST(x, params)
if (~isfield(params,'PST'))
    params.PST.firstBin_sec = 0;
    params.PST.lastBin_sec = 0.999;
    params.PST.binWidth_sec = 0.005;
    params.PST.pictSpacing = 10;
    params.PST.blockSpacing = 0;
    params.PST.nullThreshold = -1;
    params.PST.nullMask = {};
    params.PST.YLimManual = 0;
    params.PST.YLim = [0 100];
    params.PST.spontSubtract = 0;
end
multi = params.multi;
multi.offset = multi.offset + params.PST.blockSpacing*multi.newBlock;
firstBin = params.PST.firstBin_sec;
lastBin = params.PST.lastBin_sec;
binWidth = params.PST.binWidth_sec;
pst_X = [firstBin:binWidth:lastBin];
pst_Y = pst_X;
pst_Y(:) = 0;
lastBinNumber = size(pst_X, 2) - 1;
for i=1:size(x.spikes{params.spikeChan}, 1)
    spikeTime_sec = x.spikes{params.spikeChan}(i,2);
    bNullSpike = 0;
    for j = 1:size(params.PST.nullMask, 1)
        nullStart_sec = params.PST.nullMask{j}(1) - params.PST.nullMask{j}(2);
        nullStop_sec = params.PST.nullMask{j}(1) + params.PST.nullMask{j}(2);
        if (spikeTime_sec>nullStart_sec & spikeTime_sec<nullStop_sec)
            bNullSpike = 1;
        end
    end
    if (bNullSpike ~= 0)
        continue;
    end
    binNumber = floor((spikeTime_sec - firstBin) / binWidth) + 1;
    if (binNumber<1 | binNumber>(size(pst_X, 2)))  % spike time is not in an available bin.
       continue;
    else   % add a count to the appropriate bin.
        pst_Y(binNumber) = pst_Y(binNumber) + 1;
    end
end
pst_Y = pst_Y / (binWidth * x.Stimuli.fully_presented_lines);
if (params.PST.nullThreshold > 0)
   pst_Y(find(pst_Y > params.PST.nullThreshold)) = NaN;
end
if (params.PST.spontSubtract ~= 0)
   [driv spont] = compute_driv_spont_rates(x, params);
   pst_Y = pst_Y - mean(spont);
end
pst_Y = pst_Y + multi.offset;
multi.offset = multi.offset + params.PST.pictSpacing;
params.multi.offset = multi.offset;

%%#######################################################################
function [params, ep_X, ep_Y] = calc_EP(x, params)
if (~isfield(params,'EP'))
    params.EP.scaleFactor = -0.0001;
    params.EP.pictSpacing = 0;
    params.EP.blockSpacing = 0;
    params.EP.nullThreshold = -1;
end
multi = params.multi;
multi.offset = multi.offset + params.EP.blockSpacing*multi.newBlock;
if (params.EP.nullThreshold > 0)
   ep_Y(find(ep_Y > params.EP.nullThreshold)) = NaN;
end
ep_Y = (x.EP.aveY * params.EP.scaleFactor) + multi.offset;
ep_X = [1:x.EP.lineLength] * x.EP.sampleInterval;
multi.offset = multi.offset + params.EP.pictSpacing;
params.multi.offset = multi.offset;


%%#######################################################################
function [params, cap_X, cap_Y] = calc_CAP(x, params)
if (~isfield(params,'CAP'))
    params.CAP.scaleFactor = 10;  % yields µV for a system gain of 10,000.
    params.CAP.levelsOffset = 10;
    params.CAP.XLabel = 'time (msec)';
    params.CAP.YLabel = 'attenuation (dB)';
end
% multi = params.multi;
%multi.offset = multi.offset + params.CAP.blockSpacing*multi.newBlock;
% if (params.CAP.nullThreshold > 0)
%    ep_Y(find(cap_Y > params.CAP.nullThreshold)) = NaN;
% end
for i = 1:size(x.Stimuli.RunLevels_params.attenMask, 2)
    cap_Y(:,i) = ((i-1)*params.CAP.levelsOffset) + params.CAP.scaleFactor * x.CAPData.CAPAvg_V{i}';
    params.CAP.axis.YTick(i) = ((i-1)*params.CAP.levelsOffset);
    params.CAP.axis.YTickLabel(i) = x.Stimuli.atten_dB + x.Stimuli.RunLevels_params.stepdB * x.Stimuli.RunLevels_params.attenMask(i);
end
cap_X = 1000*(x.Stimuli.RunLevels_params.decimateFact/x.Stimuli.RPsamprate_Hz) * [0:size(cap_Y, 1)-1]';
% params.CAP.axis.XLabel = 'time (msec)';
% params.CAP.axis.YLabel = 'attenuation (dB)';
%     cap_Y = params.CAP.scaleFactor * x.CAPData.CAPAvg_V{:}';
% % (x.CAP.aveY * params.CAP.scaleFactor) + multi.offset;
% cap_X = [1:x.CAP.lineLength] * x.CAP.sampleInterval;
% % multi.offset = multi.offset + params.CAP.pictSpacing;
% % params.multi.offset = multi.offset;

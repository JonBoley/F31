function quick_vowel(RLFpics,CALpic)
% File quick_EHrlfs.m
% M.Heinz: 11Nov2004 (From GE: quick_vowel)
% For NOHR
%
% USAGE:quick_vowel(RLFpics,CALpic)  [e.g., quick_vowel([2 3],1)]
% Plots RLFs for a set of vowel features
%
% RLFpics: vector of vowel RLFs for different features [F1_pic, T1_pic]
% CALpic: calib file (e.g., 1)

params = [];
params.spikeChan = 1;
params.figNum = 100;
params.TriFiltWidth=5;
params.colors={'b','r','g','k','m','c','y','b','r','g','k','m','c','y'};

Npicts=length(RLFpics);
rc=cell(1,Npicts);
x=cell(1,Npicts);
driv=cell(1,Npicts);
spont=cell(1,Npicts);
legtext=cell(1,Npicts);

figure(params.figNum); clf
if length(RLFpics)==1
   nameStr = sprintf ('picture p%04d rate-level', RLFpics);
else
   nameStr = sprintf ('pictures p%04d - p%04d rate-level', RLFpics(1), RLFpics(end));
end
set(gcf, 'Name', nameStr);

[dummy CAL]=loadPic(CALpic);
CAL.CalibData(:,2)=trifilt(CAL.CalibData(:,2)',5)';  % Do some smoothing

for PICind=1:Npicts
   [rc{PICind} x{PICind}] = loadPic(RLFpics(PICind));
   if isfield(x{PICind}.Stimuli,'feature')
      legtext{PICind}=x{PICind}.Stimuli.feature;   
   else
      legtext{PICind}='TN';   
   end
   [driv{PICind} spont{PICind}] = compute_driv_spont_rates(x{PICind}, params);
   
   if isfield(x{PICind}.Stimuli,'feature')
      if strcmp(x{PICind}.Stimuli.feature(1),'F')
         linetype=strcat(params.colors{PICind},'-');
      else
         linetype=strcat(params.colors{PICind},':');
      end
   else
      linetype=strcat(params.colors{PICind},'--');
   end
   
   if isfield(x{PICind}.Stimuli,'feature')
      max_dBSPL=CalibInterp(x{1}.Stimuli.featureTargetFreq_Hz/1000,CAL.CalibData(:,1:2));
   else
      max_dBSPL=CalibInterp(x{1}.Stimuli.main.freqs/1000,CAL.CalibData(:,1:2));
   end
   
   if isfield(x{1}.Stimuli,'BASELINE')
      levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2)+x{1}.Stimuli.BASELINE.dBreTONE;
   elseif isfield(x{PICind}.Line.attens,'list')
      dBreTONE=-6.26;  % Value for vow17_MH10k.wav
      levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2)+dBreTONE;
      beep
      disp(sprintf('dBreTONE not saved in data file, using %.2f dB (from ''vow17_MH10k.wav'' - CNexps)',dBreTONE))
   else
      dBreTONE=0;  % Assume TONE
      levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.Tone(1:x{PICind}.Stimuli.fully_presented_lines,2)+dBreTONE;
%       beep
%       disp(sprintf('Assuming TONE here '))
   end
   plot(levels_dBSPL, trifilt(driv{PICind},params.TriFiltWidth),linetype);
   hold on
   
   % plot fit
   x0 = [18 100 20 10 0.5];
   RLcurveFit = lsqcurvefit(@RLcurve,x0,levels_dBSPL',trifilt(driv{PICind},params.TriFiltWidth));
   plot(levels_dBSPL, RLcurve(RLcurveFit,levels_dBSPL), 'r:');
   TwoThird = interp1(RLcurve(RLcurveFit,levels_dBSPL),levels_dBSPL,RLcurveFit(1)+2/3*(RLcurveFit(2)-RLcurveFit(1)));
%    plot([1 1]*TwoThird,[RLcurveFit(1) RLcurveFit(2)],'r:');
   text('Interpreter','tex','Position',[TwoThird,RLcurve(RLcurveFit,TwoThird)],'String',['\leftarrow ' sprintf('%s= %.1fdB',x{PICind}.Stimuli.feature,TwoThird)]);
   fprintf('Pic%d: SR=%.2f/sec; Sat=%.2f/sec; Thesh10=%.2fdB; TwoThird=%.2fdB\n',...
       PICind,RLcurveFit(1),RLcurveFit(2),RLcurveFit(3),TwoThird);
end

% set(gca, 'XDir', 'reverse');
xlabel('Overall Vowel Level (dB SPL)')
ylabel('Driven rate (sp/sec)')
legend(legtext,2)
title(nameStr)

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

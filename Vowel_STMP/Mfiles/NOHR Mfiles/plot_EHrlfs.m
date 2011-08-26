function plot_vowel(RLFpics,CALpic)
% File plot_EHrlfs.m
% Modified from: quick_EHrlfs.m (to clean up for ISH)
% M.Heinz: 11Nov2004 (From GE: quick_vowel)
% For NOHR
%
% USAGE:quick_vowel(RLFpics,CALpic)  [e.g., quick_vowel([2 3],1)]
% Plots RLFs for a set of vowel features
%
% RLFpics: vector of vowel RLFs for different features [F1_pic, T1_pic]
% CALpic: calib file (e.g., 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO (4/27/06) _ ADD excludelines from PIC %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = [];
params.spikeChan = 1;
params.figNum = 100;
params.TriFiltWidth=5;
params.colors={'b','r','g','k','m','c','y','b','r','g','k','m','c','y'};

Npicts=length(RLFpics);
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
set(gcf, 'Name', nameStr,'pos',[271   539   560   420]);

CAL=loadPic(CALpic);
CAL.CalibData(:,2)=trifilt(CAL.CalibData(:,2)',5)';  % Do some smoothing

for PICind=1:Npicts
	x{PICind} = loadPic(RLFpics(PICind));
   if isfield(x{PICind}.Stimuli,'feature')
      legtext{PICind}=x{PICind}.Stimuli.feature;   
   elseif strcmp(deblank(x{PICind}.Stimuli.short_description),'SACrlvQ')
      legtext{PICind}='SACrlvQ';
   elseif strcmp(deblank(x{PICind}.Stimuli.short_description),'TB')
      legtext{PICind}='TN'; 
   else
      legtext{PICind}='???';
   end
   [driv{PICind} spont{PICind}] = compute_driv_spont_rates(x{PICind}, params);
   
   if isfield(x{PICind}.Stimuli,'feature')
      if strcmp(x{PICind}.Stimuli.feature(1),'F')
         linetype=strcat(params.colors{PICind},'-');
      else
         linetype=strcat(params.colors{PICind},'-.');
      end
   else
      linetype=strcat(params.colors{PICind},':');
   end
   
   if isfield(x{PICind}.Stimuli,'feature')
      max_dBSPL=CalibInterp(x{PICind}.Stimuli.featureTargetFreq_Hz/1000,CAL.CalibData(:,1:2));
   elseif strcmp(deblank(x{PICind}.Stimuli.short_description),'SACrlvQ')
      %% Try to find BF, before asking
      UnitFile=dir(sprintf('UNITSdata/unit.%d.%d.mat',x{1}.General.track,x{1}.General.unit));
      if ~isempty(UnitFile)
         eval(['load ' sprintf('UNITSdata/unit.%d.%d.mat',x{1}.General.track,x{1}.General.unit)])
         if ~isempty(unit.Info.BF_kHz)
            BF_kHz=unit.Info.BF_kHz;
         end
      end
      if ~exist('BF_kHz','var')
         DataFile=dir('DataList*');
         if ~isempty(DataFile)
            eval(['load ' DataFile.name])
            if ~isempty(DataList.Units{x{1}.General.track,x{1}.General.unit})
               BF_kHz=DataList.Units{x{1}.General.track,x{1}.General.unit}.Info.BF_kHz;
            end
         end
      end
      if ~exist('BF_kHz','var')
         BF_kHz=-1;
         while ~((BF_kHz>.04)&(BF_kHz<40))
            BF_kHz=input('Enter BF (kHz): ');
         end
      end
      max_dBSPL=CalibInterp(BF_kHz,CAL.CalibData(:,1:2));
   elseif strcmp(deblank(x{PICind}.Stimuli.short_description),'TB')
      max_dBSPL=CalibInterp(x{PICind}.Stimuli.main.freqs/1000,CAL.CalibData(:,1:2));
   else
      max_dBSPL=100;
      beep
      disp('MAX_dBSPL set to 100 dB SPL arbitrarily because no Calib Data');
   end
   
   if isfield(x{PICind}.Stimuli,'BASELINE')
      levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2)+x{PICind}.Stimuli.BASELINE.dBreTONE;
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
   
   if strcmp(deblank(x{PICind}.Stimuli.short_description),'SACrlvQ')
      plot(-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2), trifilt(driv{PICind},params.TriFiltWidth),linetype);
      xlabel(sprintf('Stimulus Attenuation (dB) [0 dB = %.f dB SPL]',max_dBSPL-0+x{PICind}.Stimuli.BASELINE.dBreTONE))
   else
      plot(levels_dBSPL, trifilt(driv{PICind},params.TriFiltWidth),linetype);
      xlabel('Stimulus Level (dB SPL)')
   end
   hold on
end

% set(gca, 'XDir', 'reverse');
% xlabel('Overall Vowel Level (dB SPL)')
ylabel('Driven rate (sp/sec)')
legend(legtext,2)
title(nameStr)

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

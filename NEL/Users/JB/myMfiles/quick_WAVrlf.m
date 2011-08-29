function quick_WAVrlf(RLFpic)
% File quick_EHrlfs.m
% M.Heinz: 11Nov2004 (From GE: quick_vowel)
% For NOHR
%
% USAGE:quick_WAVrlf(RLFpic)  
% Plots RLF vs ATTEN for driv and spont to get Threshold
%

params = [];
params.spikeChan = 1;
params.figNum = 100;
params.TriFiltWidth=5;
params.colors={'b','r','g','k','m','c','y','b','r','g','k','m','c','y'};

Npicts=length(RLFpic);
rc=cell(1,Npicts);
x=cell(1,Npicts);
driv=cell(1,Npicts);
spont=cell(1,Npicts);
legtext=cell(1,Npicts);

figure(params.figNum); clf

% CAL=loadpic(CALpic);
% CAL.CalibData(:,2)=trifilt(CAL.CalibData(:,2)',5)';  % Do some smoothing

PICind=1;

[rc{PICind} x{PICind}] = loadPic(RLFpic(PICind));
legtext{PICind}={'DRIV','SPONT'};   
[driv{PICind} spont{PICind}] = compute_driv_spont_rates(x{PICind}, params);
nameStr = sprintf ('picture p%04d rate-level [%s]', RLFpic,x{1}.Stimuli.short_description);
set(gcf, 'Name', nameStr);

%    max_dBSPL=CalibInterp(current_unit_bf,CAL.CalibData(:,1:2));

%    if isfield(x{1}.Stimuli,'BASELINE')
%       levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2)+x{1}.Stimuli.BASELINE.dBreTONE;
%    elseif isfield(x{PICind}.Line.attens,'list')
%       dBreTONE=-6.26;  % Value for vow17_MH10k.wav
%       levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2)+dBreTONE;
%       beep
%       disp(sprintf('dBreTONE not saved in data file, using %.2f dB (from ''vow17_MH10k.wav'' - CNexps)',dBreTONE))
%    else
%       dBreTONE=0;  % Assume TONE
%       levels_dBSPL=max_dBSPL-x{PICind}.Line.attens.Tone(1:x{PICind}.Stimuli.fully_presented_lines,2)+dBreTONE;
% %       beep
% %       disp(sprintf('Assuming TONE here '))
%    end
ATTENS_dB=x{PICind}.Line.attens.list(1:x{PICind}.Stimuli.fully_presented_lines,2);
plot(ATTENS_dB, trifilt(driv{PICind},params.TriFiltWidth),'r-');
hold on
plot(ATTENS_dB, trifilt(spont{PICind},params.TriFiltWidth),'g-');

% set(gca, 'XDir', 'reverse');
xlabel('Attenuation (dB)')
ylabel('Driven rate (sp/sec)')
legend('DRIV','SPONT',2)
title(nameStr)
set(gca,'Xdir','rev')

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

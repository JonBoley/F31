function quick_vowel(picNum)

params = [];
params.spikeChan = 1;
params.figNum = 100;
params.TriFiltWidth=5;
params.colors={'b','r','g','k','m','c','y','b','r','g','k','m','c','y'};

Npicts=length(picNum);
rc=cell(1,Npicts);
x=cell(1,Npicts);
driv=cell(1,Npicts);
spont=cell(1,Npicts);
legtext=cell(1,Npicts);

figure(params.figNum); clf
if length(picNum)==1
   nameStr = sprintf ('picture p%04d rate-level', picNum);
else
   nameStr = sprintf ('pictures p%04d - p%04d rate-level', picNum(1), picNum(end));
end
set(gcf, 'Name', nameStr);

for PICind=1:Npicts
   [rc{PICind} x{PICind}] = loadPic(picNum(PICind));
   legtext{PICind}=x{PICind}.Stimuli.feature;   
   [driv{PICind} spont{PICind}] = compute_driv_spont_rates(x{PICind}, params);

   if strcmp(x{PICind}.Stimuli.feature(1),'F')
      linetype=strcat(params.colors{PICind},'-');
   else
      linetype=strcat(params.colors{PICind},':');
   end
   
   plot(x{PICind}.Line.attens.list(:,2), trifilt(driv{PICind},params.TriFiltWidth),linetype);
   hold on
end

set(gca, 'XDir', 'reverse');
xlabel('Atten (dB)')
legend(legtext,2)


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

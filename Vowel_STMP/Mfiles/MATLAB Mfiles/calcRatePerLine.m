function PIC=calcRatePerLine(PIC)
% M.Heinz 13Sep2004.  Taken from PSTview.m
% Calculates Rate (driv and spont) per line from PIC.
% Uses windows: driv=[10,dur+10], spont=[dur+200,period]
%
% Usage:PIC=calcRatePerLine(PIC)
% Input PIC: has PIC.x with stimulus parameters and spikes
% Output PIC: Stores calcs in PIC.RatePerLine, with paramaters used


% compute driven and spontaneous rates for each line
spikeTimes = PIC.x.spikes{1};

driven_window_sec=[10 PIC.x.Hardware.Trigger.StmOn+10]/1000;

if PIC.x.Hardware.Trigger.StmOn+PIC.x.Hardware.Trigger.StmOff-200 >250
   OFFSET=200;
else
   OFFSET=50;
end

spont_window_sec=[PIC.x.Hardware.Trigger.StmOn+OFFSET PIC.x.Hardware.Trigger.StmOn+PIC.x.Hardware.Trigger.StmOff]/1000;

for line = 1:PIC.x.Stimuli.fully_presented_lines
   spikeIndices = find((spikeTimes(:,1) == line ));
   driv(line) = length( find( (spikeTimes(spikeIndices,2) >= driven_window_sec(1)) & ...
      (spikeTimes(spikeIndices,2) < driven_window_sec(2)) ) );
   spont(line) = length( find( (spikeTimes(spikeIndices,2) >= spont_window_sec(1)) & ...
      (spikeTimes(spikeIndices,2) < spont_window_sec(2)) ) );
end
driv = driv / diff(driven_window_sec);
spont = spont / diff(spont_window_sec);

% Store calcs
PIC.RatePerLine.driv=driv;
PIC.RatePerLine.spont=spont;

% Store parameters used for calcs
PIC.RatePerLine.params.driven_window_sec=driven_window_sec;
PIC.RatePerLine.params.spont_window_sec=spont_window_sec;

return;

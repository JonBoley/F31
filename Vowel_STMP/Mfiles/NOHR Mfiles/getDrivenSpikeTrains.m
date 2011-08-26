function SpikeTrains=getDrivenSpikeTrains(spikes,excludeLines,Twin_sec)
% M.Heinz 27Sep2004.  
% Organizes SpikeTrain reps into cell array.  Takes spikes vector and returns a cell array (Nreps x 1) 
% with SpikeTrains of spikes within Twin_sec.  Allows lines to be excluded.
%
% Usage: SpikeTrains=getDrivenSpikeTrains(spikes,excludeLines,Twin_sec)
% spikes: vector of [line #, spike time]: same as NELdata format
% excludeLines: vector of lines to exclude
% Twin_sec: [StartTime_sec, EndTime_sec]
% SpikeTrains: Cell array {# reps x 1} of spike times within Twin_sec

if ~exist('excludeLines','var')
   excludeLines=[];
end
if ~exist('Twin_sec','var')
   Twin_sec=[-Inf Inf];  % Keep all spikes
end

GOODreps=setdiff(1:max(spikes(:,1)),excludeLines);
SpikeTrains=cell(1,length(GOODreps));
goodREPind=0;
for REPind=GOODreps
   goodREPind=goodREPind+1;
   SpikeTrains{goodREPind}=spikes(find((spikes(:,1)==REPind)&(spikes(:,2)>=Twin_sec(1))&(spikes(:,2)<=Twin_sec(2))),2)';
end


return;



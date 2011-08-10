function [SpikeTrain_new,Nspikes]=windowSTs(SpikeTrain_old,start_sec,end_sec,MAXspikes)
% File: [SpikeTrain_new,Nspikes]=windowSTs(SpikeTrain_old,start_sec,end_sec,MAXspikes)
%
% WINDOWS SpikeTrains given starting and ending time.
% Assumes spikes are in CCC-spiketrain format [cell_array{Nreps}]

Nspikes=0;
for  i = 1:length(SpikeTrain_old)
	spikeINDs=find((SpikeTrain_old{i}>=start_sec) & (SpikeTrain_old{i}<=end_sec));
	if (Nspikes+length(spikeINDs))<=MAXspikes
		Nspikes=Nspikes+length(spikeINDs);
		SpikeTrain_new{i} = SpikeTrain_old{i}(spikeINDs);
	else 
		break
	end
end

return;


function [STcell,Nspikes]=nelSTs2cell(NELSpikeTrains)
% File: [STcell,Nspikes]=nelSTs2cell(NELSpikeTrains)
%
% CONVERT SpikeTrains in NEL column format: [rep #,spiketime_sec] to
% CCC-spiketrain format [cell_array{Nreps}]

Nrep=max(NELSpikeTrains(:,1));
STcell=cell(1,Nrep);
Nspikes=size(NELSpikeTrains,1);

for i=1:Nrep
	repINDs=find(NELSpikeTrains(:,1)==i);
	STcell{i}=NELSpikeTrains(repINDs,2);
end

return;


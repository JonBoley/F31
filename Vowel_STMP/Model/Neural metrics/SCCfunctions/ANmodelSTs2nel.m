function [NELspikes,nspikes]=ANmodelSTs2nel(sptimes,Nreps)
% File: [NELSpikeTrains,Nspikes]=ANmodelSTs2nel(sptimes)
%
% CONVERT SpikeTrains in ANmodel row format: e.g., ARLO, ZB06/07 to 
% NEL-spiketrain format [rep #,spiketime_sec]

nspikes=length(sptimes);

NELspikes=NaN*ones(nspikes,2);
NELspikes(:,2)=sptimes;
REPendINDs(1)=0;
Nreps = 1+length(find(diff(sptimes)<0)); % needed for low amplitude signals
REPendINDs(2:Nreps)=find(diff(sptimes)<0);
REPendINDs(Nreps+1)=nspikes;
for REPind=1:Nreps
	NELspikes(REPendINDs(REPind)+1:REPendINDs(REPind+1),1)=REPind;
end


return;


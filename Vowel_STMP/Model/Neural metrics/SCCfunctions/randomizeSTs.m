function ST_rand=randomizeSTs(ST_orig,start_sec,end_sec);
% File: ST_rand=randomizeSTs(ST_orig,start_sec,end_sec);
%
% RANDOMIZES SpikeTrains given starting and ending time.  Maintains same
% number of spikes in each REP, but with uniform distribution between start
% and end times. Assumes spikes are in CCC-spiketrain format:
% [cell_array{Nreps}] 

ST_rand=cell(size(ST_orig));
for i=1:length(ST_orig)
	if max(ST_orig{i}>end_sec)|min(ST_orig{i}<start_sec)
		error('in randomizeSTs: spikes exist outside of specified time window in ST_orig!')
	end
	ST_rand{i}= sort(start_sec + (end_sec-start_sec)*rand(length(ST_orig{i}),1));
end

return;


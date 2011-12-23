% get_spikes.m
% Normal hearing
vihc = catmodel_IHC(signal_model.',CF_kHz(FiberNumber)*1e3,1,...
    1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
[sout,psth]=catmodel_Synapse(vihc,CF_kHz(FiberNumber)*1e3,1,...
    1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
[sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
SpikeTrains_plus = nelSTs2cell(NELspikes);
if length(SpikeTrains_plus)<Nreps
    SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
end
Spikes_plus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_plus;


vihc = catmodel_IHC(-signal_model.',CF_kHz(FiberNumber)*1e3,1,...
    1/ANmodel_Fs_Hz,dur_sec+1.00,Cohc,Cihc);
[sout,psth]=catmodel_Synapse(vihc,CF_kHz(FiberNumber)*1e3,1,...
    1/ANmodel_Fs_Hz,dur_sec+1.00,fibertype,1);
SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
[sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
SpikeTrains_minus = nelSTs2cell(NELspikes);
if length(SpikeTrains_minus)<Nreps
    SpikeTrains_minus = [SpikeTrains_minus cell(1,Nreps-length(SpikeTrains_minus))];
end
Spikes_minus{FiberNumber,LevelIndex,SNRindex} = SpikeTrains_minus;

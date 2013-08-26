% Initialize analysis parameters
SynOut=cell(numCFs,length(Levels),length(SNRs));
Spikes_plus=cell(numCFs,length(Levels),length(SNRs));
Spikes_minus=cell(numCFs,length(Levels),length(SNRs));
Rho=NaN*ones(numCFs,length(Levels),length(SNRs));
CD=NaN*ones(numCFs,length(Levels),length(SNRs));
Rate=NaN*ones(numCFs,length(Levels),length(SNRs));
Rho_vDeltaCF=NaN*ones(length(Levels),length(SNRs),numCFs);
CD_vDeltaCF=NaN*ones(length(Levels),length(SNRs),numCFs);

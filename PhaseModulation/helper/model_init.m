% initialize model parameters

numCFs=40;
if mod(numCFs,2)==1, numCFs=numCFs+1; end % make this an even number

spread=0.5; %number of octaves
lowCF_kHz=midCF_kHz*2^(-spread/2);
CF_kHz=(lowCF_kHz*2.^(0:spread/(numCFs-1):spread));

% Make the first fiber the one centered at the feature
CF_kHz = [midCF_kHz CF_kHz];
numCFs = length(CF_kHz);
midCF_index = floor(numCFs/2)+1;
deltaCF = log2([CF_kHz(2:midCF_index) CF_kHz(1) CF_kHz(midCF_index+1:end)]'/midCF_kHz);

fibertype=3; %[1,2,3]=[low,med,high] spontaneous rate
Cohc=1.0; %outer hair cell health
Cihc=1.0; %inner hair cell health
Nreps=120; %number of repetitions
ANmodel_Fs_Hz=100e3; %100kHz sampling

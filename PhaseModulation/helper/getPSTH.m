function PSTH = getPSTH(SpikeTrains_plus,binSize_sec)

if nargin<2
    binSize_sec = 20e-6;
end

maxTime = max(cell2mat(SpikeTrains_plus'));
PSTH = zeros(ceil(maxTime/binSize_sec),1);
for repIndex = 1:numel(SpikeTrains_plus)
    repHist = hist(SpikeTrains_plus{repIndex},...
        binSize_sec/2:binSize_sec:(maxTime+(binSize_sec/2)));
    PSTH = PSTH + repHist';
end


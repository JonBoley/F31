function NeuralDelay_sec = estimateClickLatency(spikeTimes)
% NeuralDelay_sec = estimateClickLatency(spikeTimes)
% 
% from Scheidt et al 2010:
% firstPeakDelay_sec is first peak in PSTH with discharge rate >= 70% of
% max
% drivenResponseDelay_sec is second peak > (spont rate + 2 standard
% deviations)



PSTHbinWidth_sec = 100e-6;
nreps = length(unique(spikeTimes(:,1)));
    
spontDur = 0.050; %sec

startIndices = 1+[0; find(max(0,diff(spikeTimes(:,1)==1)))];
endIndices = [find(max(0,diff(spikeTimes(:,1)==1))); length(spikeTimes(:,1))];
mreps = length(startIndices); %may have reps of reps

clear localPSTH;
for i=1:mreps
    localPSTH(:,i) = histc(spikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(spikeTimes(:,2)));
    localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
    
    psthPeakThresh(i) = 0.7 * max(localPSTH(:,i));
    
    psthDur = max(spikeTimes(startIndices(i):endIndices(i),2));
    spontSpikes = localPSTH(max(1,end-round(spontDur/PSTHbinWidth_sec)):end,i); % last N spikes
    
    psthDrivenThresh(i) = max(spontSpikes) + 2*std(spontSpikes);
end


PSTHnum = 1:mreps;
if length(PSTHnum)<2
    firstPeakDelay_sec = PSTHbinWidth_sec*...
        find(localPSTH(:,PSTHnum)>psthPeakThresh(PSTHnum),1,'first');

    tempIndex = find(localPSTH(:,PSTHnum)>psthDrivenThresh(PSTHnum),2,'first');
    drivenResponseDelay_sec = PSTHbinWidth_sec*max(tempIndex);
    
else %otherwise combine PSTHs somehow
    
    % calculate delay for each PSTH
    for i=PSTHnum
        localPSTH2 = localPSTH(:,i); %pick one PSTH

        psthPeakThresh = 0.7 * max(localPSTH2);
        psthDrivenThresh2 = psthDrivenThresh(i);

        firstPeakDelay_sec(i) = PSTHbinWidth_sec*...
            find(localPSTH2>psthPeakThresh,1,'first');

        tempIndex = find(localPSTH2>psthDrivenThresh2,2,'first');
        if isempty(tempIndex)
            drivenResponseDelay_sec(i) = NaN;
        else
            drivenResponseDelay_sec(i) = PSTHbinWidth_sec*max(tempIndex);
        end
    end
    

end


NeuralDelay_sec = drivenResponseDelay_sec;


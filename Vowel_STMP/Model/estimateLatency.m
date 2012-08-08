function NeuralDelay_sec = estimateLatency(spikeTimes)
% NeuralDelay_sec = estimateLatency(spikeTimes)
% 
% from Scheidt et al 2010:
% firstPeakDelay_sec is first peak in PSTH with discharge rate >= 70% of
% max
% drivenResponseDelay_sec is second peak > (spont rate + 2 standard
% deviations)
plotPSTH=0;

PSTHbinWidth_sec = 100e-6;
spontDur = 0.001; %sec (use first 1ms... ending times will all be different for STMP)
nreps = length(unique(spikeTimes(:,1)));

startIndices = 1+[0; find(max(0,diff(spikeTimes(:,1)==1)))];
endIndices = [find(max(0,diff(spikeTimes(:,1)==1))); length(spikeTimes(:,1))];
mreps = length(startIndices); %may have reps of reps

clear localPSTH;
for i=1:mreps
    localPSTH(:,i) = histc(spikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(spikeTimes(:,2)));
    localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
    
    psthPeakThresh(i) = 0.7 * max(localPSTH(:,i));
    
    spontSpikes = localPSTH(1:round(spontDur/PSTHbinWidth_sec),i);
    psthDrivenThresh(i) = max(spontSpikes) + 2*std(spontSpikes);
end


PSTHnum = 1:mreps;
if length(PSTHnum)<2
    firstPeakDelay_sec = PSTHbinWidth_sec*...
        find(localPSTH(:,PSTHnum)>psthPeakThresh(PSTHnum),1,'first');

    tempIndex = find(localPSTH(:,PSTHnum)>psthDrivenThresh(PSTHnum),2,'first');
    drivenResponseDelay_sec = PSTHbinWidth_sec*max(tempIndex);
    
    if plotPSTH
        bar(0:PSTHbinWidth_sec:max(spikeTimes(:,2)),localPSTH(:,PSTHnum)); hold on;
        plot([0; max(spikeTimes(:,2))],[psthPeakThresh(PSTHnum); psthPeakThresh(PSTHnum)],'b');
        plot([0; max(spikeTimes(:,2))],[psthDrivenThresh(PSTHnum); psthDrivenThresh(PSTHnum)],'g'); hold off;
    end
    
else %otherwise pool
    localPSTH = sum(localPSTH(:,PSTHnum),2);
    
    psthPeakThresh = 0.7 * max(localPSTH);
    
    spontSpikes = localPSTH(1:round(spontDur/PSTHbinWidth_sec));
    psthDrivenThresh = max(spontSpikes) + 2*std(spontSpikes);
    
    firstPeakDelay_sec = PSTHbinWidth_sec*...
        find(localPSTH>psthPeakThresh,1,'first');

    tempIndex = find(localPSTH>psthDrivenThresh,2,'first');
    drivenResponseDelay_sec = PSTHbinWidth_sec*max(tempIndex);
    
    if plotPSTH
        bar(0:PSTHbinWidth_sec:max(spikeTimes(:,2)),localPSTH); hold on;
        plot([0; max(spikeTimes(:,2))],[psthPeakThresh; psthPeakThresh],'b');
        plot([0; max(spikeTimes(:,2))],[psthDrivenThresh; psthDrivenThresh],'g'); hold off;
    end
end

NeuralDelay_sec = drivenResponseDelay_sec;


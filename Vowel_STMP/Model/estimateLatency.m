function NeuralDelay_sec = estimateLatency(spikeTimes,poolPSTH,plotPSTH)
% NeuralDelay_sec = estimateLatency(spikeTimes)
% 
% from Scheidt et al 2010:
% firstPeakDelay_sec is first peak in PSTH with discharge rate >= 70% of
% max
% drivenResponseDelay_sec is second peak > (spont rate + 2 standard
% deviations)

if nargin<2
    poolPSTH = 0;
end
if nargin<3
    plotPSTH=0;
end

PSTHbinWidth_sec = 100e-6;
nreps = length(unique(spikeTimes(:,1)));
    
if poolPSTH
    spontDur = 0.001; %sec
else
    spontDur = 0.200; %sec
end

startIndices = 1+[0; find(max(0,diff(spikeTimes(:,1)==1)))];
endIndices = [find(max(0,diff(spikeTimes(:,1)==1))); length(spikeTimes(:,1))];
mreps = length(startIndices); %may have reps of reps

clear localPSTH;
for i=1:mreps
    localPSTH(:,i) = histc(spikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(spikeTimes(:,2)));
    localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
    
    psthPeakThresh(i) = 0.7 * max(localPSTH(:,i));
    
    psthDur = max(spikeTimes(startIndices(i):endIndices(i),2));
    if poolPSTH
        spontSpikes = localPSTH(1:round(spontDur/PSTHbinWidth_sec),i); % first N spikes
    else
        spontSpikes = localPSTH(round((psthDur-spontDur)/PSTHbinWidth_sec):...
            round(psthDur/PSTHbinWidth_sec),i); % last N spikes
    end
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
    
else %otherwise combine PSTHs somehow
    
    if poolPSTH
        %pool all PSTHs
        localPSTH = sum(localPSTH(:,PSTHnum),2); %pool all PSTHs
        psthPeakThresh = 0.7 * max(localPSTH);
        
        spontSpikes = localPSTH(1:round(spontDur/PSTHbinWidth_sec));
        psthDrivenThresh = max(spontSpikes) + 2*std(spontSpikes);
        tempIndex = find(localPSTH>psthDrivenThresh,2,'first');
        drivenResponseDelay_sec = PSTHbinWidth_sec*min(tempIndex);
    else
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
    
    if plotPSTH
        bar(0:PSTHbinWidth_sec:max(spikeTimes(:,2)),localPSTH); hold on;
        plot([0; max(spikeTimes(:,2))],[psthPeakThresh; psthPeakThresh]);
        plot([0; max(spikeTimes(:,2))],[psthDrivenThresh; psthDrivenThresh]); hold off;
    end
end

NeuralDelay_sec = drivenResponseDelay_sec;


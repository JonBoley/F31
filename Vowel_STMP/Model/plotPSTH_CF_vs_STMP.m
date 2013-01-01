%% plotPSTH_CF_vs_STMP

filenames = {'STMPvsCF_CF_*NH.mat','STMPvsCF_STMP_*NH.mat','STMPvsCF_STMP_2012-12-30_131541.mat'};
colors = 'kgr';
legendText = {'CF','STMP','STMP (corrected)'};

for z=1:length(filenames)
    clear -regexp ^((?!\<filenames\>|\<z\>|\<h\>|\<colors\>|\<legendText\>).)*$ %keep CFfiles, z, colors
    home;

    CFfiles=dir(filenames{z});
    CFfilename = CFfiles.name;
    load(CFfilename,'-regexp','[^CFfilename^STMPfilename^CFfiles^filenames^z^h^colors^legendText]');
    LocalSpikeTimes = [];
    for FiberNumber=1:numCFs %pool across simulated CF
        for i=1:length(Spikes_plus{FiberNumber,LevelIndex,SNRindex})
            LocalSpikeTimes = [LocalSpikeTimes;...
                repmat(i,length(Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}),1), ...
                Spikes_plus{FiberNumber,LevelIndex,SNRindex}{i}];
        end
    end
    NeuralDelay_sec = estimateLatency(LocalSpikeTimes,0,0);

    % build PSTH
    PSTHbinWidth_sec = 100e-6;
    nreps = length(unique(LocalSpikeTimes(:,1)));

    startIndices = 1+[0; find(max(0,diff(LocalSpikeTimes(:,1)==1)))];
    endIndices = [find(max(0,diff(LocalSpikeTimes(:,1)==1))); length(LocalSpikeTimes(:,1))];
    mreps = length(startIndices); %may have reps of reps

    clear localPSTH;
    for i=1:mreps
        localPSTH(:,i) = histc(LocalSpikeTimes(startIndices(i):endIndices(i),2),0:PSTHbinWidth_sec:max(LocalSpikeTimes(:,2)));
        localPSTH(:,i) = (localPSTH(:,i)/nreps) / PSTHbinWidth_sec; % spikes per sec
    end

    sizeFactor = 50/max(max(localPSTH.^2));
    thresh = 500^2;
    figure(1), hold on;
    for i=1:size(localPSTH,2)
        localData = localPSTH(1:300,i).^2;
        localData(localData<=thresh)=NaN;
        h(z)=scatter((1:size(localData,1))*PSTHbinWidth_sec,CF_kHz(i)*ones(1,size(localData,1)),sizeFactor*localData,colors(z));
        drawnow;
    end


end

legend(h,legendText);
ylabel('CF (kHz)'); xlabel('Time (sec)');
set(gca,'YScale','log','YTick',[0.25 0.5 1 1.7 2 2.5 4],...
    'YTickLabel',{'0.25','0.5','1.0','1.7','2.0','2.5','4.0'});


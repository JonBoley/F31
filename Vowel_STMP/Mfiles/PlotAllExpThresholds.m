%% plot all AN thresholds
ExpDataDir = 'C:\Research\MATLAB\Vowel_STMP\ExpData\';
cd(ExpDataDir);

% search for all AN experiments
Experiments = dir(fullfile(ExpDataDir,'JB-2011*AN*'));

figure(1); clf;
load normema;

%%%%%%% Add lines from Kale and Heinz (2009) normal-hearing CHIN data
QlowM97(:,1)=[.4 10]';  QhighM97(:,1)=QlowM97(:,1);
% 5th and 95th percentiles from Kale and Heinz (2009)
QlowM97(:,2)=10^.1572*QlowM97(:,1).^.1416;
QhighM97(:,2)=10^.7989*QlowM97(:,1).^.1416;

% go through each experiment
Thrs_dbSPL = cell(length(Experiments),1);
BFs_kHz = cell(length(Experiments),1);
Q10s = cell(length(Experiments),1);
for ExpNum=1:length(Experiments)
    cd(fullfile(ExpDataDir,Experiments(ExpNum).name));
    
    CALIBpics=findPics('calib');
    TrackUnitList=getTrackUnitList;
    NUMunits=size(TrackUnitList,1);
    Thrs_dbSPL{ExpNum} = -999+zeros(NUMunits,1);
    BFs_kHz{ExpNum} = -999+zeros(NUMunits,1);
    Q10s{ExpNum} = -999+zeros(NUMunits,1);
    for ind = 1:NUMunits
        % Plot TC
        TCpics=findPics('tc',TrackUnitList(ind,:));
        if ~isempty(TCpics)
            Thrs_dbSPL_TEMP=[];
            BFs_kHz_TEMP=[];
            Q10s_TEMP=[];
            for TCind=1:length(TCpics)
                CALIB_PIC = CALIBpics(max(find((CALIBpics<TCpics(TCind)))));
                [Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs(TCpics(TCind),CALIB_PIC,0);
                % [Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs_BFbyhand(TCpics(TCind),CALIB_PIC,0);
            end
            [Thrs_dbSPL{ExpNum}(ind) index_TC] = min(Thrs_dbSPL_TEMP);
            BFs_kHz{ExpNum}(ind) = BFs_kHz_TEMP(index_TC);
            Q10s{ExpNum}(ind) = Q10s_TEMP(index_TC);
        end %if ~isempty(TCpics)
    end %for ind = 1:NUMunits
    
    if strfind(lower(Experiments(ExpNum).name),'normal')
        figure(1), subplot(211), hold on;
        semilogx(BFs_kHz{ExpNum},Thrs_dbSPL{ExpNum},'kx'); hold off;
        subplot(212), hold on;
        loglog(BFs_kHz{ExpNum},Q10s{ExpNum},'kx'); hold off;
    elseif (~isempty(strfind(lower(Experiments(ExpNum).name),'impaired')) || ~isempty(strfind(lower(Experiments(ExpNum).name),'obn')))
        figure(1), subplot(211), hold on;
        semilogx(BFs_kHz{ExpNum},Thrs_dbSPL{ExpNum},'rx'); hold off;
        subplot(212), hold on;
        loglog(BFs_kHz{ExpNum},Q10s{ExpNum},'rx'); hold off;
    end
    
end %for ExpNum=1:length(Experiments)

xmin=0.03; xmax=10; ymin=-20; ymax=110;
figure(1), subplot(211), hold on;
semilogx(normt(1,:),normt(2,:),'k');
hold off;
ylabel('dB SPL'); xlabel('Frequency (kHz)');
axis([xmin xmax ymin ymax]);
set(gca,'YTick',[0:10:100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
title('All Experiments');
grid on

subplot(212), hold on;
loglog(QlowM97(:,1),QlowM97(:,2),'k-');
loglog(QhighM97(:,1),QhighM97(:,2),'k-');
hold off;
ylabel('Q10'); xlabel('Frequency (kHz)');
xlim([xmin xmax]);
ylim([.1 10]);
% set(gca,'YTick',[0 20 40 60 80 100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
grid on

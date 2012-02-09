%% plot all AN thresholds
switch lower(getenv('computername'))
    case {'m4500'}
        ExpDataDir = 'C:\Research\MATLAB\Vowel_STMP\ExpData';
    otherwise
        ExpDataDir = 'C:\NEL\ExpData\';
end
cd(ExpDataDir);

% search for all AN experiments
% Experiments = dir(fullfile(ExpDataDir,'JB-2011*AN*'));
Experiments = dir(fullfile(ExpDataDir,'*AN*'));

%%% experiment dates by exposure type
% 500Hz OBN 117dB SPL (2 hrs)
SNHL_117 = {...
    'SK-2010_11_11-AN_impaired',...
    'SK-2011_02_08-AN_impaired',...
    'SK-2011_03_08-AN_impaired',...% both chins 10-27 & 10-32 ???
    'SK-2011_03_28-AN_impaired',...
    'SK-2011_05_03-AN_impaired',...
    'JB-2011_08_29-Chin1103_AN_500OBN'};
% 500Hz OBN 116dB SPL (2 hrs)
SNHL_116 = {...
    'JB-2011_05_16-AN_impaired',...
    'KH-2011_06_22-chin1119_AN_500OBN',...
    'JB-2011_06_23-Chin1119_AN_500OBN',...
    'KH-2011_07_06-chin1123_AN_500OBN',...
    'JB-2011_07_21-Chin1124_AN_500OBN',...
    'KH-2011_07_22-chin1124_AN_500OBN',...
    'JB-2011_08_01-Chin1125_AN_500OBN',...
    'KH-2011_08_02-chin1125_AN_500OBN',...
    'JB-2011_08_09-Chin1135_AN_500OBN',...
    'KH-2011_08_10-chin1135_AN_500OBN',...
    'JB-2011_08_15-Chin1136_AN_500OBN',...
    'KH-2011_08_16-chin1136_AN_500OBN',...
    'JB-2011_09_07-Chin1131_AN_500OBN',...
    'JB-2011_09_19-Chin1134_AN_500OBN'};

figure(1); clf;
load normema;

ThirdOctaveskHz = 0.2*2.^(0:1/3:6);
MinThresh = NaN*ones(length(ThirdOctaveskHz),1);
AvgThresh = NaN*ones(length(ThirdOctaveskHz),1);
MinThresh_Normal = NaN*ones(length(ThirdOctaveskHz),1);
AvgThresh_Normal = NaN*ones(length(ThirdOctaveskHz),1);

AllThresh_Normal = []; % [Thrs_dbSPL BFs_kHz]
AllThresh_Impaired = [];

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
    if isdir(fullfile(ExpDataDir,Experiments(ExpNum).name))
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
                    try
                        [Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs(TCpics(TCind),CALIB_PIC,0);
                        % [Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs_BFbyhand(TCpics(TCind),CALIB_PIC,0);
                    catch
                        Thrs_dbSPL_TEMP(TCind)=NaN;
                        BFs_kHz_TEMP(TCind)=NaN;
                        Q10s_TEMP(TCind)=NaN;
                    end
                end
                [Thrs_dbSPL{ExpNum}(ind) index_TC] = min(Thrs_dbSPL_TEMP);
                BFs_kHz{ExpNum}(ind) = BFs_kHz_TEMP(index_TC);
                Q10s{ExpNum}(ind) = Q10s_TEMP(index_TC);
            end %if ~isempty(TCpics)
        end %for ind = 1:NUMunits

        % Normal
        if strfind(lower(Experiments(ExpNum).name),'norm')
            AllThresh_Normal = [AllThresh_Normal; Thrs_dbSPL{ExpNum} BFs_kHz{ExpNum}];
            
            % Find avg normal threshold
            for i=2:length(ThirdOctaveskHz)
                localThreshes = AllThresh_Normal((AllThresh_Normal(:,2)>=ThirdOctaveskHz(i-1) & AllThresh_Normal(:,2)<ThirdOctaveskHz(i)),1);
                if ~isempty(localThreshes), AvgThresh_Normal(i-1)=mean(localThreshes); end
                if ~isempty(localThreshes), MinThresh_Normal(i-1)=min(localThreshes); end
            end
            
            figure(1), subplot(211), hold on;
            semilogx(BFs_kHz{ExpNum},Thrs_dbSPL{ExpNum},'k.','MarkerSize',3); hold off;
            xmin=0.03; xmax=10; ymin=-20; ymax=110;
            axis([xmin xmax ymin ymax]);
            subplot(212), hold on;
            loglog(BFs_kHz{ExpNum},Q10s{ExpNum},'k.','MarkerSize',3); hold off;
            xlim([xmin xmax]);
            ylim([.1 10]);
            
            if exist('hThresh3','var')
                try
                    set(hThresh3,'YDataSource','AvgThresh_Normal');
                    set(hThresh4,'YDataSource','MinThresh_Normal');
                    refreshdata;
                catch % data linking doesn't work for older versions of Matlab
                    set(hThresh3,'YData',AvgThresh_Normal);
                    set(hThresh4,'YData',MinThresh_Normal);
                end
            else
                figure(1), subplot(211), hold on;
                hThresh3 = semilogx(ThirdOctaveskHz*2^(1/3),AvgThresh_Normal,'k:','Linewidth',2);
                hThresh4 = semilogx(ThirdOctaveskHz*2^(1/3),MinThresh_Normal,'k','Linewidth',2); hold off;
            end
            
        % Impaired
        elseif (~isempty(strfind(lower(Experiments(ExpNum).name),'impaired')) | ~isempty(strfind(lower(Experiments(ExpNum).name),'obn')))
            if strmatch(Experiments(ExpNum).name,SNHL_116) % if exposure was at 116dB SPL
                marker = 'r.';
            elseif strmatch(Experiments(ExpNum).name,SNHL_117) % if exposure was at 117dB SPL
                marker = 'm.';
            else
                marker = '';
            end
            
            if ~isempty(marker)
                AllThresh_Impaired = [AllThresh_Impaired; Thrs_dbSPL{ExpNum} BFs_kHz{ExpNum}];
                
                figure(1), subplot(211), hold on;
                semilogx(BFs_kHz{ExpNum},Thrs_dbSPL{ExpNum},marker); hold off;
                xmin=0.03; xmax=10; ymin=-20; ymax=110;
                axis([xmin xmax ymin ymax]);
                subplot(212), hold on;
                loglog(BFs_kHz{ExpNum},Q10s{ExpNum},marker); hold off;
                xlim([xmin xmax]);
                ylim([.1 10]);

                % Find threshold shift
                for i=2:length(ThirdOctaveskHz)
                    %     BFs_kHz{ExpNum},Thrs_dbSPL{ExpNum}
%                     localThreshes = Thrs_dbSPL{ExpNum}(BFs_kHz{ExpNum}>=ThirdOctaveskHz(i-1) & BFs_kHz{ExpNum}<ThirdOctaveskHz(i));
                    localThreshes = AllThresh_Impaired((AllThresh_Impaired(:,2)>=ThirdOctaveskHz(i-1) & AllThresh_Impaired(:,2)<ThirdOctaveskHz(i)),1);
                    if ~isempty(localThreshes), AvgThresh(i-1)=mean(localThreshes); end
                    if ~isempty(localThreshes), MinThresh(i-1)=min(localThreshes); end
                end
                if exist('hThresh','var')
                    try
                        set(hThresh,'YDataSource','AvgThresh');
                        set(hThresh2,'YDataSource','MinThresh');
                        refreshdata;
                    catch % data linking doesn't work for older versions of Matlab
                        set(hThresh,'YData',AvgThresh);
                        set(hThresh2,'YData',MinThresh);
                    end
                    try
                        set(hThresh,'YDataSource','AvgThresh');
                        refreshdata;
                    catch % data linking doesn't work for older versions of Matlab
                        set(hThresh,'YData',AvgThresh);
                    end
                else
                    figure(1), subplot(211), hold on;
                    hThresh = semilogx(ThirdOctaveskHz*2^(1/3),AvgThresh,'r:','Linewidth',2);
                    hThresh2 = semilogx(ThirdOctaveskHz*2^(1/3),MinThresh,'r','Linewidth',2); hold off;
                end
            end %if ~isempty(marker)
        end % if normal, else impaired

        pause(0.01);
    end %if isdir(...)
end %for ExpNum=1:length(Experiments)

xmin=0.1; xmax=10; ymin=-20; ymax=110;
figure(1), subplot(211), hold on;
semilogx(normt(1,:),normt(2,:),'k');
hold off;
ylabel('dB SPL'); xlabel('Frequency (kHz)');
axis([xmin xmax ymin ymax]);
set(gca,'YTick',[0:10:100]);
set(gca,'XScale','log');
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10]);
title('All Experiments');
grid off

subplot(212), hold on;
loglog(QlowM97(:,1),QlowM97(:,2),'k-');
loglog(QhighM97(:,1),QhighM97(:,2),'k-');
hold off;
ylabel('Q10'); xlabel('Frequency (kHz)');
xlim([xmin xmax]);
ylim([.1 10]);
% set(gca,'YTick',[0 20 40 60 80 100])
set(gca,'XScale','log');
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
grid off


% Calculate audiogram
f = [0.25 0.5 1 2 3 4 6]; %kHz
audiogram = interp1(ThirdOctaveskHz*2^(1/3),AvgThresh_Normal,f) -...
    interp1(ThirdOctaveskHz*2^(1/3),AvgThresh,f);
audiogram2 = interp1(normt(1,:),normt(2,:),f) -...
    interp1(ThirdOctaveskHz*2^(1/3),MinThresh,f);
figure(2), semilogx(f,audiogram,'bo-'); hold on;
semilogx(f,audiogram2,'bo:'); hold off;
legend('Avg','Min');
xlim([min(f) max(f)]); ylim([-100 0]);
set(gca,'XTick',f,'XTickLabel',f);
title('Audiogram'); ylabel('dB HL'); xlabel('Frequency (kHz)');

Thresh_SPL_normal(1,:) = audiogram+interp1(ThirdOctaveskHz*2^(1/3),AvgThresh,f);
Thresh_SPL_normal(2,:) = audiogram2+interp1(ThirdOctaveskHz*2^(1/3),MinThresh,f);
Thresh_SPL_normal(3,:) = mean(Thresh_SPL_normal);

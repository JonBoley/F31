% clear all; close all; home;
setup_Vowel_STMP;

%% Run all units for a given experiment
if ~exist('date','var')
    date = '041811';
end
switch date
    case '041805'   % 04/18/05 fibers: (in noise, normal)
        % 1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11
        unitNums = {1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11};
    case '111804'    % 11/18/04 fibers: (in quiet, normal)
        % unitNums = {1.01, 1.02, 1.03, 1.04, 1.05, 1.08, 1.10, 1.12, 1.13,...
        %     1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.24, 1.25, 1.26, 1.27, ...
        %     1.28, 1.29, 1.30, 1.31, 1.32, 1.34, 1.36};
        % these units had >= 2levels tested:
        unitNums = {1.02, 1.22, 1.25, 1.27, 1.28, 1.30};
        % notes on 11/18/04 experiment:
        % Some units had practically no spikes for some cases (T1?)
        %   - Units 1.01, 1.10, 1.18
    case '032805'   % 03/28/05 fibers: (in quiet, normal)
        % 2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41, 2.48, 2.51, 2.53,
        % 2.54, 2.58
        unitNums = {2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41,...
            2.48, 2.51, 2.53, 2.54, 2.58};
    case '071305'   % 07/13/05 fibers: (in noise, impaired)
        % 1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27, 1.30, 1.33, 1.35,
        % 2.04, 2.07
        unitNums = {1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27,...
            1.30, 1.33, 1.35, 2.04, 2.07};
    case '121404'   % 12/14/04 fibers: (in quiet, impaired)
        % 1.03, 1.07, 1.08, 2.03, 2.04, 2.05, 2.07, 2.09, 7.01, 7.02, 7.05, 7.06,
        % 7.09, 7.10, 8.01, 9.01, 9.03, 9.05, 9.08, 9.10, 9.15, 9.16, 11.02, 13.01,
        % 13.02, 13.04

    case '011811'   % 01/18/11 fiber: (in white noise, normal but exposed)
        % 2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11
        unitNums = {2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11};
        colors = repmat({'k'},size(unitNums));
    
    case '041811'   % JB-2011_04_18-AN_Normal
        % 1.03, 2.01, 2.02, 2.07, 3.03, 3.04
        unitNums = {1.03, 2.01, 2.02, 2.07, 3.03, 3.04};
        colors = repmat({'k'},size(unitNums));
    case '051211'   % JB-2011_05_12-AN_Normal_Chin1120
        % 1.01, 1.02, 1.04
        unitNums = {1.01, 1.02, 1.04};
        colors = repmat({'k'},size(unitNums));
    case '062711'   % JB-2011_06_27-Chin1133_AN_Normal
        % 1.01, 1.02, 1.03, 1.04, 2.01, 2.02, 2.03, 2.04, 2.05, 2.08, 4.01, 4.05,
        % 4.06, 4.08, 4.09
        unitNums = {1.01, 1.02, 1.03, 1.04, 2.01, 2.02, 2.03, 2.04, 2.05, 2.08, 4.01, 4.05, 4.06, 4.08, 4.09};
        colors = repmat({'k'},size(unitNums));
    case '072011'   % JB-2011_07_20-Chin1132_AN_normal
        % 1.02, 1.03, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10
        unitNums = {1.02, 1.03, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10};
        colors = repmat({'k'},size(unitNums));

    case '062311'   % JB-2011_06_23-Chin1119_AN_500OBN
        % 1.01, 1.02, 1.09, 1.11, 1.13, 1.14
        unitNums = {1.01, 1.02, 1.09, 1.11, 1.13, 1.14};
        colors = repmat({'r'},size(unitNums));
    case '072111'   % JB-2011_07_21-Chin1124_AN_500OBN
        % 1.01, 1.02, 1.06, 1.10, 1.11, 1.12, 1.15, 1.16
        unitNums = {1.01, 1.02, 1.06, 1.10, 1.11, 1.12, 1.15, 1.16};
        colors = repmat({'r'},size(unitNums));
    case '080111'   % JB-2011_08_01-Chin1125_AN_500OBN
        % 2.02, 4.01, 4.02, 4.03, 4.05
        unitNums = {2.02, 4.01, 4.02, 4.03, 4.05};
        colors = repmat({'r'},size(unitNums));
    case '080911'   % JB-2011_08_09-Chin1135_AN_500OBN
        % 1.04, 1.05, 1.06, 1.07, 1.09, 1.10, 1.11, 1.15, 1.16
        unitNums = {1.04, 1.05, 1.06, 1.07, 1.09, 1.10, 1.11, 1.15, 1.16};
        colors = repmat({'r'},size(unitNums));
    case '081511'   % JB-2011_08_15-Chin1136_AN_500OBN
        % 1.01, 1.02, 1.05, 1.06, 1.07, 1.08, 1.10, 1.13, 3.02, 3.05, 3.06, 3.07,
        % 3.09
        unitNums = {1.01, 1.02, 1.05, 1.06, 1.07, 1.08, 1.10, 1.13, 3.02, 3.05, 3.06, 3.07, 3.09};
        colors = repmat({'r'},size(unitNums));
    case '082911'   % JB-2011_08_29-Chin1103_AN_500OBN
        unitNums = {}; %no STMP data
        colors = repmat({'r'},size(unitNums));
    case '090911'   % JB-2011_09_07-Chin1131_AN_500OBN
        unitNums = {1.04,1.06,1.09};
        colors = repmat({'r'},size(unitNums));
    case '091911'   % JB-2011_09_19-Chin1134_AN_500OBN
        unitNums = {}; % no STMP data
        colors = repmat({'r'},size(unitNums));

    case '100611'  % JB-2011_10_06-Chin1144_AN_normal
        % 1.02, 1.10, 3.06
        unitNums = {1.02, 3.06};
        colors = repmat({'k'},size(unitNums)); %normal=black
    case '101111'  % JB-2011_10_11-Chin1139_AN_normal
        % 1.01, 1.03
        unitNums = {1.01, 1.03};
        colors = repmat({'k'},size(unitNums)); %normal=black
    case '101711'  % JB-2011_10_17-Chin1151_AN_normal
        % 1.03, 2.04, 3.01, 3.02, 3.03
        unitNums = {1.03, 2.04, 3.01, 3.02, 3.03};
        colors = repmat({'k'},size(unitNums)); %normal=black

        % ultra-wide STMP (F1 & F2)
    case '112911'   % JB-2011_11_29-Chin1146_AN_normal
        % 1.09, 3.03, 4.04, 4.12, 4.13, 4.14, 4.16
        unitNums = {1.09, 3.03, 4.04, 4.12, 4.13, 4.14, 4.16};
        colors = repmat({'k'},size(unitNums));
    case '120511' %JB-2011_12_05-Chin1148_AN_normal
        % 1.09, 2.01, 2.05, 2.08, 2.11, 2.12, 2.14, 2.19, 2.23, 2.26, 2.32, 2.33, 2.34
        unitNums = {2.05, 2.08, 2.11, 2.12, 2.14, 2.19, 2.32, 2.34};
        colors = repmat({'k'},size(unitNums));

    otherwise
        error('Invalid experiment date');
end
dates=repmat({date},size(unitNums));

%% Run a combo of experiments
dates = {'041811','041811','041811','041811','041811','041811',...
    '062711','062711','062711','062711','062711','062711',...
    '072011','072011','072011','072011','072011','072011','072011',...
    '062311','062311','062311','062311',...
    '072111','072111','072111','072111','072111','072111',...
    '080111','080111','080111','080111','080111',...
    '080911','080911','080911','080911','080911','080911','080911','080911'...
    };
unitNums = {1.03,2.01,2.02,2.07,3.03,3.04,...
    1.01,1.02,1.03,1.04,2.01,2.03,...
    1.02,1.03,1.05,1.06,1.07,1.08,1.09...
    1.01,1.09,1.11,1.14,...
    1.01, 1.06, 1.11, 1.12, 1.15, 1.16,...
    2.02, 4.01, 4.02, 4.03, 4.05,...
    1.04, 1.05, 1.06, 1.07, 1.09, 1.11, 1.15, 1.16...
    };
colors = {'k','k','k','k','k','k',...
    'k','k','k','k','k','k',...
    'k','k','k','k','k','k','k',...
    'r','r','r','r',...
    'r','r','r','r','r','r',...
    'r','r','r','r','r',...
    'r','r','r','r','r','r','r','r'...
    };

%% The loop
RecalcAll = 0;  % if enabled, this automatically recalculates everything,
                % but loads previous characteristic delays
for i=1:length(unitNums)
    date = dates{i};
    disp(sprintf('%d of %d) [%s] Calculating unit %1.2f...',i,length(unitNums),date,unitNums{i}));

    % Use this for SNR conditions:
        [UnitCF(i),UnitThresh(i),UnitQ10(i),...
            Rate_failpoint(i),Rate_fail_limit(i),...
            ALSR_failpoint(i),ALSR_fail_limit(i),...
            Nscc_CD_pos_failpoint(i),Nscc_CD_pos_fail_limit(i),...
            Nscc_CD_neg_failpoint(i),Nscc_CD_neg_fail_limit(i),...
            Nscc0_pos_failpoint(i),Nscc0_pos_fail_limit(i),...
            Nscc0_neg_failpoint(i),Nscc0_neg_fail_limit(i),...
            Rho_width{i},CD_slope{i},CDatHalfOct{i},SNR{i},...
            SyncValues{i},FeatureFreqs{i}]=...
    UnitLook_EHIN_CoincDet2_JB(date,num2str(unitNums{i}),RecalcAll);

    % Use this for level conditions:
    %     [UnitCF(i),UnitThresh(i),UnitQ10(i)]=...
    %         UnitLook_EH_CoincDet2_JB(num2str(unitNums{i}));

%     if ~RecalcAll, pause; end
    close all; home;
end


if exist('CDatHalfOct') % plot CD (@0.5 oct) as a function of CF
    
end



if exist('SyncValues') % plot SyncValues as a function of CF & noise level
    % SyncValues{featureNum,snrNum} = [FeatureFreqs;BFs;Synchs];
    figure; 
    FeatureNums=[1 3 5]; %F0,F1,F2
    for i=1:length(SyncValues) % for each unit
        plotNum=0;
        for FeatureNum=FeatureNums % sync to which feature?
            for j=1:length(SNR{i})
                plotNum=plotNum+1;
                subplot(length(FeatureNums),length(SNR{i}),plotNum), hold on;
                switch j
                    case 1
                        indx_snr = find(SNR{i}==max(SNR{i}));
                    case 2
                        indx_snr = find(SNR{i}==0);
                    case 3
                        indx_snr = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
                end
                
                syncFreqsF1 = unique(SyncValues{i}{1,indx_snr}(1,:)); %frequencies for which sync was calculated
                syncFreqsT1 = unique(SyncValues{i}{2,indx_snr}(1,:));
                syncFreqsF2 = unique(SyncValues{i}{3,indx_snr}(1,:));
                syncFreqsT2 = unique(SyncValues{i}{4,indx_snr}(1,:));
                
                indx_f1 = interp1(syncFreqsF1,1:length(syncFreqsF1),...
                    FeatureFreqs{i}{1,indx_snr}(FeatureNum),'nearest','extrap');
                indx_t1 = interp1(syncFreqsT1,1:length(syncFreqsT1),...
                    FeatureFreqs{i}{2,indx_snr}(FeatureNum),'nearest','extrap');
                indx_f2 = interp1(syncFreqsF2,1:length(syncFreqsF2),...
                    FeatureFreqs{i}{3,indx_snr}(FeatureNum),'nearest','extrap');
                indx_t2 = interp1(syncFreqsT2,1:length(syncFreqsT2),...
                    FeatureFreqs{i}{4,indx_snr}(FeatureNum),'nearest','extrap');
                
                syncPlotCFsF1 = SyncValues{i}{1,indx_snr}(2,...
                    SyncValues{i}{1,indx_snr}(1,:)==syncFreqsF1(indx_f1));
                syncPlotCFsT1 = SyncValues{i}{2,indx_snr}(2,...
                    SyncValues{i}{2,indx_snr}(1,:)==syncFreqsT1(indx_t1));
                syncPlotCFsF2 = SyncValues{i}{3,indx_snr}(2,...
                    SyncValues{i}{3,indx_snr}(1,:)==syncFreqsF2(indx_f2));
                syncPlotCFsT2 = SyncValues{i}{4,indx_snr}(2,...
                    SyncValues{i}{4,indx_snr}(1,:)==syncFreqsT2(indx_t2));
                
                syncPlotValsF1 = SyncValues{i}{1,indx_snr}(3,...
                    SyncValues{i}{1,indx_snr}(1,:)==syncFreqsF1(indx_f1));
                syncPlotValsT1 = SyncValues{i}{2,indx_snr}(3,...
                    SyncValues{i}{2,indx_snr}(1,:)==syncFreqsT1(indx_t1));
                syncPlotValsF2 = SyncValues{i}{3,indx_snr}(3,...
                    SyncValues{i}{3,indx_snr}(1,:)==syncFreqsF2(indx_f2));
                syncPlotValsT2 = SyncValues{i}{4,indx_snr}(3,...
                    SyncValues{i}{4,indx_snr}(1,:)==syncFreqsT2(indx_t2));
                
                % We'll just plot sync @ CF for now
                indxCF_f1 = interp1(syncPlotCFsF1,1:length(syncPlotCFsF1),UnitCF(i),'nearest','extrap');
                indxCF_t1 = interp1(syncPlotCFsT1,1:length(syncPlotCFsT1),UnitCF(i),'nearest','extrap');
                indxCF_f2 = interp1(syncPlotCFsF2,1:length(syncPlotCFsF2),UnitCF(i),'nearest','extrap');
                indxCF_t2 = interp1(syncPlotCFsT2,1:length(syncPlotCFsT2),UnitCF(i),'nearest','extrap');
                
                hold on;
                h1=semilogx(UnitCF(i),syncPlotValsF1(indxCF_f1),[colors{i},'x']);
                if isnan(syncPlotValsF1(indxCF_f1))
                    switch colors{i}
                        case 'k', h1=semilogx(UnitCF(i),0,[colors{i},'x']); 
                        otherwise, h1=semilogx(UnitCF(i),0.01,[colors{i},'x']);
                    end
                end
                
%                 h2=semilogx(UnitCF(i),syncPlotValsT1(indxCF_t1),[colors{i},'s']);
%                 if isnan(syncPlotValsT1(indxCF_t1))
%                     switch colors{i}
%                         case 'k', h2=semilogx(UnitCF(i),0,[colors{i},'s']); 
%                         otherwise, h2=semilogx(UnitCF(i),0.01,[colors{i},'s']);
%                     end
%                 end
                
                h3=semilogx(UnitCF(i),syncPlotValsF2(indxCF_f2),[colors{i},'o']);
                if isnan(syncPlotValsF2(indxCF_f2))
                    switch colors{i}
                        case 'k', h3=semilogx(UnitCF(i),0,[colors{i},'o']); 
                        otherwise, h3=semilogx(UnitCF(i),0.01,[colors{i},'o']);
                    end
                end
                
%                 h4=semilogx(UnitCF(i),syncPlotValsT2(indxCF_t2),[colors{i},'h']); hold off;
%                 if isnan(syncPlotValsT2(indxCF_t2))
%                     switch colors{i}
%                         case 'k', h4=semilogx(UnitCF(i),0,[colors{i},'h']); 
%                         otherwise, h4=semilogx(UnitCF(i),0.01,[colors{i},'h']);
%                     end
%                 end
                
                set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
                xlim([0.35 7]); ylim([0 1]);
                
                if j==1 % apply ylabels
                    switch FeatureNum
                        case 1
                            ylabel('Sync (F0)');
                        case 3
                            ylabel('Sync (F1)');
                        case 5
                            ylabel('Sync (F2)');
                    end
                end
                switch FeatureNum % apply xlabels & titles
                    case FeatureNums(1)
                        switch j
                            case 1
                                title('In Quiet');
                            case 2
                                title('Equal SPL');
                            case 3
                                title('Equal SL');
%                                 legend([h1 h2 h3 h4],'F1 @ CF','T1 @ CF','F2 @ CF','T2 @ CF');
                                legend([h1 h3],'F1 @ CF','F2 @ CF');
%                                 legend([h2 h4],'T1 @ CF','T2 @ CF');
                        end
                    case FeatureNums(end)
                        xlabel('CF (kHz)');
                end %switch FeatureNum
            end %for i=1:length(SNR{i})
        end %for FeatureNum=FeatureNums
    end
end



if exist('Rho_width') % plot Rho & CD as a function of CF & noise level
    figure;
    FeatureNum=3;
    for i=1:length(Rho_width)
        subplot(231), hold on;
        indx = find(SNR{i}==max(SNR{i}));
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum,indx)),[colors{i},'o']);
        
        subplot(232), hold on;
        indx = find(SNR{i}==0);
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum,indx)),[colors{i},'o']);
        
        subplot(233), hold on;
        indx = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum,indx)),[colors{i},'o']);
    end
    subplot(231), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        title('In Quiet');
    ylabel('\Delta\rho@0.8 (octaves)');
    subplot(232), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        title('Equal SPL');
    subplot(233), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        title('Equal SL');
    
    
    for i=1:length(CD_slope)
        subplot(234), hold on;
        indx = find(SNR{i}==max(SNR{i}));
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum,indx)),[colors{i},'o']);
        
        subplot(235), hold on;
        indx = find(SNR{i}==0);
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum,indx)),[colors{i},'o']);
        
        subplot(236), hold on;
        indx = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum,indx)),[colors{i},'o']);
    end
    subplot(234), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        xlabel('CF (kHz)');
    ylabel('CD Slope (CF cycles/octave)');
    subplot(235), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        xlabel('CF (kHz)');
    subplot(236), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
        xlabel('CF (kHz)');
      
end



if exist('UnitCF') % Plot Thresholds/Q10s & Failpoints
    figure, subplot(211),
    load normema % load normal thresholds
    semilogx(normt(1,:),normt(2,:),'k'); hold on;
    for i=1:length(UnitCF)
        semilogx(UnitCF(i),UnitThresh(i),[colors{i} 'x']);
    end
    hold off;
    legend('Miller et al. (1997)','This Experiment');
    ylabel('Threshold (dB SPL)'); xlabel('Best Frequency (kHz)');
    axis([0.1 20 -10 110]);
    % Plot Q10
    subplot(212),
    QlowM97(:,1)=[.27 10]';
    QhighM97(:,1)=QlowM97(:,1); % Add lines from Miller et al 1997 Normal data
    QlowM97(:,2)=10^.2779*QlowM97(:,1).^.4708;
    QhighM97(:,2)=10^.6474*QlowM97(:,1).^.4708; % 5th and 95th percentiles from Miller et al, 1997
    for i=1:length(UnitCF)
        loglog(UnitCF(i),UnitQ10(i),[colors{i} 'x']); hold on;
    end
    loglog(QlowM97(:,1),QlowM97(:,2),'k-');
    loglog(QhighM97(:,1),QhighM97(:,2),'k-');
    hold off;
    ylabel('Q_{10}');  xlabel('Best Frequency (kHz)');
    axis([0.1 20 0.3 20]);

    if exist('Rate_failpoint')
        figure,
        subplot(411), hold on;
        for i=1:length(UnitCF)
            plot(UnitCF(i),Rate_failpoint(i),[colors{i} 'o']);
        end
        for i=1:length(Rate_fail_limit)
            if Rate_fail_limit(i)>0
                text(UnitCF(i),Rate_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif Rate_fail_limit(i)<0
                text(UnitCF(i),Rate_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        hold off;
        axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
        xlabel('CF');
        ylabel('SNR (dB)');
        title('Rate Failure Point');

        subplot(412), hold on;
        for i=1:length(UnitCF)
            plot(UnitCF(i),ALSR_failpoint(i),[colors{i} 'o']);
        end
        for i=1:length(ALSR_fail_limit)
            if ALSR_fail_limit(i)>0
                text(UnitCF(i),ALSR_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif ALSR_fail_limit(i)<0
                text(UnitCF(i),ALSR_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        hold off;
        axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
        xlabel('CF');
        ylabel('SNR (dB)');
        title('ALSR Failure Point');

        subplot(413), hold on;
        for i=1:length(UnitCF)
            hLine=plot(UnitCF(i),Nscc_CD_pos_failpoint(i),[colors{i} 'o']);
            if i>1
                set(get(get(hLine,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
            end
        end
        for i=1:length(Nscc_CD_pos_fail_limit)
            if Nscc_CD_pos_fail_limit(i)>0
                text(UnitCF(i),Nscc_CD_pos_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif Nscc_CD_pos_fail_limit(i)<0
                text(UnitCF(i),Nscc_CD_pos_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        for i=1:length(UnitCF)
            hLine=plot(UnitCF(i),Nscc_CD_neg_failpoint(i),[colors{i} 'o'],'MarkerFaceColor',colors{i});
            if i>1
                set(get(get(hLine,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
            end
        end
        for i=1:length(Nscc_CD_neg_fail_limit)
            if Nscc_CD_neg_fail_limit(i)>0
                text(UnitCF(i),Nscc_CD_neg_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif Nscc_CD_neg_fail_limit(i)<0
                text(UnitCF(i),Nscc_CD_neg_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        hold off;
        axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
        legend('+ slope','- slope');
        xlabel('CF');
        ylabel('SNR (dB)');
        title('Correlation Failure Point');

        subplot(414), hold on;
        for i=1:length(UnitCF)
            hLine=plot(UnitCF(i),Nscc0_pos_failpoint(i),[colors{i} 'o']);
            if i>1
                set(get(get(hLine,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
            end
        end
        for i=1:length(Nscc0_pos_fail_limit)
            if Nscc0_pos_fail_limit(i)>0
                text(UnitCF(i),Nscc0_pos_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif Nscc0_pos_fail_limit(i)<0
                text(UnitCF(i),Nscc0_pos_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        for i=1:length(UnitCF)
            hLine=plot(UnitCF(i),Nscc0_neg_failpoint(i),[colors{i} 'o'],'MarkerFaceColor',colors{i});
            if i>1
                set(get(get(hLine,'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off'); % Exclude line from legend
            end
        end
        for i=1:length(Nscc0_neg_fail_limit)
            if Nscc0_neg_fail_limit(i)>0
                text(UnitCF(i),Nscc0_neg_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center','color',colors{i});
            elseif Nscc0_neg_fail_limit(i)<0
                text(UnitCF(i),Nscc0_neg_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center','color',colors{i});
            end
        end
        hold off;
        axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
        legend('+ slope','- slope');
        xlabel('CF');
        ylabel('SNR (dB)');
        title('Delay Failure Point');
    end
end


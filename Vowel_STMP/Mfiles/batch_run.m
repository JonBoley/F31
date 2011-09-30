% clear all; close all; home;
setup_Vowel_STMP;

%% Run all units for a given experiment
% date = '041811';
% switch date
%     case '041505'   % 04/18/05 fibers: (in noise, normal)
%         % 1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11
%         unitNums = {1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11};
%     case '111804'    % 11/18/04 fibers: (in quiet, normal)
%         % unitNums = {1.01, 1.02, 1.03, 1.04, 1.05, 1.08, 1.10, 1.12, 1.13,...
%         %     1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.24, 1.25, 1.26, 1.27, ...
%         %     1.28, 1.29, 1.30, 1.31, 1.32, 1.34, 1.36};
%         % these units had >= 2levels tested:
%         unitNums = {1.02, 1.22, 1.25, 1.27, 1.28, 1.30};
%         % notes on 11/18/04 experiment:
%         % Some units had practically no spikes for some cases (T1?)
%         %   - Units 1.01, 1.10, 1.18
%     case '032805'   % 03/28/05 fibers: (in quiet, normal)
%         % 2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41, 2.48, 2.51, 2.53,
%         % 2.54, 2.58
%         unitNums = {2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41,...
%             2.48, 2.51, 2.53, 2.54, 2.58};
%     case '071305'   % 07/13/05 fibers: (in noise, impaired)
%         % 1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27, 1.30, 1.33, 1.35,
%         % 2.04, 2.07
%         unitNums = {1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27,...
%             1.30, 1.33, 1.35, 2.04, 2.07};
%     case '121404'   % 12/14/04 fibers: (in quiet, impaired)
%         % 1.03, 1.07, 1.08, 2.03, 2.04, 2.05, 2.07, 2.09, 7.01, 7.02, 7.05, 7.06,
%         % 7.09, 7.10, 8.01, 9.01, 9.03, 9.05, 9.08, 9.10, 9.15, 9.16, 11.02, 13.01,
%         % 13.02, 13.04    
%         
%     case '011811'   % 01/18/11 fiber: (in noise, normal but exposed)
%         % 2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11
%         unitNums = {2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11};
%         colors = repmat({'k'},size(unitNums));
%     case '041811'   % JB-2011_04_18-AN_Normal
%         % 1.03, 2.01, 2.02, 2.07, 3.03, 3.04
%         unitNums = {1.03, 2.01, 2.02, 2.07, 3.03, 3.04};
%         colors = repmat({'k'},size(unitNums));
%     case '051211'   % JB-2011_05_12-AN_Normal_Chin1120
%         % 1.01, 1.02, 1.04
%         unitNums = {1.01, 1.02, 1.04};
%         colors = repmat({'k'},size(unitNums));
%     case '062711'   % JB-2011_06_27-Chin1133_AN_Normal
%         % 1.01, 1.02, 1.03, 1.04, 2.01, 2.02, 2.03, 2.04, 2.05, 2.08, 4.01, 4.05,
%         % 4.06, 4.08, 4.09
%         unitNums = {1.01, 1.02, 1.03, 1.04, 2.01, 2.02, 2.03, 2.04, 2.05, 2.08, 4.01, 4.05, 4.06, 4.08, 4.09};
%         colors = repmat({'k'},size(unitNums));
%     case '072011'   % JB-2011_07_20-Chin1132_AN_normal
%         % 1.02, 1.03, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10
%         unitNums = {1.02, 1.03, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10};
%         colors = repmat({'k'},size(unitNums));
%         
%     case '062311'   % JB-2011_06_23-Chin1119_AN_500OBN
%         % 1.01, 1.02, 1.09, 1.11, 1.13, 1.14
%         unitNums = {1.01, 1.02, 1.09, 1.11, 1.13, 1.14};
%         colors = repmat({'r'},size(unitNums));
%     case '072111'   % JB-2011_07_21-Chin1124_AN_500OBN
%         % 1.01, 1.02, 1.06, 1.10, 1.11, 1.12, 1.15, 1.16
%         unitNums = {1.01, 1.02, 1.06, 1.10, 1.11, 1.12, 1.15, 1.16};
%         colors = repmat({'r'},size(unitNums));
%     case '080111'   % JB-2011_08_01-Chin1125_AN_500OBN
%         % 2.02, 4.01, 4.02, 4.03, 4.05
%         unitNums = {2.02, 4.01, 4.02, 4.03, 4.05};
%         colors = repmat({'r'},size(unitNums));
%     case '080911'   % JB-2011_08_09-Chin1135_AN_500OBN
%         % 1.04, 1.05, 1.06, 1.07, 1.09, 1.10, 1.11, 1.15, 1.16
%         unitNums = {1.04, 1.05, 1.06, 1.07, 1.09, 1.10, 1.11, 1.15, 1.16};
%         colors = repmat({'r'},size(unitNums));
%     case '081511'   % JB-2011_08_15-Chin1136_AN_500OBN
%         % 1.01, 1.02, 1.05, 1.06, 1.07, 1.08, 1.10, 1.13, 3.02, 3.05, 3.06, 3.07,
%         % 3.09
%         unitNums = {1.01, 1.02, 1.05, 1.06, 1.07, 1.08, 1.10, 1.13, 3.02, 3.05, 3.06, 3.07, 3.09};
%         colors = repmat({'r'},size(unitNums));
%         
%     otherwise
%         error('Invalid experiment date');
% end

%% Run a combo of experiments
dates = {'041811','041811','062711','062711','062711',...
    '072111','080911','081511','081511','081511'};
unitNums = {2.01,2.02,1.01,2.02,2.04,...
    1.11,1.11,1.07,1.08,3.09};
colors = {'k','k','k','k','k',...
    'r','r','r','r','r'};

%% The loop
for i=1:length(unitNums)
    date = dates{i};
    disp(sprintf('%d) [%s] Calculating unit %1.2f...',i,date,unitNums{i}));

    % Use this for SNR conditions:
    [UnitCF(i),UnitThresh(i),UnitQ10(i),...
        Rate_failpoint(i),Rate_fail_limit(i),...
        ALSR_failpoint(i),ALSR_fail_limit(i),...
        Nscc_CD_pos_failpoint(i),Nscc_CD_pos_fail_limit(i),...
        Nscc_CD_neg_failpoint(i),Nscc_CD_neg_fail_limit(i),...
        Nscc0_pos_failpoint(i),Nscc0_pos_fail_limit(i),...
        Nscc0_neg_failpoint(i),Nscc0_neg_fail_limit(i)]=...
        UnitLook_EHIN_CoincDet2_JB(date,num2str(unitNums{i}));

    % Use this for level conditions:
    %     [UnitCF(i),UnitThresh(i),UnitQ10(i)]=...
    %         UnitLook_EH_CoincDet2_JB(num2str(unitNums{i}));

    pause; close all; home;
end

% Plot Threshold
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



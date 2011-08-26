% clear all; close all; home;
setup_Vowel_STMP;

% 04/18/05 fibers: (in noise, normal)
% 1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11
% unitNums = {1.01, 1.04, 1.10, 2.04, 3.03, 4.06, 4.09, 4.11};

% 11/18/04 fibers: (in quiet, normal)
% unitNums = {1.01, 1.02, 1.03, 1.04, 1.05, 1.08, 1.10, 1.12, 1.13,...
%     1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.24, 1.25, 1.26, 1.27, ...
%     1.28, 1.29, 1.30, 1.31, 1.32, 1.34, 1.36};
% these units had >= 2levels tested:
% unitNums = {1.02, 1.22, 1.25, 1.27, 1.28, 1.30};
% notes on 11/18/04 experiment:
% Some units had practically no spikes for some cases (T1?)
%   - Units 1.01, 1.10, 1.18

% 03/28/05 fibers: (in quiet, normal)
% 2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41, 2.48, 2.51, 2.53,
% 2.54, 2.58
% unitNums = {2.11, 2.15, 2.18, 2.21, 2.28, 2.29, 2.31, 2.38, 2.41,...
%     2.48, 2.51, 2.53, 2.54, 2.58};

% 07/13/05 fibers: (in noise, impaired)
% 1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27, 1.30, 1.33, 1.35,
% 2.04, 2.07
% unitNums = {1.02, 1.04, 1.14, 1.16, 1.17, 1.20, 1.22, 1.26, 1.27,...
%     1.30, 1.33, 1.35, 2.04, 2.07};

% 12/14/04 fibers: (in quiet, impaired)
% 1.03, 1.07, 1.08, 2.03, 2.04, 2.05, 2.07, 2.09, 7.01, 7.02, 7.05, 7.06,
% 7.09, 7.10, 8.01, 9.01, 9.03, 9.05, 9.08, 9.10, 9.15, 9.16, 11.02, 13.01,
% 13.02, 13.04

% 01/18/11 fiber: (in noise, normal but exposed)
% 2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11
unitNums = {2.01, 2.02, 2.06, 2.08, 3.02, 3.08, 3.10, 3.11};

    
for i=4:8%length(unitNums)
    disp(sprintf('Calculating unit %1.2f...',unitNums{i}));
%     UnitLook_EH_CoincDet2_JB('1.07');%num2str(unitNums{i}));

    % Use this for SNR conditions:
    [UnitCF(i),UnitThresh(i),UnitQ10(i),...
        Rate_failpoint(i),Rate_fail_limit(i),...
        ALSR_failpoint(i),ALSR_fail_limit(i),...
        Nscc_CD_pos_failpoint(i),Nscc_CD_pos_fail_limit(i),...
        Nscc_CD_neg_failpoint(i),Nscc_CD_neg_fail_limit(i),...
        Nscc0_pos_failpoint(i),Nscc0_pos_fail_limit(i),...
        Nscc0_neg_failpoint(i),Nscc0_neg_fail_limit(i)]=...
        UnitLook_EHIN_CoincDet2_JB(num2str(unitNums{i}));
    
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
    semilogx(UnitCF(i),UnitThresh(i),'xk');
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
    loglog(UnitCF(i),UnitQ10(i),'xk'); hold on;
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-');
loglog(QhighM97(:,1),QhighM97(:,2),'k-');
hold off;
ylabel('Q_{10}');  xlabel('Best Frequency (kHz)');
axis([0.1 20 0.3 20]);

if exist('Rate_failpoint')
    figure,
    subplot(411),
    plot(UnitCF,Rate_failpoint,'ok'); hold on;
    for i=1:length(Rate_fail_limit)
        if Rate_fail_limit(i)>0
            text(UnitCF(i),Rate_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif Rate_fail_limit(i)<0
            text(UnitCF(i),Rate_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    hold off;
    axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
    xlabel('CF');
    ylabel('SNR (dB)');
    title('Rate Failure Point');
    
    subplot(412),
    plot(UnitCF,ALSR_failpoint,'ok'); hold on;
    for i=1:length(ALSR_fail_limit)
        if ALSR_fail_limit(i)>0
            text(UnitCF(i),ALSR_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif ALSR_fail_limit(i)<0
            text(UnitCF(i),ALSR_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    hold off;
    axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
    xlabel('CF');
    ylabel('SNR (dB)');
    title('ALSR Failure Point');
    
    subplot(413),
    plot(UnitCF,Nscc_CD_pos_failpoint,'ok'); hold on;
    for i=1:length(Nscc_CD_pos_fail_limit)
        if Nscc_CD_pos_fail_limit(i)>0
            text(UnitCF(i),Nscc_CD_pos_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif Nscc_CD_pos_fail_limit(i)<0
            text(UnitCF(i),Nscc_CD_pos_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    plot(UnitCF,Nscc_CD_neg_failpoint,'ok','MarkerFaceColor','k');
    for i=1:length(Nscc_CD_neg_fail_limit)
        if Nscc_CD_neg_fail_limit(i)>0
            text(UnitCF(i),Nscc_CD_neg_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif Nscc_CD_neg_fail_limit(i)<0
            text(UnitCF(i),Nscc_CD_neg_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    hold off;
    axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
    legend('+ slope','- slope');
    xlabel('CF');
    ylabel('SNR (dB)');
    title('Correlation Failure Point');
    
    subplot(414),
    plot(UnitCF,Nscc0_pos_failpoint,'ok'); hold on;
    for i=1:length(Nscc0_pos_fail_limit)
        if Nscc0_pos_fail_limit(i)>0
            text(UnitCF(i),Nscc0_pos_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif Nscc0_pos_fail_limit(i)<0
            text(UnitCF(i),Nscc0_pos_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    plot(UnitCF,Nscc0_neg_failpoint,'ok','MarkerFaceColor','k');
    for i=1:length(Nscc0_neg_fail_limit)
        if Nscc0_neg_fail_limit(i)>0
            text(UnitCF(i),Nscc0_neg_failpoint(i),'\uparrow','VerticalAlign','bottom','HorizontalAlign','center');
        elseif Nscc0_neg_fail_limit(i)<0
            text(UnitCF(i),Nscc0_neg_failpoint(i),'\downarrow','VerticalAlign','top','HorizontalAlign','center');
        end
    end
    hold off;
    axis([0.8*min(UnitCF) 1.2*max(UnitCF) -30 30]);
    legend('+ slope','- slope');
    xlabel('CF');
    ylabel('SNR (dB)');
    title('Delay Failure Point');
end



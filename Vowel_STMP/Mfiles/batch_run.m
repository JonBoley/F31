% clear all; close all; home;
setup_Vowel_STMP;

%% Run all units for a given experiment
if ~exist('date','var')
    date = '072112';
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

    case '050412' %JB-2012_05_04-Chin1149_AN_500OBN (with linear amplification)
        % 1.01, 1.02, 1.03, 1.05, 2.01, 2.03, 2.10
        unitNums = {1.01, 1.02, 1.03, 1.05, 2.01, 2.03};
        colors = repmat({'m'},size(unitNums));
    case '062312' %JB-2012_06_23-Chin1202_AN_500OBN (with linear amplification)
        % 1.04, 1.05, 1.07, 1.08, 1.09
        unitNums = {1.04, 1.05, 1.07, 1.08, 1.09};
        colors = repmat({'m'},size(unitNums));
    case '072112' %JB-2012_07_21-Chin1206_AN_500OBN (with linear & nonlinear amplification)
        % 1.04, 1.08, 1.12, 1.13, 3.01
        unitNums = {1.04, 1.08, 1.12, 1.13, 3.01};
        colors = repmat({'m'},size(unitNums));

    otherwise
        error('Invalid experiment date');
end
dates=repmat({date},size(unitNums));

%% Run a combo of experiments
if isempty(date) % if not any particular date, do a bunch
    dates = {'041811','041811','041811','041811','041811','041811',...
        '062711','062711','062711','062711','062711','062711',...
        '072011','072011','072011','072011','072011','072011','072011',...
        '062311','062311','062311','062311',...
        '072111','072111','072111','072111','072111','072111',...
        '080111','080111','080111','080111','080111',...
        '080911','080911','080911','080911','080911','080911','080911','080911',...
        '050412','050412','050412','050412','050412','050412'...
        };
    unitNums = {1.03,2.01,2.02,2.07,3.03,3.04,...
        1.01,1.02,1.03,1.04,2.01,2.03,...
        1.02,1.03,1.05,1.06,1.07,1.08,1.09...
        1.01,1.09,1.11,1.14,...
        1.01, 1.06, 1.11, 1.12, 1.15, 1.16,...
        2.02, 4.01, 4.02, 4.03, 4.05,...
        1.04, 1.05, 1.06, 1.07, 1.09, 1.11, 1.15, 1.16,...
        1.01, 1.02, 1.03, 1.05, 2.01, 2.03...
        };
    colors = {'k','k','k','k','k','k',...
        'k','k','k','k','k','k',...
        'k','k','k','k','k','k','k',...
        'r','r','r','r',...
        'r','r','r','r','r','r',...
        'r','r','r','r','r',...
        'r','r','r','r','r','r','r','r',...
        'm','m','m','m','m','m'...
        };
end

%% The loop
RecalcAll = 0;  % if enabled, this automatically recalculates everything,
                % but loads previous characteristic delays
LoadMAT = 1; % load mat file instead of running UnitLook_EHIN_CoincDet2_JB

if LoadMAT
    [FileName,PathName,FilterIndex] = uigetfile('*.mat','Pick a file',...
        'C:\Research\MATLAB\Vowel_STMP\ExpData\batch_053112.mat');
    if FileName
        LoadMATbackup = LoadMAT;
        load(fullfile(PathName,FileName),'-regexp','[^h1]');
        LoadMAT = LoadMATbackup;
    else
        LoadMAT=0;
    end
end
if ~LoadMAT
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
            Rho_width{i},CD_slope{i},CD_slope2{i},CDatHalfOct{i},SNR{i},...
            SyncValues{i},FeatureFreqs{i}]=...
            UnitLook_EHIN_CoincDet2_JB(date,num2str(unitNums{i}),RecalcAll);

        % Use this for level conditions:
        %     [UnitCF(i),UnitThresh(i),UnitQ10(i)]=...
        %         UnitLook_EH_CoincDet2_JB(num2str(unitNums{i}));

        %     if ~RecalcAll, pause; end
        close all; home;
    end
end

% use CD_slope2
CD_slope_bak = CD_slope;
CD_slope = CD_slope2;

if 0%exist('CDatHalfOct') % plot CD (@0.5 oct) as a function of CF
    figure;
    FeatureNums = [1 3]; %F1,F2
    
    arrayCDatHalfOct.nh = cell(length(FeatureNums),length(SNR{1})); %normal
    arrayCDatHalfOct.snhl = cell(length(FeatureNums),length(SNR{1})); %SNHL
    arrayCDatHalfOct.linear = cell(length(FeatureNums),length(SNR{1})); %linear
    arrayCDatHalfOct.nonlinear = cell(length(FeatureNums),length(SNR{1})); %nonlinear
    for i=1:length(FeatureNums)
        for j=1:length(SNR{1})
            arrayCDatHalfOct.nh{i,j} = NaN*ones(length(CDatHalfOct),2);
            arrayCDatHalfOct.snhl{i,j} = NaN*ones(length(CDatHalfOct),2);
            arrayCDatHalfOct.linear{i,j} = NaN*ones(length(CDatHalfOct),2);
            arrayCDatHalfOct.nonlinear{i,j} = NaN*ones(length(CDatHalfOct),2);
        end
    end
    
    
    
    for i=1:length(CDatHalfOct) % for each unit
        plotNum=0;
        for FeatureNum=FeatureNums % sync to which feature?
            for j=1:length(SNR{i})
                plotNum=plotNum+1;
                subplot(length(FeatureNums),length(SNR{i}),plotNum), hold on;
                switch j % find indx_snr
                    case 1
                        indx_snr = find(SNR{i}==max(SNR{i}));
                    case 2
                        indx_snr = find(SNR{i}==0);
                    case 3
                        indx_snr = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
                end
                
                % convert from us to CF cycles
                tempCD = CDatHalfOct{i}{find(FeatureNum==FeatureNums),indx_snr};
                CF_cycle_usec = 1e6/(UnitCF(i)*1e3);
                tempCD = tempCD./CF_cycle_usec;
                
                
                switch colors{i}
                    case {'k'} % normal
                        arrayCDatHalfOct.nh{find(FeatureNum==FeatureNums),indx_snr}(i,:) =...
                            [UnitCF(i), tempCD];
                    case {'r'} % impaired
                        arrayCDatHalfOct.snhl{find(FeatureNum==FeatureNums),indx_snr}(i,:) =...
                            [UnitCF(i), tempCD];
                    case {'m'} % impaired + linear amplification
                        arrayCDatHalfOct.linear{find(FeatureNum==FeatureNums),indx_snr}(i,:) =...
                            [UnitCF(i), tempCD];
                    case {'g'} % impaired + nonlinear amplification
                        arrayCDatHalfOct.nonlinear{find(FeatureNum==FeatureNums),indx_snr}(i,:) =...
                            [UnitCF(i), tempCD];
                end
                
                if exist('h1'), cla(h1); cla(h2); cla(h3); cla(h4); end
                h1=semilogx(arrayCDatHalfOct.nh{find(FeatureNum==FeatureNums),indx_snr}(:,1),...
                    arrayCDatHalfOct.nh{find(FeatureNum==FeatureNums),indx_snr}(:,2),'k.');
                h2=semilogx(arrayCDatHalfOct.snhl{find(FeatureNum==FeatureNums),indx_snr}(:,1),...
                    arrayCDatHalfOct.snhl{find(FeatureNum==FeatureNums),indx_snr}(:,2),'r.');
                h3=semilogx(arrayCDatHalfOct.linear{find(FeatureNum==FeatureNums),indx_snr}(:,1),...
                    arrayCDatHalfOct.linear{find(FeatureNum==FeatureNums),indx_snr}(:,2),'m.');
                h4=semilogx(arrayCDatHalfOct.nonlinear{find(FeatureNum==FeatureNums),indx_snr}(:,1),...
                    arrayCDatHalfOct.nonlinear{find(FeatureNum==FeatureNums),indx_snr}(:,2),'g.');
                
                set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
                xlim([0.1 10]); ylim([0 10]);
                
                if i==1
                    h_text(plotNum)=text(1,9,'[empty]');
                    set(h_text(plotNum),'HandleVisibility','off');
                end
                
                % smooth the data
                smoothOct = 0.7; % number of octaves to smooth over
                xData = arrayCDatHalfOct.nh{find(FeatureNum==FeatureNums),indx_snr}(:,1);
                yData = arrayCDatHalfOct.nh{find(FeatureNum==FeatureNums),indx_snr}(:,2);
                indexNaN = isnan(xData); xData = xData(~indexNaN); yData = yData(~indexNaN);
                indexNaN = isnan(yData); xData = xData(~indexNaN); yData = yData(~indexNaN);
                [xData,IX]=unique(xData); yData=yData(IX); %warning: this throws away duplicate CFs
                if length(xData)>2
                    clear xData_smoothed yData_smoothed;
                    xData_smoothed = 100e-3*2.^(0:smoothOct/4:6); yData_smoothed=NaN*ones(size(xData_smoothed));
                    for q=1:length(xData_smoothed)
                        temp_ind = xData<=(xData_smoothed(q)*2^(smoothOct/2)) & xData>=(xData_smoothed(q)*2^(-smoothOct/2));
                        distance_oct = abs(log2(xData(temp_ind)/xData_smoothed(q)));
                        triangular_weights=-(abs(distance_oct)-(smoothOct/2))/(smoothOct/2);
                        yData_smoothed(q) = median(yData(temp_ind));
                    end
                    semilogx(xData_smoothed,yData_smoothed,'k');
                    
%                     cfun1 = fit(xData,yData,'linearinterp');
%                     plot(cfun1,'k'); drawnow;
                end
                xData2 = arrayCDatHalfOct.snhl{find(FeatureNum==FeatureNums),indx_snr}(:,1);
                yData2 = arrayCDatHalfOct.snhl{find(FeatureNum==FeatureNums),indx_snr}(:,2);
                indexNaN = isnan(xData2); xData2 = xData2(~indexNaN); yData2 = yData2(~indexNaN);
                indexNaN = isnan(yData2); xData2 = xData2(~indexNaN); yData2 = yData2(~indexNaN);
                [xData2,IX]=unique(xData2); yData2=yData2(IX); %warning: this throws away duplicate CFs
                if length(xData2)>2
                    clear xData2_smoothed yData2_smoothed;
                    xData2_smoothed = 100e-3*2.^(0:smoothOct/4:6); yData2_smoothed=NaN*ones(size(xData2_smoothed));
                    for q=1:length(xData2_smoothed)
                        temp_ind = xData2<=(xData2_smoothed(q)*2^(smoothOct/2)) & xData2>=(xData2_smoothed(q)*2^(-smoothOct/2));
                        distance_oct = abs(log2(xData2(temp_ind)/xData2_smoothed(q)));
                        triangular_weights=-(abs(distance_oct)-(smoothOct/2))/(smoothOct/2);
                        yData2_smoothed(q) = median(yData2(temp_ind));
                    end
                    semilogx(xData2_smoothed,yData2_smoothed,'r');
                    
                    % avg CD within a half-octave of 1.4kHz
                    avgCD_nh = yData_smoothed(abs(log2(xData2_smoothed/1.414))<0.5);
                    avgCD_snhl = yData2_smoothed(abs(log2(xData2_smoothed/1.414))<0.5);
                    deltaCD = avgCD_snhl - avgCD_nh;
                    set(h_text(plotNum),'Interpreter','latex');
                    set(h_text(plotNum),'String',[sprintf('From 1kHz to 2kHz: \n'),...
                        '$\overline{CD_{NH}}$=' sprintf('%1.2f cycles\n',mean(avgCD_nh(~isnan(avgCD_nh)))),...
                        '$\overline{CD_{SNHL}}$=' sprintf('%1.2f cycles\n',mean(avgCD_snhl(~isnan(avgCD_snhl)))),...
                        '$\overline{\Delta CD}$=' sprintf('%1.2f cycles',mean(deltaCD(~isnan(deltaCD))))]);
                    set(h_text(plotNum),'HorizontalAlignment','center');
                    
                    % avg [of trend line]...
                    avgCD2_nh(FeatureNum,j) = mean(avgCD_nh(~isnan(avgCD_nh)));
                    avgCD2_snhl(FeatureNum,j) = mean(avgCD_snhl(~isnan(avgCD_snhl)));
                    
                    % stderr [of actual points]
                    temp = yData(abs(log2(xData/1.414))<0.5);
                    temp2 = yData(abs(log2(xData2/1.414))<0.5);
                    stderrCD_nh(FeatureNum,j) = ...
                        sqrt(sum((temp(~isnan(temp))-mean(temp(~isnan(temp)))).^2/length(temp)))/sqrt(length(temp));
                    stderrCD_snhl(FeatureNum,j) = ...
                        sqrt(sum((temp2(~isnan(temp2))-mean(temp2(~isnan(temp2)))).^2/length(temp2)))/sqrt(length(temp2));

%                     cfun2 = fit(xData,yData,'linearinterp');
%                     plot(cfun2,'r'); drawnow;
                end
                
                xData3 = arrayCDatHalfOct.linear{find(FeatureNum==FeatureNums),indx_snr}(:,1);
                yData3 = arrayCDatHalfOct.linear{find(FeatureNum==FeatureNums),indx_snr}(:,2);
                indexNaN = isnan(xData3); xData3 = xData3(~indexNaN); yData3 = yData3(~indexNaN);
                indexNaN = isnan(yData3); xData3 = xData3(~indexNaN); yData3 = yData3(~indexNaN);
                [xData3,IX]=unique(xData3); yData3=yData3(IX); %warning: this throws away duplicate CFs
                if length(xData3)>2
                    clear xData3_smoothed yData3_smoothed;
                    xData3_smoothed = 100e-3*2.^(0:smoothOct/4:6); yData3_smoothed=NaN*ones(size(xData3_smoothed));
                    for q=1:length(xData3_smoothed)
                        temp_ind = xData3<=(xData3_smoothed(q)*2^(smoothOct/2)) & xData3>=(xData3_smoothed(q)*2^(-smoothOct/2));
                        distance_oct = abs(log2(xData3(temp_ind)/xData3_smoothed(q)));
                        triangular_weights=-(abs(distance_oct)-(smoothOct/2))/(smoothOct/2);
                        yData3_smoothed(q) = median(yData3(temp_ind));
                    end
                    semilogx(xData3_smoothed,yData3_smoothed,'m');
                end
                
                xData4 = arrayCDatHalfOct.nonlinear{find(FeatureNum==FeatureNums),indx_snr}(:,1);
                yData4 = arrayCDatHalfOct.nonlinear{find(FeatureNum==FeatureNums),indx_snr}(:,2);
                indexNaN = isnan(xData4); xData4 = xData4(~indexNaN); yData4 = yData4(~indexNaN);
                indexNaN = isnan(yData4); xData4 = xData4(~indexNaN); yData4 = yData4(~indexNaN);
                [xData4,IX]=unique(xData4); yData4=yData4(IX); %warning: this throws away duplicate CFs
                if length(xData4)>2
                    clear xData4_smoothed yData4_smoothed;
                    xData4_smoothed = 100e-3*2.^(0:smoothOct/4:6); yData4_smoothed=NaN*ones(size(xData4_smoothed));
                    for q=1:length(xData4_smoothed)
                        temp_ind = xData4<=(xData4_smoothed(q)*2^(smoothOct/2)) & xData4>=(xData4_smoothed(q)*2^(-smoothOct/2));
                        distance_oct = abs(log2(xData4(temp_ind)/xData4_smoothed(q)));
                        triangular_weights=-(abs(distance_oct)-(smoothOct/2))/(smoothOct/2);
                        yData4_smoothed(q) = median(yData4(temp_ind));
                    end
                    semilogx(xData4_smoothed,yData4_smoothed,'g');
                end
                
                
                if j==1 % apply ylabels
                    switch FeatureNum
                        case 1
                            ylabel(sprintf('CD (CF cycles)\n[0.5 oct re F1]'));
                        case 3
                            ylabel(sprintf('CD (CF cycles)\n[0.5 oct re F2]'));
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
%                                 legend([h1 h3],'F1 @ CF','F2 @ CF');
%                                 legend([h2 h4],'T1 @ CF','T2 @ CF');
                        end
                    case FeatureNums(end)
                        xlabel('CF (kHz)');
                end %switch FeatureNum
                
                
            end
        end 
    end
    save('C:\Research\MATLAB\Vowel_STMP\ExpData\regressionTest.mat','arrayCDatHalfOct');
end


if 0%exist('SyncValues') % plot SyncValues as a function of CF & noise level
    % SyncValues{featureNum,snrNum} = [FeatureFreqs;BFs;Synchs];
    array_syncPlotValsF1.nh = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsT1.nh = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsF2.nh = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsT2.nh = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsF1.snhl = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsT1.snhl = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsF2.snhl = cell(length(FeatureNums),length(SNR{1}));
    array_syncPlotValsT2.snhl = cell(length(FeatureNums),length(SNR{1}));
    
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
                
                switch colors{i}
                    case 'k'
                        array_syncPlotValsF1.nh{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsF1(indxCF_f1);
                        array_syncPlotValsT1.nh{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsT1(indxCF_t1);
                        array_syncPlotValsF2.nh{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsF2(indxCF_f2);
                        array_syncPlotValsT2.nh{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsT2(indxCF_t2);
                        
                        array_syncPlotValsF1.snhl{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsT1.snhl{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsF2.snhl{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsT2.snhl{find(FeatureNum==FeatureNums),j}(i) = NaN;
                    otherwise
                        array_syncPlotValsF1.nh{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsT1.nh{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsF2.nh{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        array_syncPlotValsT2.nh{find(FeatureNum==FeatureNums),j}(i) = NaN;
                        
                        array_syncPlotValsF1.snhl{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsF1(indxCF_f1);
                        array_syncPlotValsT1.snhl{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsT1(indxCF_t1);
                        array_syncPlotValsF2.snhl{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsF2(indxCF_f2);
                        array_syncPlotValsT2.snhl{find(FeatureNum==FeatureNums),j}(i) = syncPlotValsT2(indxCF_t2);
                end
                array_syncPlotValsF1.CFs(i) = UnitCF(i);
                array_syncPlotValsT1.CFs(i) = UnitCF(i);
                array_syncPlotValsF2.CFs(i) = UnitCF(i);
                array_syncPlotValsT2.CFs(i) = UnitCF(i);

                
                hold on;
                h1=semilogx(UnitCF(i),syncPlotValsF1(indxCF_f1),[colors{i},'+']);
                if isnan(syncPlotValsF1(indxCF_f1))
                    switch colors{i}
                        case 'k', h1=semilogx(UnitCF(i),0,[colors{i},'+']); 
                        otherwise, h1=semilogx(UnitCF(i),0.01,[colors{i},'+']);
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


if 0%exist('array_syncPlotValsF1') % group sync by CF range
    % array_syncPlotValsF1.snhl{find(FeatureNum==FeatureNums),j}(i)
    LocalSymbols = ['x','s','h'];
    centerCFs = [0.5 1 2 4]; %kHz
    CFrange = 1; % octaves
    CFranges = [centerCFs'*2^(-CFrange/2), centerCFs'*2^(CFrange/2)];
    
    figure; plotNum=0;
    for FeatureIndex=1:size(array_syncPlotValsF1.nh,1)
        for SNRIndex=1:size(array_syncPlotValsF1.nh,2)
            plotNum=plotNum+1;
            subplot(length(FeatureNums),length(SNR{i}),plotNum), hold on;

            for centerCFIndex=1:length(centerCFs)
                CFindices = find(array_syncPlotValsF1.CFs>=CFranges(centerCFIndex,1) &...
                    array_syncPlotValsF1.CFs<CFranges(centerCFIndex,2));
                
                %%%% F1
                % replace NaN's with zeros
                NaN_NHindices = find(isnan(array_syncPlotValsF1.nh{FeatureIndex,SNRIndex}));
                NaN_SNHLindices = find(isnan(array_syncPlotValsF1.snhl{FeatureIndex,SNRIndex}));
                array_syncPlotValsF1.nh{FeatureIndex,SNRIndex}(intersect(CFindices,NaN_NHindices)) = 0;
                array_syncPlotValsF1.snhl{FeatureIndex,SNRIndex}(intersect(CFindices,NaN_SNHLindices)) = 0;
                % save mean & std error
                syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,1)=...
                    mean(array_syncPlotValsF1.nh{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,2)=...
                    std(array_syncPlotValsF1.nh{FeatureIndex,SNRIndex}(CFindices))/...
                    numel(array_syncPlotValsF1.nh{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1)=...
                    mean(array_syncPlotValsF1.snhl{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2)=...
                    std(array_syncPlotValsF1.snhl{FeatureIndex,SNRIndex}(CFindices))/...
                    numel(array_syncPlotValsF1.snhl{FeatureIndex,SNRIndex}(CFindices));
                % plot
                h_syncCF(1)=errorbar(centerCFs(centerCFIndex),...
                    syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,1),...
                    syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,2),...
                    'k+');
                errorbar(centerCFs(centerCFIndex),...
                    max(0.01,syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1)),...
                    syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2),...
                    'r+');
                
                %%%% F2
                % replace NaN's with zeros
                NaN_NHindices = find(isnan(array_syncPlotValsF2.nh{FeatureIndex,SNRIndex}));
                NaN_SNHLindices = find(isnan(array_syncPlotValsF2.snhl{FeatureIndex,SNRIndex}));
                array_syncPlotValsF2.nh{FeatureIndex,SNRIndex}(intersect(CFindices,NaN_NHindices)) = 0;
                array_syncPlotValsF2.snhl{FeatureIndex,SNRIndex}(intersect(CFindices,NaN_SNHLindices)) = 0;
                % save mean & std error
                syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,1)=...
                    mean(array_syncPlotValsF2.nh{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,2)=...
                    std(array_syncPlotValsF2.nh{FeatureIndex,SNRIndex}(CFindices))/...
                    numel(array_syncPlotValsF2.nh{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1)=...
                    mean(array_syncPlotValsF2.snhl{FeatureIndex,SNRIndex}(CFindices));
                syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2)=...
                    std(array_syncPlotValsF2.snhl{FeatureIndex,SNRIndex}(CFindices))/...
                    numel(array_syncPlotValsF2.snhl{FeatureIndex,SNRIndex}(CFindices));
                % plot
                h_syncCF(1)=errorbar(centerCFs(centerCFIndex),...
                    syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,1),...
                    syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,2),...
                    'ko');
                errorbar(centerCFs(centerCFIndex),...
                    max(0.01,syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1)),...
                    syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2),...
                    'ro');
                
            end
            
            set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
            xlim([0.35 7]); ylim([0 1]);
            
            if SNRIndex==1 % apply ylabels
                switch FeatureIndex
                    case 1
                        ylabel('Sync (F0)');
                    case 2
                        ylabel('Sync (F1)');
                    case 3
                        ylabel('Sync (F2)');
                end
            end
            switch FeatureIndex % apply xlabels & titles
                case 1
                    switch SNRIndex
                        case 1
                            title('In Quiet');
                        case 2
                            title('Equal SPL');
                        case 3
                            title('Equal SL');
                            legend(h_syncCF,'F1 @ CF','F2 @ CF');
                    end
                case FeatureNums(end)
                    xlabel('CF (kHz)');
            end %switch FeatureNum
        end
    end
end
% set(gcf,'visible','off');

if 0%exist('syncCFgroups') % plot sync as a function noise level (grouped by CF)
    % syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex)
    figure; 
    plotNum=0;
    for FeatureIndex=1:size(syncCFgroups.F1.nh,1)
        for centerCFIndex=1:length(syncCFgroups.F2.snhl{FeatureIndex,1})
            plotNum=plotNum+1;
            subplot(length(FeatureNums),length(syncCFgroups.F2.snhl{FeatureIndex,1}),plotNum), hold on;
            for SNRIndex=1:size(syncCFgroups.F1.nh,2)
                tmpF1_nh(SNRIndex,1) = syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,1);
                tmpF1_snhl(SNRIndex,1) = syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1);
                tmpF2_nh(SNRIndex,1) = syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,1);
                tmpF2_snhl(SNRIndex,1) = syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,1);
                
                tmpF1_nh(SNRIndex,2) = syncCFgroups.F1.nh{FeatureIndex,SNRIndex}(centerCFIndex,2);
                tmpF1_snhl(SNRIndex,2) = syncCFgroups.F1.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2);
                tmpF2_nh(SNRIndex,2) = syncCFgroups.F2.nh{FeatureIndex,SNRIndex}(centerCFIndex,2);
                tmpF2_snhl(SNRIndex,2) = syncCFgroups.F2.snhl{FeatureIndex,SNRIndex}(centerCFIndex,2);
            end %SNRIndex
            h1=errorbar(tmpF1_nh(:,1),tmpF1_nh(:,2),'k+-'); hold on;
            errorbar(tmpF1_snhl(:,1),tmpF1_snhl(:,2),'r+-');
            h2=errorbar(tmpF2_nh(:,1),tmpF2_nh(:,2),'ko-');
            errorbar(tmpF2_snhl(:,1),tmpF2_snhl(:,2),'ro-'); hold off;
            
            xlim([0.5 3.5]); ylim([0 1]);
            set(gca,'XTick',[1 2 3]);
            set(gca,'XTickLabel',{'Quiet','SPL','SL'});
            
            if centerCFIndex==1 % apply ylabels
                switch FeatureIndex
                    case 1
                        ylabel('Sync to F0');
                    case 2
                        ylabel('Sync to F1');
                    case 3
                        ylabel('Sync to F2');
                end
            end
            switch FeatureIndex % apply xlabels & titles
                case 1
                    switch centerCFIndex
                        case 1
                            title('CF\approx500Hz (\pm0.5oct)');
                        case 2
                            title('CF\approx1kHz (\pm0.5oct)');
                        case 3
                            title('CF\approx2kHz (\pm0.5oct)');
                        case 4
                            title('CF\approx4kHz (\pm0.5oct)');
                            legend([h1,h2],'F1 @ CF','F2 @ CF');
                    end
                case FeatureNums(end)
                    %N/A
            end %switch FeatureNum
            
        end %centerCFIndex
    end %FeatureIndex
    
end


if exist('Rho_width') % plot Rho & CD as a function of CF & noise level
    figure;
    FeatureNum=3; %F2
    for i=1:length(Rho_width)
        FeatureNum2 = min(FeatureNum,size(Rho_width{i},1));
        
        subplot(231), hold on;
        indx = find(SNR{i}==max(SNR{i}));
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum2,indx)),[colors{i},'o']);

        subplot(232), hold on;
        indx = find(SNR{i}==0);
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum2,indx)),[colors{i},'o']);

        subplot(233), hold on;
        indx = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
        semilogx(UnitCF(i),cell2mat(Rho_width{i}(FeatureNum2,indx)),[colors{i},'o']);
    end
    subplot(231), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    title('In Quiet');
    ylabel('\Delta\rho@0.8 (octaves)');
    subplot(232), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    title('Equal SPL');
    subplot(233), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    title('Equal SL');


    for i=1:length(CD_slope)
        FeatureNum2 = min(FeatureNum,size(CD_slope{i},1));
        
        subplot(234), hold on;
        indx = find(SNR{i}==max(SNR{i}));
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum2,indx)),[colors{i},'o']);

        subplot(235), hold on;
        indx = find(SNR{i}==0);
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum2,indx)),[colors{i},'o']);

        subplot(236), hold on;
        indx = find(SNR{i}~=max(SNR{i}) & SNR{i}~=0);
        semilogx(UnitCF(i),(CD_slope{i}(FeatureNum2,indx)),[colors{i},'o']);
    end
    subplot(234), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    xlabel('CF (kHz)');
    ylabel('CD Slope (CF cycles/octave)');
    subplot(235), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    xlabel('CF (kHz)');
    subplot(236), set(gca,'XScale','log'); set(gca,'XTick',[1 2 4]);
    xlabel('CF (kHz)');

    %%% save rho & CD
    array_rho.nh=cell(2,length(SNR{1}));
    array_cd.nh=cell(2,length(SNR{1}));
    array_rho.snhl=cell(2,length(SNR{1}));
    array_cd.snhl=cell(2,length(SNR{1}));
    nhCFs=UnitCF(cell2mat(colors)=='k'); 
    snhlCFs=UnitCF(cell2mat(colors)~='k');
    FeatureIndex=0;
    for FeatureNum=[1 3]; %F1,F2
        FeatureIndex=FeatureIndex+1;
        for CFindex=1:length(Rho_width)
            FeatureNum2 = min(FeatureNum,size(Rho_width{CFindex},1));
            indx = find(SNR{CFindex}==max(SNR{CFindex}));
            if colors{CFindex}=='k'
                array_rho.nh{FeatureIndex,1}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.nh{FeatureIndex,1}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);
            else
                array_rho.snhl{FeatureIndex,1}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.snhl{FeatureIndex,1}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);                
            end

            indx = find(SNR{CFindex}==0);
            if colors{CFindex}=='k'
                array_rho.nh{FeatureIndex,2}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.nh{FeatureIndex,2}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);
            else
                array_rho.snhl{FeatureIndex,2}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.snhl{FeatureIndex,2}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);
            end

            indx = find(SNR{CFindex}~=max(SNR{CFindex}) & SNR{CFindex}~=0);
            if colors{CFindex}=='k'
                array_rho.nh{FeatureIndex,3}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.nh{FeatureIndex,3}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);
            else
                array_rho.snhl{FeatureIndex,3}(end+1) = cell2mat(Rho_width{CFindex}(FeatureNum2,indx));
                array_cd.snhl{FeatureIndex,3}(end+1) = CD_slope{CFindex}(FeatureNum2,indx);
            end
        end
    end
end


if exist('array_rho') % plot Rho & CD as a function noise level (grouped by CF)
    centerCFs = [0.5 1 2 4]; %kHz
    CFrange = 1; % octaves
    CFranges = [centerCFs'*2^(-CFrange/2), centerCFs'*2^(CFrange/2)];
    
    figure(998); 
    figure(999);
    plotNum=0;
    for FeatureNum=1:size(array_rho.nh,1)
        for centerCFIndex=1:length(centerCFs)
            plotNum=plotNum+1;
            
            for SNRindex=1:size(array_rho.nh,2)
                %NH
                CFindices1 = find(nhCFs>=CFranges(centerCFIndex,1) &...
                    nhCFs<CFranges(centerCFIndex,2));
                
                rhoCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,1) = ...
                    mean(array_rho.nh{FeatureNum,SNRindex}(CFindices1));
                rhoCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,2) = ...
                    std(array_rho.nh{FeatureNum,SNRindex}(CFindices1))/...
                    numel(array_rho.nh{FeatureNum,SNRindex}(CFindices1));
                tmpRho_nh(SNRindex,:) = rhoCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,:);
                
                cdCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,1) = ...
                    mean(array_cd.nh{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.nh{FeatureNum,SNRindex})),CFindices1)));
                cdCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,2) = ...
                    std(array_cd.nh{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.nh{FeatureNum,SNRindex})),CFindices1)))/...
                    numel(array_cd.nh{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.nh{FeatureNum,SNRindex})),CFindices1)));
                tmpCD_nh(SNRindex,:) = cdCFgroups.nh{FeatureNum,SNRindex}(centerCFIndex,:);
                
                %SNHL
                CFindices2 = find(snhlCFs>=CFranges(centerCFIndex,1) &...
                    snhlCFs<CFranges(centerCFIndex,2));
                
                rhoCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,1) = ...
                    mean(array_rho.snhl{FeatureNum,SNRindex}(CFindices2));
                rhoCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,2) = ...
                    std(array_rho.snhl{FeatureNum,SNRindex}(CFindices2))/...
                    numel(array_rho.snhl{FeatureNum,SNRindex}(CFindices2));
                tmpRho_snhl(SNRindex,:) = rhoCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,:);
                
                cdCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,1) = ...
                    mean(array_cd.snhl{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.snhl{FeatureNum,SNRindex})),CFindices2)));
                cdCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,2) = ...
                    std(array_cd.snhl{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.snhl{FeatureNum,SNRindex})),CFindices2)))/...
                    numel(array_cd.snhl{FeatureNum,SNRindex}(intersect(find(~isnan(array_cd.snhl{FeatureNum,SNRindex})),CFindices2)));
                tmpCD_snhl(SNRindex,:) = cdCFgroups.snhl{FeatureNum,SNRindex}(centerCFIndex,:);
            end
            
            figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum), hold on;
            errorbar(tmpRho_nh(:,1),tmpRho_nh(:,2),'k-');
            errorbar(tmpRho_snhl(:,1),tmpRho_snhl(:,2),'r-');
            set(gca,'XTick',[1 2 3]);
            set(gca,'XTickLabel',{'Quiet','SPL','SL'});
            ylim([0 0.5]);
            legend(sprintf('%d NH units',length(CFindices1)),...
                sprintf('%d SNHL units',length(CFindices2)));
            
            figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum), hold on;
            errorbar(tmpCD_nh(:,1),tmpCD_nh(:,2),'k-');
            errorbar(tmpCD_snhl(:,1),tmpCD_snhl(:,2),'r-');
            set(gca,'XTick',[1 2 3]);
            set(gca,'XTickLabel',{'Quiet','SPL','SL'});
            ylim([0 50]);
            legend(sprintf('%d NH units',length(CFindices1)),...
                sprintf('%d SNHL units',length(CFindices2)));
            
            if centerCFIndex==1 % apply ylabels
                switch FeatureNum
                    case 1
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        ylabel('\rho width (octaves re F1)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        ylabel('CD slope (cycles/oct, re F1)');
                    case 2
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        ylabel('\rho width (octaves re F2)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        ylabel('CD slope (cycles/oct, re F2)');
                end
            end
            if FeatureNum==1 % apply xlabels & titles
                switch centerCFIndex
                    case 1
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx500Hz (\pm0.5oct)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx500Hz (\pm0.5oct)');
                    case 2
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx1kHz (\pm0.5oct)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx1kHz (\pm0.5oct)');
                    case 3
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx2kHz (\pm0.5oct)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx2kHz (\pm0.5oct)');
                    case 4
                        figure(998), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx4kHz (\pm0.5oct)');
                        figure(999), subplot(size(array_rho.nh,1),length(centerCFs),plotNum),
                        title('CF\approx4kHz (\pm0.5oct)');
                end
            end
            
        end
    end
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

figure, subplot(121),
errorbar([1:3],avgCD2_nh(1,:),stderrCD_nh(1,:),'ko-'); hold on;
errorbar([1:3],avgCD2_snhl(1,:),stderrCD_snhl(1,:),'ro-');
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Quiet','SPL','SL'});
ylabel('CD re Normal (CF cycles) [0.5 oct re F1]');
title('Strength of Spatiotemporal Coding (CF=1-2kHz)');
legend('Normal','Impaired'); ylim([0 3]);
subplot(122),
errorbar([1:3],avgCD2_nh(3,:),stderrCD_nh(3,:),'ko-'); hold on;
errorbar([1:3],avgCD2_snhl(3,:),stderrCD_snhl(3,:),'ro-');
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Quiet','SPL','SL'});
ylabel('CD (CF cycles) [0.5 oct re F2]');
title('Strength of Spatiotemporal Coding (CF=1-2kHz)');
legend('Normal','Impaired'); ylim([0 3]);

figure
plot([1:3],avgCD2_snhl(1,:)-avgCD2_nh(1,:),'ro-'); hold on;
plot([1:3],avgCD2_snhl(3,:)-avgCD2_nh(3,:),'rsq-'); hold off;
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Quiet','SPL','SL'});
ylabel('CD re Normal (CF cycles) [0.5 oct re F1]');
title('Reduction of Spatiotemporal Coding (CF=1-2kHz)');
legend('F1','F2'); ylim([-1 0]);


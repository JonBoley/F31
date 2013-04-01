%% findDelay.m
% plots the PST histograms for a collection of neural data
% so that a curve can be fit to the conduction delay vs. CF

% clear all; close all; home;
setup_Vowel_STMP;

config = 'amplified'; %{'normal','impaired','amplified'}

global ROOT_dir
global SavedPICS SavedPICnums SavedPICSuse
global FeaturesText FormsAtHarmonicsText InvertPolarityText

path(fullfile(ROOT_dir,'MFiles','NewMfiles'),path);

switch config
    case {'normal'}
        dates = {'041811','041811','041811','041811','041811','041811',...
            '062711','062711','062711','062711','062711','062711',...
            '072011','072011','072011','072011','072011','072011','072011'...
            };
        unitNums = {1.03,2.01,2.02,2.07,3.03,3.04,...
            1.01,1.02,1.03,1.04,2.01,2.03,...
            1.02,1.03,1.05,1.06,1.07,1.08,1.09...
            };
        colors = {'k','k','k','k','k','k',...
            'k','k','k','k','k','k',...
            'k','k','k','k','k','k','k'
            };
    case {'impaired'}
        dates = {'062311','062311','062311','062311',...
            '072111','072111','072111','072111','072111','072111',...
            '080111','080111','080111','080111','080111',...
            '080911','080911','080911','080911','080911','080911','080911','080911'...
            };
        unitNums = {1.01,1.09,1.11,1.14,...
            1.01, 1.06, 1.11, 1.12, 1.15, 1.16,...
            2.02, 4.01, 4.02, 4.03, 4.05,...
            1.04, 1.05, 1.06, 1.07, 1.09, 1.11, 1.15, 1.16
            };
        colors = {'r','r','r','r',...
            'r','r','r','r','r','r',...
            'r','r','r','r','r',...
            'r','r','r','r','r','r','r','r'
            };
    case {'amplified'}
        dates = {'050412','050412','050412','050412','050412','050412'...
            };
        unitNums = {1.01, 1.02, 1.03, 1.05, 2.01, 2.03...
            };
        colors = {'m','m','m','m','m','m'...
            };
end

binWidth_sec = 100e-6;
PSTdur_sec = 10.0;
numPSTsamples = round(PSTdur_sec/binWidth_sec);
PSTH = NaN*ones(10*length(unitNums),numPSTsamples);

j=0;
for i=1:length(unitNums)
    ExpDate = dates{i};
    UnitName = num2str(unitNums{i});
    
    if strcmp(ExpDate,'011811')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_01_18-AN-norm-exposed') ''''])
    elseif strcmp(ExpDate,'041811')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_04_18-AN_Normal') ''''])
    elseif strcmp(ExpDate,'051211')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_05_12-AN_Normal_Chin1120') ''''])
    elseif strcmp(ExpDate,'062711')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_06_27-Chin1133_AN_Normal') ''''])
    elseif strcmp(ExpDate,'072011')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_07_20-Chin1132_AN_normal') ''''])
    elseif strcmp(ExpDate,'062311')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_06_23-Chin1119_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'072111')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_07_21-Chin1124_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'080111')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_08_01-Chin1125_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'080911')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_08_09-Chin1135_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'081511')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_08_15-Chin1136_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'100611')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_10_06-Chin1144_AN_normal') ''''])
    elseif strcmp(ExpDate,'101111')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_10_11-Chin1139_AN_normal') ''''])
    elseif strcmp(ExpDate,'101711')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_10_17-Chin1151_AN_normal') ''''])
    elseif strcmp(ExpDate,'112911')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_11_29-Chin1146_AN_normal') ''''])
    elseif strcmp(ExpDate,'120511')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2011_12_05-Chin1148_AN_normal') ''''])
    elseif strcmp(ExpDate,'050412')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2012_05_04-Chin1149_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'062312')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2012_06_23-Chin1202_AN_500OBN') ''''])
    elseif strcmp(ExpDate,'072112')
        eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'JB-2012_07_21-Chin1206_AN_500OBN') ''''])
    else
        error('BAD dir')
    end
    
    [p,ExpName,e]=fileparts(pwd);
    data_dir=fullfile(ROOT_dir,'ExpData',ExpName);
    unitdata_dir=fullfile(data_dir,'UNITSdata');
    data_dir_bak = data_dir;
    
    %%%% Parse out the Track and Unit Number
    TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
    UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
    
    SAVECALCSfilename=sprintf('UnitLook_EHIN.%d.%02d.mat',TrackNum,UnitNum);
    eval(['load ''' fullfile(unitdata_dir,SAVECALCSfilename) ''' -regexp ^(?!data_dir_bak$|ROOT_dir$|i$).']);
    
    ATTind = find(yTEMP.Nattens_dB==max(yTEMP.Nattens_dB)); %quiet
    
%     [~,BFind]=min(abs(yTEMP.BFs_kHz-unit.Info.BF_kHz)); % center BF
    for BFind=1:length(yTEMP.BFs_kHz)
        j=j+1;
        PIC=concatPICS_NOHR(yTEMP.picNums{ATTind,BFind},yTEMP.excludeLines{ATTind,BFind});
        PIC=simFF_PICshift(PIC);
        PIC=calcPST(PIC,binWidth_sec);
        PSTH(j,:)=[PIC.PST.pst_Y_sps zeros(1,numPSTsamples-length(PIC.PST.pst_Y_sps))]; %pad to PSTdur_sec
        PSTH_BF_kHz(j) = yTEMP.BFs_kHz(BFind);%unit.Info.BF_kHz;
    end
    
    data_dir = data_dir_bak;  % restore in case it was overwritten by mat file
    cd(data_dir);
end
%
figure(1), clf; 
for i=1:length(PSTH_BF_kHz)
    hold on;
    sizeFactor = 50/max(max(PSTH(i,:).^2));
    thresh = 0;
    PSTH(PSTH<=thresh)=NaN;
    scatter((1:numPSTsamples)*binWidth_sec,...
        PSTH_BF_kHz(i)*ones(1,numPSTsamples),sizeFactor*PSTH(i,:).^2,...
        'k');%,'MarkerFaceColor',color);
%     plot((1:numPSTsamples)*binWidth_sec,PSTH_BF_kHz(i)+sizeFactor*PSTH(i,:));
    hold off;
    
end

% Compare to line fit via click response
% Greenwood Function:
A = 163.5; k = 0.85; a = 2.1; % Chinchilla
x = 0:0.001:1; %proportion of cochlear length
F_Hz = A * (10.^(a*x) - k);
[B,IX] = sort(PSTH_BF_kHz);
dist = interp1(F_Hz,x,B*1e3);
NDfit = 0.005228*dist.^2 - 0.01203*dist + 0.008404; %from clickResponse.m

switch config % set NDoffset
    case {'normal'}
        NDoffset = 14.0e-3; %ms
    case {'impaired'}
        NDoffset = 13.5e-3; %ms
    otherwise % amplified
        NDoffset = 11.5e-3; %ms
end
hold on; plot(NDfit+NDoffset,B,'r'); hold off;



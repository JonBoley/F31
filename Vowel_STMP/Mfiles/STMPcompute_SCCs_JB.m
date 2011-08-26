function STMPcompute_SCCs_JB(ExpDate,UnitName,STIMtype,DATAtype)

global STMP_dir STMP_ExpList
global FeaturesText

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=1;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03),
%%%% Specify ExpDate if not provided
if ~exist('ExpDate','var')
    if TESTunitNUM==1
        ExpDate='111804'; UnitName='1.28'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';
    elseif TESTunitNUM==2
        ExpDate='041805'; UnitName='2.04';
    elseif TESTunitNUM==3
        ExpDate='111606'; UnitName='2.09';
    elseif TESTunitNUM==4
        ExpDate='041307'; UnitName='1.03';
    end
end

%%%% Find the full Experiment Name
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(STMP_ExpList)
    if ~isempty(strfind(STMP_ExpList{i},ExpDateText))
        ExpName=STMP_ExpList{i};
        break;
    end
end
if ~exist('ExpName','var')
    disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
    disp(strvcat(STMP_ExpList))
    beep
    error('STOPPED');
end

%%%% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Setup related directories
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)


global DataList
if isempty(DataList)
    if strcmp(ExpDate,'041805')
        %%% Load DataList for 041805
        disp(' ... Loading DataList for 04_18-05');
        load DataList_2005_04_18
    elseif strcmp(ExpDate,'111804')
        %%% Load DataList for 111804
        disp(' ... Loading DataList for 11_18-04');
        load DataList_2004_11_18
    elseif strcmp(ExpDate,'071305')
        %%% Load DataList for 071305
        disp(' ... Loading DataList for 07_13-05');
        load DataList_2005_07_13
    elseif strcmp(ExpDate,'032805')
        %%% Load DataList for 032805
        disp(' ... Loading DataList for 03_28-05');
        load DataList_2005_03_28
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Verify that there is data for this unit
if isempty(DataList.Units{TrackNum,UnitNum})
    error('NO DATA FOR THIS UNIT!!');
end
Datafields=fieldnames(DataList.Units{TrackNum,UnitNum});
if ~sum(strcmp(Datafields,{'EHrlv'}))
    error('NO ''EHrlv'' DATA FOR THIS UNIT, THEREFORE STOPPING!!');
end

cd(data_dir)
disp(sprintf('Looking at Basic EH for Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile(unitdata_dir,UnitFileName) ''');'])
% If UNITSdata file does not exist, load create from DataList, ow/ load it
if isempty(ddd)
    unit=DataList.Units{TrackNum,UnitNum};
    eval(['save ''' fullfile(unitdata_dir,UnitFileName) ''' unit'])
else
    eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
end

% If EHINvNreBFi_simFF analysis is not completed, run here (this will also
% check EHINvNreBFi analysis and run if needed!)
if ~isfield(unit,'EH_reBF_simFF')
    UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);
    eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
end


numSCCs=0; %later set to max #of SCCs
% SCC_octOFFSET=0.25; parameter of the offset between BF and other AN fiber

% ASSUME OFFSET1 is below OFFSET2
SCC_octOFFSET1=NaN; %previously -0.05; % parameter of the offset between lower BF and Nominal BF
SCC_octOFFSET2=NaN; %previously 0.05; % parameter of the offset between higher BF and Nominal BF

if 	(~isnan(SCC_octOFFSET1) && ~isnan(SCC_octOFFSET2)) && SCC_octOFFSET1>SCC_octOFFSET2
    error('SCC OFFSET 1 > OFFSET2')
end
disp(sprintf('   *** SCC_octOFFSET1 = %.3f',SCC_octOFFSET1))
disp(sprintf('   *** SCC_octOFFSET2 = %.3f',SCC_octOFFSET2))

clear paramsIN
paramsIN.SCC.StartTime_sec=.02;  % Take 20-400(scaled to higher BF TimeFact) ms as stimulus window
paramsIN.SCC.EndTime_sec=.400; % scaled for each BF_Hz later
%	disp(sprintf('   * SCC window: %.3f - %.3f sec',paramsIN.SCC.StartTime_sec,paramsIN.SCC.EndTime_sec))

paramsIN.SCC.DELAYbinwidth_sec=50e-6;  % 50e-6 is what Joris used
%paramsIN.SCC.Duration_sec=paramsIN.SCC.EndTime_sec-paramsIN.SCC.StartTime_sec;

NSCCs=cell(NUMrows,length(levels_dBSPL));
NSCCs_sps=cell(NUMrows,length(levels_dBSPL));
NSCC_delays_usec=cell(NUMrows,length(levels_dBSPL));
NSCC_BFs_kHz=cell(NUMrows,length(levels_dBSPL));
NSCC_avgrates=cell(NUMrows,length(levels_dBSPL));
NSCC_nsps=cell(NUMrows,length(levels_dBSPL));
if ~exist('NSCC_CDs_usec')
    NSCC_CDs_usec=cell(NUMrows,length(levels_dBSPL));
end
if ~exist('NSCC_peaks')
    NSCC_peaks=cell(NUMrows,length(levels_dBSPL));
end
NSCC_0delay=cell(NUMrows,length(levels_dBSPL));
NSCC_Rho=cell(NUMrows,length(levels_dBSPL));
NSCC_ARBdelay=cell(NUMrows,length(levels_dBSPL));
SCCsMAX=0;
NSACs=cell(NUMrows,length(levels_dBSPL));
NSACs_sps=cell(NUMrows,length(levels_dBSPL));
NSAC_delays_usec=cell(NUMrows,length(levels_dBSPL));
NSAC_BFs_kHz=cell(NUMrows,length(levels_dBSPL));
NSAC_avgrates=cell(NUMrows,length(levels_dBSPL));
NSAC_nsps=cell(NUMrows,length(levels_dBSPL));
NSAC_CDs_usec=cell(NUMrows,length(levels_dBSPL));
NSAC_peaks=cell(NUMrows,length(levels_dBSPL));


%%%%%%%%%%%%%%%%
% Compute SCCs
%%%%%%%%%%%%%%%%
NSCCs{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCCs_sps{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_delays_usec{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_avgrates{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_nsps{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_BFs_kHz{ROWind,ATTind}=cell(size(NSCC_BFinds));
if isempty(NSCC_CDs_usec{ROWind,ATTind})
    NSCC_CDs_usec{ROWind,ATTind}=cell(size(NSCC_BFinds));
end
if isempty(NSCC_peaks{ROWind,ATTind})
    NSCC_peaks{ROWind,ATTind}=cell(size(NSCC_BFinds));
end
NSCC_0delay{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_Rho{ROWind,ATTind}=cell(size(NSCC_BFinds));
NSCC_ARBdelay{ROWind,ATTind}=cell(size(NSCC_BFinds));

for SCCind=1:length(NSCC_BFinds)  % index of SCC to calculate
    % Find SpikeTrains needed for this SCC
    emptySCC=0;
    for i=1:2
        SpikeTrains{i}=SCC_allSpikeTrains{NSCC_BFinds{SCCind}(i)};
        if isempty(SpikeTrains{i})
            emptySCC=1;
        end
    end
    if ~emptySCC
        disp(sprintf('Feature: %s; SNR: %.f dB  --  Computing SCC # %d between BFs %d and %d ........', ...
            FeaturesText{FeatIND},levels_dBSPL(ATTind),SCCind,NSCC_BFinds{SCCind}))
        FreqFact =yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(2))/unit.Info.BF_kHz;
        myDuration_sec = paramsIN.SCC.EndTime_sec/FreqFact-paramsIN.SCC.StartTime_sec;
        [NSCCs{ROWind,ATTind}{SCCind},NSCC_delays_usec{ROWind,ATTind}{SCCind},NSCC_avgrates{ROWind,ATTind}{SCCind},NSCC_nsps{ROWind,ATTind}{SCCind}] ...
            = ShufCrossCorr(SpikeTrains,paramsIN.SCC.DELAYbinwidth_sec,myDuration_sec);
        
        
        NSCCbinwidth_sec=diff(NSCC_delays_usec{ROWind,ATTind}{SCCind}(1:2))*1e-6;
        
        
        % Determine Maximum of ALL plotted SCCs (i.e., post-Smoothing)
        if max(NSCCs{ROWind,ATTind}{SCCind})>SCCsMAX
            SCCsMAX=max(NSCCs{ROWind,ATTind}{SCCind});
        end
        NSCC_BFs_kHz{ROWind,ATTind}{SCCind}=[yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(1)),yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(2))];
        NSCC_0delay{ROWind,ATTind}{SCCind}=NSCCs{ROWind,ATTind}{SCCind}(find(NSCC_delays_usec{ROWind,ATTind}{SCCind}==0));
        
        F0per_us=1/yTEMP.FeatureFreqs_Hz{1}(1)*1e6;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Manually select characteristic delay
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        loadDELAYS=0;
        if exist(fullfile(unitdata_dir,SAVECALCSfilename),'file')
            loadDELAYS=input('Do you want to load the existing characteristic delays (0: no; [1]: yes)?? ');
            if isempty(loadDELAYS) % if user just pressed [enter]
                loadDELAYS=1;
            end
            if loadDELAYS
                disp(sprintf('Loading characteristic delays from ''%s''!!',SAVECALCSfilename))
                thisPRINTyes=PRINTyes;
                eval(['load ''' fullfile(unitdata_dir,SAVECALCSfilename) ''' NSCC_CDs_usec NSCC_peaks'])
                PRINTyes=thisPRINTyes;
            end
        else
            [NSCC_CDs_usec{ROWind,ATTind}{SCCind},NSCC_peaks{ROWind,ATTind}{SCCind}] =...
                calcCD_manual(NSCCs{ROWind,ATTind}{SCCind},NSCC_delays_usec{ROWind,ATTind}{SCCind},F0per_us,NSCC_CDs_usec{ROWind,ATTind},NSCC_BFs_kHz{ROWind,ATTind},SCCind);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%
        % calculate rho
        %%%%%%%%%%%%%%%%%%%%%
        if ~isempty(NSCC_peaks{ROWind,ATTind}{SCCind})&&...
                ~isempty(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(1)})&&...
                ~isempty(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(2)})
            SCCpeak=max(NSCC_peaks{ROWind,ATTind}{SCCind});
            SACpeak1=max(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(1)});
            SACpeak2=max(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(2)});
            NSCC_Rho{ROWind,ATTind}{SCCind}=SCCpeak/sqrt(SACpeak1*SACpeak2);
        else
            fprintf('Cannot compute NSCC_Rho. No peaks in Feat %s, SNR %.f dB...\n',...
                FeaturesText{FeatIND},levels_dBSPL(ATTind))
            if isempty(NSCC_peaks{ROWind,ATTind}{SCCind})
                fprintf('\tSCC #%d (BFs %d, %d)\n',SCCind,NSCC_BFinds{SCCind})
            end
            if isempty(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(1)})
                fprintf('\tSAC BF %d\n',NSCC_BFs_kHz{ROWind,ATTind}{SCCind}(1))
            end
            if isempty(NSAC_peaks{ROWind,ATTind}{NSCC_BFinds{SCCind}(2)})
                fprintf('\tSAC BF %d\n',NSCC_BFs_kHz{ROWind,ATTind}{SCCind}(2))
            end
            NSCC_Rho{ROWind,ATTind}{SCCind}=NaN;
        end
    else  % No Data to compute SCCs
        F0per_us=1/yTEMP.FeatureFreqs_Hz{1}(1)*1e6;
        
        fprintf('Feature: %s; SNR: %.f dB  --  CAN''T compute SCC # %d between BFs %d and %d (MISSING DATA) ........', ...
            FeaturesText{FeatIND},levels_dBSPL(ATTind),SCCind,NSCC_BFinds{SCCind})
        NSCCs{ROWind,ATTind}{SCCind}=NaN*ones(1,3);
        NSCC_delays_usec{ROWind,ATTind}{SCCind}=[-F0per_us 0 F0per_us];
        NSCC_avgrates{ROWind,ATTind}{SCCind}=NaN;
        NSCC_nsps{ROWind,ATTind}{SCCind}=0;
        NSCC_BFs_kHz{ROWind,ATTind}{SCCind}=[yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(1)),yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(2))];
        NSCC_0delay{ROWind,ATTind}{SCCind}=NaN;
        NSCC_ARBdelay{ROWind,ATTind}{SCCind}=NaN;
        if isempty(NSCC_CDs_usec{ROWind,ATTind}{SCCind})
            NSCC_CDs_usec{ROWind,ATTind}{SCCind}=NaN;
        end
        NSCC_Rho{ROWind,ATTind}{SCCind}=NaN;
        NSCC_peaks{ROWind,ATTind}{SCCind}=NaN;
    end %if ~emptySCC
end %for SCCind=...

return; %eof


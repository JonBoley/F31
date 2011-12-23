function Compare_Experiment_vs_Model()
% Compare vowel STMP coding to model
fprintf('Comparing AN Experiment to Model Data ...\n');

ExpDate = '062711';%'120511';
UnitName = '1.01';%'2.08';

% parameters
PlotSACSCC = 1;

%% Get spikes, BF, SR, etc from experiment
ExpData = getSTMPspikes_AN(ExpDate,UnitName);
SNR = Inf;%max(ExpData.SNRs)
ATTind = 4;

%% Generate spikes from model
formants =  [500 1700 2500 3300 3700];
featureNum = 2; % assume we are centered on F2 for now

% parameters
[c,BFind] = min(abs(log2(ExpData.unitBF./ExpData.BFs_kHz))); % find nearest to unitBF
BFind2 = BFind+1; % just use the next one for now
ModelParams.SNR = SNR;
ModelParams.stimDur_on = ExpData.StmOn(BFind);
ModelParams.stimDur_off = ExpData.StmOff(BFind);
ModelParams.CFs = ExpData.BFs_kHz([BFind,BFind2]);
% ModelParams.CFs = formants(featureNum)*ExpData.BFs_kHz([BFind,BFind2])/ExpData.BFs_kHz(BFind)/1e3;
ModelParams.SR = ExpData.SR;
ModelParams.SPL = ExpData.SPLvowel;

% do 2 versions:
%   - same number of reps as AN data
%   - lots of reps
ModelParams.Nreps = max(ExpData.SpikeTrains{ATTind,BFind}(:,1));
ModelData1 = getSTMPspikes_model(ModelParams);

ModelParams.Nreps = 120; % this should be plenty
ModelData2 = getSTMPspikes_model(ModelParams);

%% Analysis

% calculate period histogram
PERhistParams.spikeTimes = ExpData.SpikeTrains{ATTind,BFind}; %just one CF
PERhistParams.StimDur = ExpData.StmOn(BFind);
PERhistParams.F0 = ExpData.F0s(BFind);
PERhistParams.TimeFact = ExpData.TimeFact(BFind);
PERhistParams.NumLines = ExpData.NumLines(BFind);
PerHist_exp = calcPERhist(PERhistParams);

PERhistParams.F0 = 100;
PERhistParams.TimeFact = 1;
PERhistParams.StimDur = ModelParams.stimDur_on;
PERhistParams.spikeTimes = ModelData1.SpikeTrains{1}; %just one CF
PERhistParams.NumLines = max(PERhistParams.spikeTimes(:,1));
PerHist_model1 = calcPERhist(PERhistParams);

PERhistParams.spikeTimes = ModelData2.SpikeTrains{1}; %just one CF
PERhistParams.NumLines = max(PERhistParams.spikeTimes(:,1));
PerHist_model2 = calcPERhist(PERhistParams);

% sync (to F0 and each harmonic)
[PerHistSync_exp Rate_exp]= CalcSync(PerHist_exp, ExpData.F0s(BFind));
[PerHistSync_model1 Rate_model1] = CalcSync(PerHist_model1, ExpData.F0s(BFind));
[PerHistSync_model2 Rate_model2] = CalcSync(PerHist_model2, ExpData.F0s(BFind));

% SAC/SCC (use CCCanal)
paramsIN.SACbinwidth_msec = 50e-6;
paramsIN.ignoreONSET_msec = 0.020;
paramsIN.durA_msec = 2.500e3; % 2.5sec for long duration stimuli
paramsIN.durB_msec = paramsIN.durA_msec;
paramsIN.CF_A_Hz = ExpData.BFs_kHz(BFind)*1000;
paramsIN.CF_B_Hz = ExpData.BFs_kHz(BFind2)*1000;


% Need to format spike trains, each in [cell_array{Nreps}] format
%       A+: SpikeTrains{1,1}
%       A-: SpikeTrains{1,2}
%       B+: SpikeTrains{2,1}
%       B-: SpikeTrains{2,2}
SpikeTrains_exp{1,1} = getDrivenSpikeTrains(ExpData.SpikeTrains{ATTind,BFind});
SpikeTrains_exp{1,2} = [];
SpikeTrains_exp{2,1} = getDrivenSpikeTrains(ExpData.SpikeTrains{ATTind,BFind2});
SpikeTrains_exp{2,2} = [];
fprintf('   Calculating SACs/SCCs for AN data ...\n');
[SACSCCfunctions_exp,SACSCCmetrics_exp,paramsOUT_exp] = ...
    CCCanal_6(SpikeTrains_exp,paramsIN,0);
% Note: last cell is avg of previous bootstraps

SpikeTrains_model1{1,1} = getDrivenSpikeTrains(ModelData1.SpikeTrains{1,1});
SpikeTrains_model1{1,2} = [];
SpikeTrains_model1{2,1} = getDrivenSpikeTrains(ModelData1.SpikeTrains{1,2});
SpikeTrains_model1{2,2} = [];
fprintf('   Calculating SACs/SCCs for model data (1 of 2) ...\n');
[SACSCCfunctions_model1,SACSCCmetrics_model1,paramsOUT_model1] = ...
    CCCanal_6(SpikeTrains_model1,paramsIN,0);
% Note: last cell is avg of previous bootstraps

SpikeTrains_model2{1,1} = getDrivenSpikeTrains(ModelData2.SpikeTrains{1,1});
SpikeTrains_model2{1,2} = [];
SpikeTrains_model2{2,1} = getDrivenSpikeTrains(ModelData2.SpikeTrains{1,2});
SpikeTrains_model2{2,2} = [];
fprintf('   Calculating SACs/SCCs for model data (2 of 2) ...\n');
[SACSCCfunctions_model2,SACSCCmetrics_model2,paramsOUT_model2] = ...
    CCCanal_6(SpikeTrains_model2,paramsIN,0);
% Note: last cell is avg of previous bootstraps

if PlotSACSCC
    figure(111), subplot(311), hold on;
    plot(SACSCCfunctions_exp{end}.delays_usec,SACSCCfunctions_exp{end}.SAC_A_avg, 'g');
    plot(SACSCCfunctions_exp{end}.delays_usec,SACSCCfunctions_exp{end}.SAC_B_avg, 'r');
    plot(SACSCCfunctions_exp{end}.delays_usec,SACSCCfunctions_exp{end}.SCC_AB_avg, 'k');
    hold off;
    legend('SAC_A','SAC_B','SCC_A_B');
    xlabel('Delay (\mus)'); title('SAC & SCC (AN Data)');
    
    subplot(312), hold on;
    plot(SACSCCfunctions_model1{end}.delays_usec,SACSCCfunctions_model1{end}.SAC_A_avg, 'g');
    plot(SACSCCfunctions_model1{end}.delays_usec,SACSCCfunctions_model1{end}.SAC_B_avg, 'r');
    plot(SACSCCfunctions_model1{end}.delays_usec,SACSCCfunctions_model1{end}.SCC_AB_avg, 'k');
    hold off;
    title('SAC & SCC (Model, limited reps)');
    
    subplot(313), hold on;
    plot(SACSCCfunctions_model2{end}.delays_usec,SACSCCfunctions_model2{end}.SAC_A_avg, 'g');
    plot(SACSCCfunctions_model2{end}.delays_usec,SACSCCfunctions_model2{end}.SAC_B_avg, 'r');
    plot(SACSCCfunctions_model2{end}.delays_usec,SACSCCfunctions_model2{end}.SCC_AB_avg, 'k');
    hold off;
    title('SAC & SCC (Model, max reps)');
end

%% plots to compare AN & model data
plotData.CFs_kHz = ModelParams.CFs;
plotData.PerHists{1} = PerHist_exp;
plotData.PerHists{2} = PerHist_model1;
plotData.PerHists{3} = PerHist_model2;
plotData.Rates{1} = Rate_exp;
plotData.Rates{2} = Rate_model1;
plotData.Rates{3} = Rate_model2;
plotData.PerHistSyncs{1} = PerHistSync_exp;
plotData.PerHistSyncs{2} = PerHistSync_model1;
plotData.PerHistSyncs{3} = PerHistSync_model2;
plotData.SCC_Rhos{1} = max(SACSCCmetrics_exp{end}.SCCpeak_AB);
plotData.SCC_Rhos{2} = max(SACSCCmetrics_model1{end}.SCCpeak_AB);
plotData.SCC_Rhos{3} = max(SACSCCmetrics_model2{end}.SCCpeak_AB);
plotData.SCC_CDs{1} = abs(mod(SACSCCmetrics_exp{end}.CDscc_usec,1e6/100)); %wrap to within one period
plotData.SCC_CDs{2} = abs(mod(SACSCCmetrics_model1{end}.CDscc_usec,1e6/100));
plotData.SCC_CDs{3} = abs(mod(SACSCCmetrics_model2{end}.CDscc_usec,1e6/100));
plotSTMPcomparison(plotData);


function EXPdata = getSTMPspikes_AN(ExpDate,UnitName)
setup_Vowel_STMP;

DesiredFeature = 'F2'; %hardcode for now

% cdd
if strcmp(ExpDate,'041805')
    eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'MH-2005_04_18-ANnorm') ''''])
elseif strcmp(ExpDate,'111804')
    eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'MH-2004_11_18-NOHRnorm') ''''])
elseif strcmp(ExpDate,'071305')
    eval(['cd ''' fullfile(ROOT_dir,filesep,'ExpData',filesep,'MH-2005_07_13-ANdeafcat') ''''])
    %%%%%%%%%%%%%%%%
    % Normal: 011811, 041811, 051211, 062711, 072011
    % Impaired: 062311, 072111, 080111, 080911, 081511
elseif strcmp(ExpDate,'011811')
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
else
    error('BAD dir')
end

% Set directory names
[p,ExpName,e,v]=fileparts(pwd);
data_dir=fullfile(ROOT_dir,'ExpData',ExpName);
unitdata_dir=fullfile(data_dir,'UNITSdata');

% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));

% load data
SAVECALCSfilename=sprintf('UnitLook_EHIN.%d.%02d.mat',TrackNum,UnitNum);
data_dir_bak = data_dir;  % backup in case it is overwritten by mat file
eval(['load ''' fullfile(unitdata_dir,SAVECALCSfilename) '''']);
data_dir = data_dir_bak;  % restore in case it was overwritten by mat file
cd(data_dir)

% Plot tuning curve
% UnitVerify_TC(ExpDate,UnitName,1); drawnow;

% Get pictures
EHvNreBFi_piclist=findPics('EHvLTASSrBFi',[TrackNum UnitNum]);

HarmonicsIND=2;
PolarityIND=1;
Nattens_dB=[]; BFsTEMP_kHz=[];
if isfield(unit,'EHvLTASS_reBF_simFF')
    EHfeats=fieldnames(unit.EHvN_reBF);
    inds=find(~strcmp(EHfeats,'interleaved'))';
    for FeatIND=inds   % eliminate 'interleaved' field
        FeatINDs(FeatIND)=find(strcmp(FeaturesText,EHfeats{FeatIND}));
    end
    FeatINDs=FeatINDs(find(FeatINDs~=0));

    for FeatIND=FeatINDs
        if ~isempty(unit.EHvLTASS_reBF_simFF.F2{HarmonicsIND,PolarityIND})
            eval(['yTEMP=unit.EHvLTASS_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
        end

        if ~isempty(yTEMP)
            Nattens_dB=union(Nattens_dB,yTEMP.Nattens_dB);
            BFsTEMP_kHz=union(BFsTEMP_kHz,yTEMP.BFs_kHz);
        end %if ~isempty(yTEMP)

    end %for FeatIND=FeatINDs
end %if isfield(unit,'EHvLTASS_reBF_simFF')

%%% 
% just look at one feature
FeatINDs = find(strcmp(FeaturesText,DesiredFeature));
%%%

fprintf('   Getting spikes from experiment %s, unit %s ...\n',ExpDate,UnitName);
SpikeTrains = cell(length(Nattens_dB),length(BFsTEMP_kHz));
for ATTEN=Nattens_dB
    if isfield(unit,'EHvLTASS_reBF_simFF')
        for FeatIND=FeatINDs
            eval(['yTEMP=unit.EHvLTASS_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};']);
            if ~isempty(yTEMP)
                ATTind=find(yTEMP.Nattens_dB==ATTEN);
                EXPdata.SNRs(ATTind)= Nattens_dB(ATTind) - Nattens_dB(yTEMP.EqualSPL_index);
                for BFind=1:length(yTEMP.BFs_kHz)
                    EXPdata.BFs_kHz(BFind) = yTEMP.BFs_kHz(BFind);
                    EXPdata.F0s(BFind) = yTEMP.FeatureFreqs_Hz{BFind}(1);

                    if ~isempty(yTEMP.picNums{ATTind,BFind})
                        PIC=concatPICS_NOHR(yTEMP.picNums{ATTind,BFind},yTEMP.excludeLines{ATTind,BFind});

                        % need to correct for incorrect hardware trigger times
                        baselineStmOn = PIC.x.Hardware.Trigger.StmOn/1000*yTEMP.TimeFact{1,1};
                        baselineStmOff = PIC.x.Hardware.Trigger.StmOff/1000;
                        
                        EXPdata.StmOn(BFind) = baselineStmOn/yTEMP.TimeFact{ATTind,BFind};
                        EXPdata.StmOff(BFind) = baselineStmOff/yTEMP.TimeFact{ATTind,BFind};
                        EXPdata.NumLines(BFind) = PIC.x.Stimuli.fully_presented_lines;

                        % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                        PIC=simFF_PICshift(PIC);
                        EXPdata.TimeFact(BFind) = PIC.simFF_PICshift.TimeFact;
                    end %if ~isempty(yTEMP.picNums{ATTind,BFind})

                    %%%%%%%%%%%%%%%%
                    SpikeTrains{ATTind,BFind}=PIC.x.spikes{1};
                    %%%%%%%%%%%%%%%%

                    clear PIC;
                end %for BFind=1:length(yTEMP.BFs_kHz)

            end %if ~isempty(yTEMP)

        end %for FeatIND=FeatINDs

    else
        warning('EHvLTASS_reBF_simFF is not a valid field!');
    end %if isfield(unit,'EHvLTASS_reBF_simFF')
end %for ATTEN=Nattens_dB

EXPdata.SpikeTrains = SpikeTrains;
EXPdata.Nattens_dB = Nattens_dB;
EXPdata.unitBF = unit.Info.BF_kHz;
EXPdata.SR = unit.Info.SR_sps;
EXPdata.TC = unit.TC{1};
EXPdata.SPLvowel = unit.EHvLTASS_reBF_simFF.F2{2,1}.levels_dBSPL;

%%


function ModelData = getSTMPspikes_model(ModelParams)
addpath(genpath('C:\Research\MATLAB\Vowel_STMP\Model'));

% hard-coded params (temporary)
impaired = 0; % yes/no
amplified = 0; % [0 1 2] = [off linear nonlinear]
featureNum = 2; % [1 2 3 ...] = [F1 F2 F3 ...]
SNRs = Inf;%[Inf 6 0 -6]; % Is the noise actually getting loud enough in the speech band?

% get model params from struct
Levels = ModelParams.SPL;
LevelIndex = 1; %just do one presentation level
OALevel_dBSPL = Levels(LevelIndex);

% Generate vowel ('eh')
dur = ModelParams.stimDur_on;
F0=100;
Fs=24414.062500;
formants =  [500 1700 2500 3300 3700];
%%% center F2 on CF
formants = formants*ModelParams.CFs(1)*1000/1700;
%%%
BWs= [200 200 200 250 200];
FeaturesText={'F1','F2','F3','F4','F5'};
[time, vowel] = dovowl(formants,BWs,F0,dur,Fs);
vowel=vowel./max(abs(vowel))*0.99; % normalize

%% Initialize model parameters
Nreps = ModelParams.Nreps; %number of repetitions
CF_kHz = ModelParams.CFs;
numCFs = length(CF_kHz);

if ModelParams.SR < 1
    fibertype=1;    %[1,2,3]=[low,med,high] spontaneous rate
elseif ModelParams.SR < 18
    fibertype=2;
else
    fibertype=3;
end
Cohc=1.0; %outer hair cell health
Cihc=1.0; %inner hair cell health
ANmodel_Fs_Hz=100e3; %100kHz sampling

% Initialize analysis variables
analysis_init; % default settings (paste the script content here before modifying)

% Add noise
SNRindex = 1; % just do one noise level for now
LTASS=GenLTASS; %generate once, use across conditions
[signal,delta_dB]=AddNoise(vowel,LTASS,ModelParams.SNR);
fprintf('   Added LTASS noise at %ddB SNR (increased total by %ddB) ...\n',ModelParams.SNR,delta_dB);

% apply hearing aid amplification here
% [see modeling code at C:\Research\MATLAB\Vowel_STMP\Model\STMP_Vowel.m]

% Adjust the sample rate & level for the model
% [for STMP, adjust Fs first (?)]
refit_stim;

% Adjust stimulus length and apply window
signal_model = refit_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);
signal_model = window_waveform(signal_model,ANmodel_Fs_Hz,dur_sec*1000);

% Set impairment here
% [see modeling code at C:\Research\MATLAB\Vowel_STMP\Model\STMP_Vowel.m]

%% Run model
fprintf('   Running model (%d reps) ...\n',ModelParams.Nreps);
for FiberNumber=1:numCFs
    vihc = catmodel_IHC(signal_model.',CF_kHz(FiberNumber)*1e3,1,...
        1/ANmodel_Fs_Hz,dur_sec+ModelParams.stimDur_off,Cohc,Cihc);
    [sout,psth]=catmodel_Synapse(vihc,CF_kHz(FiberNumber)*1e3,1,...
        1/ANmodel_Fs_Hz,dur_sec+ModelParams.stimDur_off,fibertype,1);
    SynOut{FiberNumber,LevelIndex,SNRindex}=sout; % save the synapse output
    [sptimes nspikes] = SGfast([1/ANmodel_Fs_Hz, Nreps], sout);
    NELspikes=ANmodelSTs2nel(sptimes,Nreps); % convert to NEL formatting
    SpikeTrains_plus = nelSTs2cell(NELspikes);
    if length(SpikeTrains_plus)<Nreps
        SpikeTrains_plus = [SpikeTrains_plus cell(1,Nreps-length(SpikeTrains_plus))];
    end
    
    ModelData.SpikeTrains{LevelIndex,FiberNumber} = [];
    for i=1:length(SpikeTrains_plus)
        ModelData.SpikeTrains{LevelIndex,FiberNumber} = ...
            [ModelData.SpikeTrains{LevelIndex,FiberNumber};...
            i*ones(length(SpikeTrains_plus{i}),1),SpikeTrains_plus{i}];
    end
    
end

%%


function PERhistStruct = calcPERhist(PERhistParams)
spikeTimes = PERhistParams.spikeTimes;
StimDur = PERhistParams.StimDur;
F0 = PERhistParams.F0;
TimeFact = PERhistParams.TimeFact;
NumLines = PERhistParams.NumLines;

PERhist_window_sec=[0.020 StimDur]; % matches Wong et al 1998, Miller et al 1999a,b

%%%% PERIOD histogram usage for getting Synch/Phase %%%%%%%%%%%
% K=16 used for all frequencies: estimates 1st 7 harmonics, which is not enough for vowels.
% Johnson (1980) used 30-200 bins/cycle for the PERhist;
% Anderson et al (1971) used 10 bins/cycle for PERhist.
% With K bins/cycle for all frequencies, Synch/Phase from PERhist DFT gives values at harmonics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=256;  %% # bins/cycle
binWidth_sec=1/F0/K;
M=floor(diff(PERhist_window_sec)*F0);  % Integer number of cycles to include in driven-spike window
PERhist_window_sec(2)=PERhist_window_sec(1)+M/F0; % Reset EndTime to limit to integer number of cycles of F0

% Find driven spikes to use
drivenSpikeIndices = find( (spikeTimes(:,2) >= PERhist_window_sec(1)) & (spikeTimes(:,2) <= PERhist_window_sec(2)) );
drivenSpikes=spikeTimes(drivenSpikeIndices,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Period Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivenSpikes_BINS=rem(floor(drivenSpikes/binWidth_sec),K)+1; %Convert times to PERhist bins (1:K)
[PERhist,xxx]=hist(drivenSpikes_BINS,(1:K));  % Make Histogram from BINS (1:K)
% This is the actual number of recorded spikes used to create this PERhist
% (use this for stats, etc ...)
NumDrivenSpikes=sum(PERhist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADJUST for temporal scaling bias (if this is a simFF condition)!!!
PERhist = PERhist * TimeFact;
SCALED_NumDrivenSpikes=sum(PERhist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Convert PERhist to spikes per second
PERhist=PERhist/NumLines/M/binWidth_sec; % Convert to sp/sec

% store in a struct
PERhistStruct.PERhist = PERhist;
PERhistStruct.NumDrivenSpikes = SCALED_NumDrivenSpikes;
PERhistStruct.binWidth_sec = binWidth_sec;

%%


function [PerHistSync Rate] = CalcSync(PERhist,F0)

%%%% Calculate DFT from PERIOD histogram - to get Synch and phase
SynchRate_PERhist=fft(PERhist.PERhist)/length(PERhist.PERhist);
FFTfreqs=(0:length(SynchRate_PERhist)-1)*(1/PERhist.binWidth_sec)/length(SynchRate_PERhist);

Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
RayleighCRIT=chi2inv(1-Rayleigh_P,2);

% calculate synch to each harmonic
HarmonicFreqs_Hz = F0*(1:floor(5e3/F0)); % up to 5kHz
HarmonicINDs=NaN+zeros(size(HarmonicFreqs_Hz));
HarmonicSynchs=NaN+zeros(size(HarmonicFreqs_Hz));
HarmonicPhases=NaN+zeros(size(HarmonicFreqs_Hz));
HarmonicRaySig=zeros(size(HarmonicFreqs_Hz));

for i=1:length(HarmonicFreqs_Hz)
    %%%% Because the FFTfreqs are not exactly equal to the harmonics, we need to pick the closest one
    [yyy,HARMind]=min(abs(FFTfreqs-HarmonicFreqs_Hz(i)));
    if ~isempty(HARMind)
        HarmonicINDs(i)=HARMind;
        HarmonicSynchs(i)=abs(SynchRate_PERhist(HarmonicINDs(i)))/SynchRate_PERhist(1);
        HarmonicPhases(i)=angle(SynchRate_PERhist(HarmonicINDs(i)));
        % NOTE: for sim_FF conditions, the SCALED PERhist is used, but the
        % original number of spikes is used for statistics!!
        RayleighStat=2*PERhist.NumDrivenSpikes*HarmonicSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
        if RayleighStat>RayleighCRIT
            HarmonicRaySig(i)=1;
        end
    end
end

% save as a struct
PerHistSync.SynchRate_PERhist = SynchRate_PERhist;
PerHistSync.HarmonicFreqs_Hz = HarmonicFreqs_Hz;
PerHistSync.HarmonicSynchs = HarmonicSynchs;
PerHistSync.HarmonicPhases = HarmonicPhases;
PerHistSync.HarmonicRaySig = HarmonicRaySig;
Rate = SynchRate_PERhist(1); % DC component

% figure,
% bar(PerHistSync.HarmonicFreqs_Hz,PerHistSync.HarmonicSynchs,'b');
% set(gca,'XScale','log')
% xlim([100 5e3]); ylim([0 1]);
% set(gca,'XTick',[100 200 500 1000 2000 5000]);
% xlabel('Frequency (Hz)'); ylabel('Sync Coef');
% title('Sync to individual harmonics');
%%


function plotSTMPcomparison(plotData)
NUMcols=5; % time, rate, synchrony, correlation, delay
NUMrows=1;
figure(502); clf
set(gcf,'Name','Comparison of Neural Coding Strategies');
set(gcf,'units','norm','pos',[0.1,0.1,0.8,0.8]);
ALLlevelsTriFiltTONE=9;
ALLlevelsTriFilt=3;
colors = ['b','g','r','k'];
symbols = ['o','x','*','+'];

%%%% Spatio-Temporal Plots
PLOTnum=1;
h1=subplot(NUMrows,NUMcols,PLOTnum);
for i=1:length(plotData.PerHists)
    plot(trifilt(plotData.PerHists{i}.PERhist,ALLlevelsTriFilt),...
        plotData.PerHists{i}.binWidth_sec*(1:length(plotData.PerHists{i}.PERhist))*1000,...
        colors(i))
    hold on;
end
hold off; ylabel('Time (ms)');
xlabel('Rate (spikes/sec)');
set(h1,'XDir','reverse')
ylim([0 1e3/plotData.PerHistSyncs{1}.HarmonicFreqs_Hz(1)]);
legend('AN data','Model (limited reps)','Model (max reps)');

%%%% Rate Plot
PLOTnum=2;
h2=subplot(NUMrows,NUMcols,PLOTnum);
for i=1:length(plotData.Rates)
    semilogy(plotData.Rates{i},plotData.CFs_kHz(1)*1e3,['o' colors(i)]); hold on;
end
hold off;
set(h2,'XDir','reverse'); xlabel('Rate (spikes/sec)');
set(h2,'YTick',[100 200 500 1000 2000 5000]); ylabel('CF (Hz)');
ylim([100 5e3])

%%%% Synch Plot
PLOTnum=3;
h3=subplot(NUMrows,NUMcols,PLOTnum);
for i=1:length(plotData.PerHistSyncs)
    temp = plotData.PerHistSyncs{i}.HarmonicSynchs;
    temp(~plotData.PerHistSyncs{i}.HarmonicRaySig)=NaN; %plot sig data only
    semilogy(temp,...
        plotData.PerHistSyncs{i}.HarmonicFreqs_Hz,['x-' colors(i)]);
    hold on;
end
hold off;
xlabel(sprintf('Synch Coefs'));
xlim([0 1])
set(h3,'XDir','reverse')
set(h3,'YTick',[100 200 500 1000 2000 5000])
set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
ylim([100 5e3])


%%%% Rho Plot
PLOTnum=4;
h4=subplot(NUMrows,NUMcols,PLOTnum);
symbols = {'',''};
for i=1:length(plotData.SCC_Rhos)
    semilogy(plotData.SCC_Rhos{i},plotData.CFs_kHz(1)*1e3,['o' colors(i)]); hold on;
end
hold off;
xlabel(sprintf('Corr Coefs'))
% xlim([0 5])
set(h4,'XDir','reverse')
set(h4,'YTick',[100 200 500 1000 2000 5000]); %ylabel('CF (Hz)');
ylim([100 5e3])


%%%% CD Plot
PLOTnum=5;
h5=subplot(NUMrows,NUMcols,PLOTnum);
for i=1:length(plotData.SCC_CDs)
    semilogy(plotData.SCC_CDs{i},plotData.CFs_kHz(1)*1e3,['o' colors(i)]); hold on;
end
hold off;
XLIMITS_cd = [-2 2];
% xlim(XLIMITS_cd)
set(h5,'XDir','reverse'); xlabel('CD (\mus)');
set(h5,'YTick',[100 200 500 1000 2000 5000]); %ylabel('CF (Hz)');
ylim([100 5e3])


Xcorner=0.05;
Xwidth1=.4;  Xshift1=0.05;
Xwidth2=.1;  Xshift2=0.03;

Ycorner=0.05;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

TICKlength=0.02;

set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h5,'Position',[Xcorner+Xwidth1+Xshift1+3*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

orient landscape
keyboard;
%%


%%

% eof

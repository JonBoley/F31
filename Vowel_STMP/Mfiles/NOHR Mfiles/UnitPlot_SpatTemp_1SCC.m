function UnitPlot_SpatTemp_1SCC(ExpDate,UnitName)
% File: UnitPlot_SpatTemp.m
%
% Modified Date: 07Jan2005 (M Heinz)
% Modified From:  UnitPlot_SpatTemp.m and UnitPlot_SpatTemp.m to use either Interleaved Data or regular in 1 file
% 
% Date: 01Nov2004 (M. Heinz) (Modified from UnitPlot_EHrBF_simFF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '080204' (converted later)
% UnitName: '1.29' (converted later)
%
% Plots Simulated spatio-temporal response pattern (ala Shamma 1985) from EH_reBF_simFF and/or T_reBF_simFF data 
% for a given experiment and unit.  Loads 'UNITSdata/unit.T.U.mat' file.
% UnitAnal_EHrBF_simFF.m and/or UnitAnal_TonerBF_simFF.m performs the relevant analysis.
%
% 11/2/04 TO DO
% 1) only plots Fn/Tn{1,1} for now, need to add rest later
% *2) Separate Calcs and Plots, to allow smart plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NOHR_dir NOHR_ExpList
global SavedPICS SavedPICnums SavedPICSuse
global FeaturesText FormsAtHarmonicsText InvertPolarityText
SavedPICSuse=1; SavedPICS=[]; SavedPICnums=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
   %%% HARD CODE FOR NOW
%    ExpDate='080204'
   ExpDate='111804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''080204''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
   
   %%% HARD CODE FOR NOW
   if strcmp(ExpDate,'080204')
      UnitName='1.29'
   elseif strcmp(ExpDate,'111804')
      UnitName='1.28'
   end
   
   while ~ischar(UnitName)
      UnitName=input('Enter Unit Name (e.g., ''1.29''): ');
   end
end

%%%% Find the full Experiment Name 
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(NOHR_ExpList)
   if ~isempty(strfind(NOHR_ExpList{i},ExpDateText))
      ExpName=NOHR_ExpList{i};
      break;
   end
end
if ~exist('ExpName','var')
   disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
   disp(strvcat(NOHR_ExpList))
   beep
   error('STOPPED');
end

%%%% Parse out the Track and Unit Number 
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
disp(sprintf('Plotting Spatio-Temporal Patterns for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify UNITSdata/unit analyses are all done ahead of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
% If UNITSdata file does not exist, load DataList to see if interleaved or not, and then run necessary analyses
if isempty(ddd)
   FileName=strcat('DataList_',ExpDateText,'.mat');
   disp(['   *** Loading file: "' FileName '" because UNITSdata/unit does not exist yet'])
   eval(['load ' FileName])
   
   if isfield(DataList.Units{TrackNum,UnitNum},'Tone_reBF')   
      UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   end
   if isfield(DataList.Units{TrackNum,UnitNum},'EH_reBF')   
      UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);
   end
   % Verify TC
   UnitVerify_TC(ExpDate,UnitName);   
   % Estimate SR
   UnitCalc_SR(ExpDate,UnitName);   
   
end
eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

% Make sure there is either Tone_reBF or EH_reBF data
if (~isfield(unit,'EH_reBF'))&(~isfield(unit,'Tone_reBF'))
   disp(sprintf('*** No Tone_reBF or EH_reBF data for this unit!!!'))
   beep
   return;
end

% If EH_reBF_simFF analysis is not completed, run here
if isfield(unit,'EH_reBF')&~isfield(unit,'EH_reBF_simFF')
   UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);   
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

% If Tone_reBF_simFF analysis is not completed, run here
if isfield(unit,'Tone_reBF')&~isfield(unit,'Tone_reBF_simFF')
   UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

% Make sure TC is verfified (i.e., Q10 known)
if isempty(unit.Info.Q10)
   UnitVerify_TC(ExpDate,UnitName);
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

% Make sure SR is estimated
if isempty(unit.Info.SR_sps)
   UnitCalc_SR(ExpDate,UnitName);   
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Find out if there are multiple data sets for this unit (i.e., can't be shown together)
%    and if so, ask which to show
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) CHECK Tone_reBF and EH_reBF
if isfield(unit,'Tone_reBF')&isfield(unit,'EH_reBF')
   unit
   TONElevels=unit.Tone_reBF_simFF.levels_dBSPL;
   EHfeats=fieldnames(unit.EH_reBF_simFF);
   for FeatIND=1:length(EHfeats)
      FeatINDs(FeatIND)=find(strcmp(FeaturesText,EHfeats{FeatIND}));
   end
   EHlevels=[];
   for FeatIND=FeatINDs
      for HarmonicsIND=1:2
         for PolarityIND=1:2
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               EHlevels=union(EHlevels,yTEMP.levels_dBSPL);
            end
         end
      end
   end
   % if ALL levels don't overlap, then ask which to show
   if max(size(intersect(EHlevels,TONElevels)))~=min([length(TONElevels),length(EHlevels)])
      beep
      pause(.5)
      beep
      while (1==1)
         TEMP=input(sprintf('*****\n  Tone levels are: %s dB SPL\n  EH-vowel levels are: %s dB SPL\nDo you want to show Tones [1] or EH-Vowel [2] ',num2str(TONElevels),num2str(EHlevels)));
         if TEMP==1
            unit=rmfield(unit,{'EHrlv','EH_reBF','EH_reBF_simFF'});
            break
         elseif TEMP==2
            unit=rmfield(unit,{'Tone_rlv','Tone_reBF','Tone_reBF_simFF'});
            break
         end
      end
   end
   unit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this unit, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Find number of features to plot (tone, F1, T1, ...)
%%%% Find all levels and BFs to plot
NUMcols=4;
NUMrows=0;
levels_dBSPL=[]; BFsTEMP_kHz=[];
if isfield(unit,'Tone_reBF_simFF')
   NUMrows=NUMrows+1;
   levels_dBSPL=union(levels_dBSPL,unit.Tone_reBF_simFF.levels_dBSPL);
   BFsTEMP_kHz=union(BFsTEMP_kHz,unit.Tone_reBF_simFF.BFs_kHz);
end
if isfield(unit,'EH_reBF_simFF')
   EHfeats=fieldnames(unit.EH_reBF_simFF);
   NUMrows=NUMrows+length(EHfeats);
   %    clear FeatINDs
   for FeatIND=1:length(EHfeats)
      FeatINDs(FeatIND)=find(strcmp(FeaturesText,EHfeats{FeatIND}));
   end
   
   F0min=Inf;  
   for FeatIND=FeatINDs
      for HarmonicsIND=1:2
         for PolarityIND=1:2
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               levels_dBSPL=union(levels_dBSPL,yTEMP.levels_dBSPL);
               if FeatIND==FeatINDs(1)
                  FeatureLevels_dB=yTEMP.FeatureLevels_dB;
               else   
                  if sum(FeatureLevels_dB(~isnan(FeatureLevels_dB))-yTEMP.FeatureLevels_dB(~isnan(yTEMP.FeatureLevels_dB)))
                     error('FeatureLevels_dB do not match across Features for this unit');
                  end
               end
               BFsTEMP_kHz=union(BFsTEMP_kHz,yTEMP.BFs_kHz);
               % Find minimum F0 for PERhist XMAX
               for i=1:length(yTEMP.FeatureFreqs_Hz)
                  if ~isempty(yTEMP.FeatureFreqs_Hz{i})
                     if yTEMP.FeatureFreqs_Hz{i}(1)<F0min
                        F0min=yTEMP.FeatureFreqs_Hz{i}(1);
                     end
                  end
               end
            end
         end
      end
   end
end
lowBF=min(BFsTEMP_kHz);
highBF=max(BFsTEMP_kHz);
clear BFsTEMP_kHz;
TFiltWidth=1;   % What is a good number here?? is 1 OK, ow, you get major smoothing for low F0, and not much smoothing for high F0s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL CALCS (BEFORE PLOTTING), e.g., PERhists, DFTs, SCCs, ...
%%%%   - runs through all BFs, levels and saves ALL calcs prior to PLOTTING (allows amart plotting based on ALL data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BFs_kHz=cell(NUMrows,length(levels_dBSPL));
PERhists=cell(NUMrows,length(levels_dBSPL));
PERhists_Smoothed=cell(NUMrows,length(levels_dBSPL));
PERhistXs_sec=cell(NUMrows,length(levels_dBSPL));  % for plotting
DFTs=cell(NUMrows,length(levels_dBSPL));
DFTfreqs_Hz=cell(NUMrows,length(levels_dBSPL));
Nsps=cell(NUMrows,length(levels_dBSPL));
Rates=cell(NUMrows,length(levels_dBSPL));
Synchs=cell(NUMrows,length(levels_dBSPL));
Phases=cell(NUMrows,length(levels_dBSPL));
PERhistsMAX=0;
PERhistsYCHANS=0;
DFTsMAX=0;
SMP_rate=cell(1,length(levels_dBSPL));
for LEVind=1:length(levels_dBSPL)   
   SMP_rate{LEVind}=NaN*ones(1,length(FeaturesText));
end
ALSRs=cell(NUMrows,length(levels_dBSPL));
SMP_alsr=cell(1,length(levels_dBSPL));
for LEVind=1:length(levels_dBSPL)
   SMP_alsr{LEVind}=NaN*ones(1,length(FeaturesText));
end

%%%% SCC variables and params
% Changed 16Feb2005 for ARO2005
numSCCs=1; % Right now set to BF+SCC_octOFFSET vs BF-SCC_octOFFSET octaves
% SCC_octOFFSET=0.25; parameter of the offset between BF and other AN fiber
SCC_octOFFSET=0.05; % parameter of the offset between BF and other AN fiber
disp(sprintf('   *** SCC_octOFFSET = %.3f',SCC_octOFFSET))
clear paramsIN
paramsIN.SCC.StartTime_sec=.02;  % Take 20-400 ms as stimulus window
paramsIN.SCC.EndTime_sec=0.40;
paramsIN.SCC.DELAYbinwidth_sec=50e-6;  % 50e-6 is what Joris used
paramsIN.SCC.Duration_sec=paramsIN.SCC.EndTime_sec-paramsIN.SCC.StartTime_sec;

NSCCs=cell(NUMrows,length(levels_dBSPL));
NSCC_delays_usec=cell(NUMrows,length(levels_dBSPL));
NSCC_BFs_kHz=cell(NUMrows,length(levels_dBSPL));
NSCC_avgrates=cell(NUMrows,length(levels_dBSPL));
NSCC_nsps=cell(NUMrows,length(levels_dBSPL));
NSCC_CDs_usec=cell(NUMrows,length(levels_dBSPL));
NSCC_peaks=cell(NUMrows,length(levels_dBSPL));
NSCC_0delay=cell(NUMrows,length(levels_dBSPL));
NSCC_ARBdelay=cell(NUMrows,length(levels_dBSPL));
SCCsMAX=0;
SMP_NSCC_0delay=cell(numSCCs,length(levels_dBSPL));
SMP_NSCC_ARBdelay=cell(numSCCs,length(levels_dBSPL));
SMP_NSCC_CD=cell(numSCCs,length(levels_dBSPL));
SMP_NSCC_peak=cell(numSCCs,length(levels_dBSPL));
for i=1:numSCCs
   for LEVind=1:length(levels_dBSPL)
      SMP_NSCC_0delay{i,LEVind}=NaN*ones(1,length(FeaturesText));
      SMP_NSCC_ARBdelay{i,LEVind}=NaN*ones(1,length(FeaturesText));
      SMP_NSCC_CD{i,LEVind}=NaN*ones(1,length(FeaturesText));
      SMP_NSCC_peak{i,LEVind}=NaN*ones(1,length(FeaturesText));
   end
end


for LEVEL=levels_dBSPL
% beep
% disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
% for LEVEL=levels_dBSPL(end)
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Calcs
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      PERhists{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      PERhistXs_sec{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      
      %%%% Decide which BFs needed for SCCs
      [y,BF_INDEX]=min(abs(unit.Tone_reBF_simFF.BFs_kHz-unit.Info.BF_kHz));  % Finds index of BF from unit.Tone_reBF_simFF.BFs_kHz
      SCC_allBFinds=[]; % Lists all BFinds needed for SCCs
      % SCC is between BF+SCC_octOFFSET and BF-SCC_octOFFSET octaves
      SCCind=1;  % index of SCC to calculate
      [y,BFind1]=min(abs(unit.Tone_reBF_simFF.BFs_kHz-unit.Info.BF_kHz*2^SCC_octOFFSET));
      [y,BFind2]=min(abs(unit.Tone_reBF_simFF.BFs_kHz-unit.Info.BF_kHz*2^-SCC_octOFFSET));
      if abs(log2(unit.Tone_reBF_simFF.BFs_kHz(BFind1)/unit.Tone_reBF_simFF.BFs_kHz(BFind2))-2*SCC_octOFFSET)>0.03
         beep
         disp('*** WARNING ***: The two BFs used for the SCC (BF+SCC_octOFFset & BF-SCC_octOFFset) are more than 0.03 octaves different than the desired distance')
      end
      NSCC_BFinds{SCCind}=[BFind1 BFind2];
      SCC_allBFinds=[SCC_allBFinds NSCC_BFinds{SCCind}];
      % Tally all BFs needed
      SCC_allBFinds=unique(SCC_allBFinds);
      SCC_allSpikeTrains=cell(size(unit.Tone_reBF_simFF.BFs_kHz));  % Store Spike Trains (driven spikes), but on;y those we need for computing SCCs later
      
      for BFind=1:length(unit.Tone_reBF_simFF.BFs_kHz)
         if ~isempty(unit.Tone_reBF_simFF.picNums{LEVind,BFind})
            
            PIC=concatPICS_NOHR(unit.Tone_reBF_simFF.picNums{LEVind,BFind},unit.Tone_reBF_simFF.excludeLines{LEVind,BFind});
            % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
            PIC=simFF_PICshift(PIC);
            PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
            
            BFs_kHz{ROWind,LEVind}(BFind)=unit.Tone_reBF_simFF.BFs_kHz(BFind);
            PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
            PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;
            
            % Filter PERhists with Triangular Filter: filter 3 reps, take middle rep to avoid edge effects and keep periodic
            N=length(PERhists{ROWind,LEVind}{BFind});
            SmoothedPERhist=trifilt(repmat(PERhists{ROWind,LEVind}{BFind},1,3),TFiltWidth);
            PERhists_Smoothed{ROWind,LEVind}{BFind}=SmoothedPERhist(N+1:2*N); % Take middle 3rd
            % Determine Maximum of ALL plotted PERhists (i.e., post-Smoothing)
            if max(PERhists_Smoothed{ROWind,LEVind}{BFind})>PERhistsMAX
               PERhistsMAX=max(PERhists_Smoothed{ROWind,LEVind}{BFind});
            end
            
            %%%%%%%%%%%%% SKIP ALL DFT stuff
            
            Nsps{ROWind,LEVind}(BFind)=PIC.PERhist.NumDrivenSpikes;
            Rates{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);
            if PIC.SynchRate_PERhist.FeatureRaySig
               Synchs{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeatureSynchs;
               Phases{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeaturePhases;
            else
               Synchs{ROWind,LEVind}(BFind)=NaN;
               Phases{ROWind,LEVind}(BFind)=NaN;
            end
            
            %%%% Save SpikeTrains for this BF if it is used for SCCs
            SCCind=find(SCC_allBFinds==BFind);
            if ~isempty(SCCind)
               SCC_allSpikeTrains{BFind}=getDrivenSpikeTrains(PIC.x.spikes{1},[],[paramsIN.SCC.StartTime_sec paramsIN.SCC.EndTime_sec]);
            end
            
         else
            BFs_kHz{ROWind,LEVind}(BFind)=unit.Tone_reBF_simFF.BFs_kHz(BFind);
            Nsps{ROWind,LEVind}(BFind)=NaN;
            Rates{ROWind,LEVind}(BFind)=NaN;
            Synchs{ROWind,LEVind}(BFind)=NaN;
            Phases{ROWind,LEVind}(BFind)=NaN;
            
            %%%% Warn that no SpikeTrains for this BF if it is used for SCCs
            SCCind=find(SCC_allBFinds==BFind);
            if ~isempty(SCCind)
               %                            uiwait(warndlg(sprintf('SCC_allSpikeTrains{SCCind} is set to EMPTY because no data for this BF!\n\n   Feature: TONE; LEVEL=%.f dB SPL, BFind=%d', ...
               %                               levels_dBSPL(LEVind),BFind),'SCC WARNING'))
               disp(sprintf('SCC_allSpikeTrains{SCCind} is set to EMPTY because no data for this BF!\n   Feature: TONE; LEVEL=%.f dB SPL, BFind=%d', ...
                  levels_dBSPL(LEVind),BFind));
            end                        
            
         end
      end % BFinds
      
      %%%%%%%%%%%%%%%%
      % Compute SCCs
      %%%%%%%%%%%%%%%%
      NSCCs{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_delays_usec{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_avgrates{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_nsps{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_BFskHz{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_CDs_usec{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_peaks{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_0delay{ROWind,LEVind}=cell(size(NSCC_BFinds));
      NSCC_ARBdelay{ROWind,LEVind}=cell(size(NSCC_BFinds));
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
            disp(sprintf('Feature: TONE; Level: %.f  --  Computing SCC # %d between BFs %d and %d ........', ...
               levels_dBSPL(LEVind),SCCind,NSCC_BFinds{SCCind}))                     
            [NSCCs{ROWind,LEVind}{SCCind},NSCC_delays_usec{ROWind,LEVind}{SCCind},NSCC_avgrates{ROWind,LEVind}{SCCind},NSCC_nsps{ROWind,LEVind}{SCCind}] ...
               = ShufCrossCorr(SpikeTrains,paramsIN.SCC.DELAYbinwidth_sec,paramsIN.SCC.Duration_sec);
            
            % Determine Maximum of ALL plotted SCCs (i.e., post-Smoothing)
            if max(NSCCs{ROWind,LEVind}{SCCind})>SCCsMAX
               SCCsMAX=max(NSCCs{ROWind,LEVind}{SCCind});
            end
            NSCC_BFskHz{ROWind,LEVind}{SCCind}=[unit.Tone_reBF_simFF.BFs_kHz(NSCC_BFinds{SCCind}(1)),unit.Tone_reBF_simFF.BFs_kHz(NSCC_BFinds{SCCind}(2))];                     
            NSCC_0delay{ROWind,LEVind}{SCCind}=NSCCs{ROWind,LEVind}{SCCind}(find(NSCC_delays_usec{ROWind,LEVind}{SCCind}==0));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Finding CD is TOUGH (like below)!!!! %%%%                                   %
            %                                                                                %
            % Need a robust way to do it                                                     %
            % 1) Find Local Maxima within 5 periods, based on 3-pt smoothed NSCC              %
            % 2) Restrict to only those within 15% of max                                    %
            % 3) Take local max closest to 0                                                 %
            % *** THIS SOMETIMES CHOOSES CD close to 0, when the NSCC is fairly periodic     %
            %     BUT, these are just tough ones to get right without a lot of "knowledge"   %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            F0per_us=1/(unit.Info.BF_kHz*1000)*1e6;
            SCCperiodsTONE=5;
            delaysTEMP_usec=NSCC_delays_usec{ROWind,LEVind}{SCCind}(find((NSCC_delays_usec{ROWind,LEVind}{SCCind}>=-SCCperiodsTONE*F0per_us)&(NSCC_delays_usec{ROWind,LEVind}{SCCind}<=SCCperiodsTONE*F0per_us)));
            nsccTEMP=NSCCs{ROWind,LEVind}{SCCind}(find((NSCC_delays_usec{ROWind,LEVind}{SCCind}>=-SCCperiodsTONE*F0per_us)&(NSCC_delays_usec{ROWind,LEVind}{SCCind}<=SCCperiodsTONE*F0per_us)));
            
            %%%%%%%%%%% LOCAL MAXIMA
            NSCC_TFiltWidth=3;
            %%%% find local max within 1 period, based on smoothed SCCs
            diff_LtoR=diff(trifilt(nsccTEMP,NSCC_TFiltWidth));
            diff_RtoL=diff(fliplr(trifilt(nsccTEMP,NSCC_TFiltWidth)));
            diffDelays_LtoR=delaysTEMP_usec(1:end-1);
            diffDelays_RtoL=fliplr(delaysTEMP_usec(2:end));
            LocalMaxDelays1=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
            LocalMaxDelays=intersect(LocalMaxDelays1,delaysTEMP_usec(find(nsccTEMP>=1)));
            nsccLocalMax=zeros(size(LocalMaxDelays));
            for i=1:length(LocalMaxDelays)
               nsccLocalMax(i)=nsccTEMP(find(delaysTEMP_usec==LocalMaxDelays(i)));
            end
            
            %% Restrict to local max within 15% of peak
            PercCRIT=0.15;
            RESTRICTinds=find(nsccLocalMax>max(nsccLocalMax)*(1-PercCRIT));
            LocalMaxDelaysRESTRICT=LocalMaxDelays(RESTRICTinds);
            nsccLocalMaxRESTRICT=nsccLocalMax(RESTRICTinds);
            
            [y,i]=min(abs(LocalMaxDelaysRESTRICT));
            NSCC_CDs_usec{ROWind,LEVind}{SCCind}=LocalMaxDelaysRESTRICT(i);
            NSCC_peaks{ROWind,LEVind}{SCCind}=nsccLocalMaxRESTRICT(i);  % This is not actually the peak, but the NSCC at the CD
         else  % No Data to compute SCCs
            F0per_us=1/(unit.Info.BF_kHz*1000)*1e6;
            disp(sprintf('Feature: TONE; Level: %.f  --  CAN''T compute SCC # %d between BFs %d and %d (MISSING DATA) ........', ...
               levels_dBSPL(LEVind),SCCind,NSCC_BFinds{SCCind}))                     
            NSCCs{ROWind,LEVind}{SCCind}=NaN*ones(1,3);
            NSCC_delays_usec{ROWind,LEVind}{SCCind}=[-F0per_us 0 F0per_us];
            NSCC_avgrates{ROWind,LEVind}{SCCind}=NaN;
            NSCC_nsps{ROWind,LEVind}{SCCind}=0;
            NSCC_BFskHz{ROWind,LEVind}{SCCind}=[unit.Tone_reBF_simFF.BFs_kHz(NSCC_BFinds{SCCind}(1)),unit.Tone_reBF_simFF.BFs_kHz(NSCC_BFinds{SCCind}(2))];                     
            NSCC_0delay{ROWind,LEVind}{SCCind}=NaN;
            NSCC_ARBdelay{ROWind,LEVind}{SCCind}=NaN;
            NSCC_CDs_usec{ROWind,LEVind}{SCCind}=NaN;
            NSCC_peaks{ROWind,LEVind}{SCCind}=NaN;
         end
      end
     
      % Determine how many CHANNELS to plot
      if length(unit.Tone_reBF_simFF.BFs_kHz)>PERhistsYCHANS
         PERhistsYCHANS=length(unit.Tone_reBF_simFF.BFs_kHz);
      end
   end   % if Tone data
   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Calcs
   if isfield(unit,'EH_reBF_simFF')
      for FeatIND=FeatINDs
         ROWind=ROWind+1;
         for HarmonicsIND=1:1
            for PolarityIND=1:1
               eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
               if ~isempty(yTEMP)
                  %%%% EH_reBF plots
                  LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                  PERhists{ROWind,LEVind}=cell(size(yTEMP.BFs_kHz));
                  PERhistXs_sec{ROWind,LEVind}=cell(size(yTEMP.BFs_kHz));
                  
                  %%%% Decide which BFs needed for SCCs
                  [y,BF_INDEX]=min(abs(yTEMP.BFs_kHz-unit.Info.BF_kHz));  % Finds index of BF from unit.Tone_reBF_simFF.BFs_kHz
                  SCC_allBFinds=[]; % Lists all BFinds needed for SCCs
                  % SCC is between BF+SCC_octOFFSET and BF-SCC_octOFFSET octaves
                  SCCind=1;  % index of SCC to calculate
                  [y,BFind1]=min(abs(yTEMP.BFs_kHz-unit.Info.BF_kHz*2^SCC_octOFFSET));
                  [y,BFind2]=min(abs(yTEMP.BFs_kHz-unit.Info.BF_kHz*2^-SCC_octOFFSET));
                  if abs(log2(yTEMP.BFs_kHz(BFind1)/yTEMP.BFs_kHz(BFind2))-2*SCC_octOFFSET)>0.03
                     beep
                     disp('*** WARNING ***: The two BFs used for the SCC (BF+SCC_octOFFset & BF-SCC_octOFFset) are more than 0.03 octaves different than the desired distance')
                  end
                  NSCC_BFinds{SCCind}=[BFind1 BFind2];
                  SCC_allBFinds=[SCC_allBFinds NSCC_BFinds{SCCind}];
                  % Tally all BFs needed
                  SCC_allBFinds=unique(SCC_allBFinds);
                  SCC_allSpikeTrains=cell(size(yTEMP.BFs_kHz));  % Store Spike Trains (driven spikes), but on;y those we need for computing SCCs later
                  
                  for BFind=1:length(yTEMP.BFs_kHz)
                     if ~isempty(yTEMP.picNums{LEVind,BFind})
                        PIC=concatPICS_NOHR(yTEMP.picNums{LEVind,BFind},yTEMP.excludeLines{LEVind,BFind});
                        % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                        PIC=simFF_PICshift(PIC);
                        PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
                        
                        BFs_kHz{ROWind,LEVind}(BFind)=yTEMP.BFs_kHz(BFind);
                        PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
                        PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;

                        % Filter PERhists with Triangular Filter: filter 3 reps, take middle rep to avoid edge effects and keep periodic
                        N=length(PERhists{ROWind,LEVind}{BFind});
                        SmoothedPERhist=trifilt(repmat(PERhists{ROWind,LEVind}{BFind},1,3),TFiltWidth);
                        PERhists_Smoothed{ROWind,LEVind}{BFind}=SmoothedPERhist(N+1:2*N);
                        % Determine Maximum of ALL plotted PERhists (i.e., post-Smoothing)
                        if max(PERhists_Smoothed{ROWind,LEVind}{BFind})>PERhistsMAX
                           PERhistsMAX=max(PERhists_Smoothed{ROWind,LEVind}{BFind});
                        end                        
                        
                        % Save DFTs as well of PERhists
                        DFTs{ROWind,LEVind}{BFind}=PIC.SynchRate_PERhist.SynchRate_PERhist;
                        DFTfreqs_Hz{ROWind,LEVind}{BFind}=PIC.SynchRate_PERhist.FFTfreqs;
                        % Determine Maximum of ALL plotted DFTs
                        if max(abs(DFTs{ROWind,LEVind}{BFind}))>DFTsMAX
                           DFTsMAX=max(abs(DFTs{ROWind,LEVind}{BFind}));
                        end                        

                        Nsps{ROWind,LEVind}(BFind)=PIC.PERhist.NumDrivenSpikes;
                        Rates{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);
                        if PIC.SynchRate_PERhist.FeatureRaySig(FeatIND)
                           Synchs{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeatureSynchs(FeatIND);
                           Phases{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeaturePhases(FeatIND);
                        else
                           Synchs{ROWind,LEVind}(BFind)=NaN;
                           Phases{ROWind,LEVind}(BFind)=NaN;
                        end
                        
                        %%%% Save SpikeTrains for this BF if it is used for SCCs
                        SCCind=find(SCC_allBFinds==BFind);
                        if ~isempty(SCCind)
                           SCC_allSpikeTrains{BFind}=getDrivenSpikeTrains(PIC.x.spikes{1},[],[paramsIN.SCC.StartTime_sec paramsIN.SCC.EndTime_sec]);
                        end
                        
                     else
                        BFs_kHz{ROWind,LEVind}(BFind)=yTEMP.BFs_kHz(BFind);
                        Nsps{ROWind,LEVind}(BFind)=NaN;
                        Rates{ROWind,LEVind}(BFind)=NaN;
                        Synchs{ROWind,LEVind}(BFind)=NaN;
                        Phases{ROWind,LEVind}(BFind)=NaN;
                        
                        PERhists{ROWind,LEVind}{BFind}=[];
                        PERhistXs_sec{ROWind,LEVind}{BFind}=[];
                        DFTs{ROWind,LEVind}{BFind}=[];
                        DFTfreqs_Hz{ROWind,LEVind}{BFind}=[];

                        %%%% Warn that no SpikeTrains for this BF if it is used for SCCs
                        SCCind=find(SCC_allBFinds==BFind);
                        if ~isempty(SCCind)
                           %                            uiwait(warndlg(sprintf('SCC_allSpikeTrains{SCCind} is set to EMPTY because no data for this BF!\n\n   Feature: %s; LEVEL=%.f dB SPL, BFind=%d', ...
                           %                               FeaturesText{FeatIND},levels_dBSPL(LEVind),BFind),'SCC WARNING'))
                           disp(sprintf('SCC_allSpikeTrains{SCCind} is set to EMPTY because no data for this BF!\n   Feature: %s; LEVEL=%.f dB SPL, BFind=%d', ...
                              FeaturesText{FeatIND},levels_dBSPL(LEVind),BFind));
                        end                        
                     end
                  end % BFinds

                  %%%%%%%%%%%%%%%%
                  % Calculate ALSR data
                  %%%%%%%%%%%%%%%%
                  if ~isempty(DFTfreqs_Hz{ROWind,LEVind})
                     ALSR_OCTrange=0.28;  % Slight slop for slight mismatches in sampling rates
                     ALSRinds=find((yTEMP.BFs_kHz>=unit.Info.BF_kHz*2^-ALSR_OCTrange)&(yTEMP.BFs_kHz<=unit.Info.BF_kHz*2^ALSR_OCTrange));
                     if length(ALSRinds)<length(yTEMP.BFs_kHz)
                        %                         warndlg('Not all channels taken in ALSR')  % TAKE OUT eventually, just want to know when we're not getting whole story
                        disp('*** WARNING ***: Not all channels taken in ALSR')  % TAKE OUT eventually, just want to know when we're not getting whole story
                     end
                     SynchRatesTEMP=NaN*ones(1,length(yTEMP.BFs_kHz));
                     for BFind=ALSRinds
                        [y,DFT_INDEX]=min(abs(DFTfreqs_Hz{ROWind,LEVind}{BFind}-unit.Info.BF_kHz*1000));
                        if ~isempty(DFT_INDEX)
                           SynchRatesTEMP(BFind)=abs(DFTs{ROWind,LEVind}{BFind}(DFT_INDEX));
                        end
                     end
                     ALSRs{ROWind,LEVind}=mean(SynchRatesTEMP(~isnan(SynchRatesTEMP)));
                  else
                     ALSRs{ROWind,LEVind}=NaN;
                  end
                  
                  %%%%%%%%%%%%%%%%
                  % Compute SCCs
                  %%%%%%%%%%%%%%%%
                  NSCCs{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_delays_usec{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_avgrates{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_nsps{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_BFskHz{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_CDs_usec{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_peaks{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_0delay{ROWind,LEVind}=cell(size(NSCC_BFinds));
                  NSCC_ARBdelay{ROWind,LEVind}=cell(size(NSCC_BFinds));
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
                        disp(sprintf('Feature: %s; Level: %.f  --  Computing SCC # %d between BFs %d and %d ........', ...
                           FeaturesText{FeatIND},levels_dBSPL(LEVind),SCCind,NSCC_BFinds{SCCind}))                     
                        [NSCCs{ROWind,LEVind}{SCCind},NSCC_delays_usec{ROWind,LEVind}{SCCind},NSCC_avgrates{ROWind,LEVind}{SCCind},NSCC_nsps{ROWind,LEVind}{SCCind}] ...
                           = ShufCrossCorr(SpikeTrains,paramsIN.SCC.DELAYbinwidth_sec,paramsIN.SCC.Duration_sec);
                        
                        % Determine Maximum of ALL plotted SCCs (i.e., post-Smoothing)
                        if max(NSCCs{ROWind,LEVind}{SCCind})>SCCsMAX
                           SCCsMAX=max(NSCCs{ROWind,LEVind}{SCCind});
                        end
                        NSCC_BFskHz{ROWind,LEVind}{SCCind}=[yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(1)),yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(2))];                     
                        NSCC_0delay{ROWind,LEVind}{SCCind}=NSCCs{ROWind,LEVind}{SCCind}(find(NSCC_delays_usec{ROWind,LEVind}{SCCind}==0));
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%% Finding CD is TOUGH!!!! %%%%                                                %
                        %                                                                                %
                        % Need a robust way to do it                                                     %
                        % 1) Find Local Maxima within 1 period, based on 3-pt smoothed NSCC              %
                        % 2) Restrict to only those within 15% of max                                    %
                        % 3) Take local max closest to 0                                                 %
                        % *** THIS SOMETIMES CHOOSES CD close to 0, when the NSCC is fairly periodic     %
                        %     BUT, these are just tough ones to get right without a lot of "knowledge"   %
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        F0per_us=1/yTEMP.FeatureFreqs_Hz{1}(1)*1e6;
                        delaysTEMP_usec=NSCC_delays_usec{ROWind,LEVind}{SCCind}(find((NSCC_delays_usec{ROWind,LEVind}{SCCind}>=-F0per_us)&(NSCC_delays_usec{ROWind,LEVind}{SCCind}<=F0per_us)));
                        nsccTEMP=NSCCs{ROWind,LEVind}{SCCind}(find((NSCC_delays_usec{ROWind,LEVind}{SCCind}>=-F0per_us)&(NSCC_delays_usec{ROWind,LEVind}{SCCind}<=F0per_us)));
                        
                        %%%%%%%%%%% LOCAL MAXIMA
                        NSCC_TFiltWidth=3;
                        %%%% find local max within 1 period, based on smoothed SCCs
                        diff_LtoR=diff(trifilt(nsccTEMP,NSCC_TFiltWidth));
                        diff_RtoL=diff(fliplr(trifilt(nsccTEMP,NSCC_TFiltWidth)));
                        diffDelays_LtoR=delaysTEMP_usec(1:end-1);
                        diffDelays_RtoL=fliplr(delaysTEMP_usec(2:end));
                        LocalMaxDelays1=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
                        LocalMaxDelays=intersect(LocalMaxDelays1,delaysTEMP_usec(find(nsccTEMP>=1)));
                        nsccLocalMax=zeros(size(LocalMaxDelays));
                        for i=1:length(LocalMaxDelays)
                           nsccLocalMax(i)=nsccTEMP(find(delaysTEMP_usec==LocalMaxDelays(i)));
                        end
                        
                        %% Restrict to local max within 15% of peak
                        PercCRIT=0.15;
                        RESTRICTinds=find(nsccLocalMax>max(nsccLocalMax)*(1-PercCRIT));
                        LocalMaxDelaysRESTRICT=LocalMaxDelays(RESTRICTinds);
                        nsccLocalMaxRESTRICT=nsccLocalMax(RESTRICTinds);
                        
                        [y,i]=min(abs(LocalMaxDelaysRESTRICT));
                        NSCC_CDs_usec{ROWind,LEVind}{SCCind}=LocalMaxDelaysRESTRICT(i);
                        NSCC_peaks{ROWind,LEVind}{SCCind}=nsccLocalMaxRESTRICT(i);  % This is not actually the peak, but the NSCC at the CD
                     else  % No Data to compute SCCs
                        F0per_us=1/yTEMP.FeatureFreqs_Hz{1}(1)*1e6;
                        disp(sprintf('Feature: %s; Level: %.f  --  CAN''T compute SCC # %d between BFs %d and %d (MISSING DATA) ........', ...
                           FeaturesText{FeatIND},levels_dBSPL(LEVind),SCCind,NSCC_BFinds{SCCind}))                     
                        NSCCs{ROWind,LEVind}{SCCind}=NaN*ones(1,3);
                        NSCC_delays_usec{ROWind,LEVind}{SCCind}=[-F0per_us 0 F0per_us];
                        NSCC_avgrates{ROWind,LEVind}{SCCind}=NaN;
                        NSCC_nsps{ROWind,LEVind}{SCCind}=0;
                        NSCC_BFskHz{ROWind,LEVind}{SCCind}=[yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(1)),yTEMP.BFs_kHz(NSCC_BFinds{SCCind}(2))];                     
                        NSCC_0delay{ROWind,LEVind}{SCCind}=NaN;
                        NSCC_ARBdelay{ROWind,LEVind}{SCCind}=NaN;
                        NSCC_CDs_usec{ROWind,LEVind}{SCCind}=NaN;
                        NSCC_peaks{ROWind,LEVind}{SCCind}=NaN;
                     end
                     
                     %                   % TO COMPUTE
                     %                   *NSCCs
                     %                   *NSCC_delays_usec
                     %                   *NSCC_BFs_kHz
                     %                   *NSCC_avgrates
                     %                   *NSCC_nsps
                     %                   *NSCC_CDs_usec
                     %                   *NSCC_peaks
                     %                   *NSCC_0delay
                     %                   *(LATER in SMP plot) NSCC_ARBdelay
                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
                     %%%%%%%%%%%%%%%% 2/5/05 TODO %%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % *1) Compute peak, CD (use 1/F0 to limit search?)
                     % *2) Setup *NSCC and SCC plots,
                     %         *and plot SCC and NSCC values
                     % *3) Set up SMP plots
                     % 3.5) Verify all calcs look OK, decide how to handle NSCC/SCC
                     % 4) Start to looking at data to see effects!!!!
                     %      - look for CNL
                     %      - look for robust SCCs
                     %              - !!! LOOKS LIKE D&G story holds, largest Cross-Corr is at troughs, drops at formants
                     %      - DO WE NEED 2 SCCs, or just one BF+epsilon, BF-epsilon???  WOULD BE SIMPLER!!!
                     %      - DEVELOP STORY
                     
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     %%%%%%%%%%%%%%%%%%%% TODO  2/2/05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     % Compute these 2 SCCs for demo, and get all things working
                     % TRY to setup general, simple implementation of these things 
                     % Then can start looking for real issues
                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                     
                  end
                  
                  
                  %%%%%%%%%%%%%%%%
                  % Store SMP data
                  %%%%%%%%%%%%%%%%
                  SMP_rate{LEVind}(FeatIND)=Rates{ROWind,LEVind}(BF_INDEX);
                  SMP_alsr{LEVind}(FeatIND)=ALSRs{ROWind,LEVind};
                  
                  for i=1:numSCCs
                     SMP_NSCC_0delay{i,LEVind}(FeatIND)=NSCC_0delay{ROWind,LEVind}{i};
                     SMP_NSCC_CD{i,LEVind}(FeatIND)=NSCC_CDs_usec{ROWind,LEVind}{i};
                     SMP_NSCC_peak{i,LEVind}(FeatIND)=NSCC_peaks{ROWind,LEVind}{i};
                  end
                  
                  % Determine how many CHANNELS to plot
                  if length(yTEMP.BFs_kHz)>PERhistsYCHANS
                     PERhistsYCHANS=length(yTEMP.BFs_kHz);
                  end
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end % If EHrBF data
end % Levels


close all
UnitVerify_TC(ExpDate,UnitName,1);  % Re-show TC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL PERhist PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot data
LEVELcolors={'b','r','g','c','k','y'};
FEATmarkers={'.','x','s','^','*','<','^','>'};
FIG.FontSize=8;

FeatureColors={'r','g'};
%%% Show multiple periods of tone PERhist
NperiodsTONE=5;
if isfield(unit,'EH_reBF_simFF')
   PERhist_XMAX=1/F0min*1000;  % take the PERhist_Xlims on vowels (in ms)
   XLIMITS_perhist=[0 PERhist_XMAX];
else
   % Need to find PERhist_XMAX by hand
   PERhist_XMAX=0;
   for ROWind=1:NUMrows
      for LEVind=1:length(levels_dBSPL)
         for BFind=1:length(PERhistXs_sec{1,1})
            if max(PERhistXs_sec{ROWind,LEVind}{BFind})>PERhist_XMAX
               PERhist_XMAX=max(PERhistXs_sec{ROWind,LEVind}{BFind});
            end   
         end
      end
   end
   PERhist_XMAX=PERhist_XMAX*1000; % in ms
   XLIMITS_perhist=[0 PERhist_XMAX*NperiodsTONE];
end
XLIMITS_rate=[0 300];
XLIMITS_synch=[0 1];
XLIMITS_phase=[-pi pi];
%% Find PERhist_YMAX
PERhistGAIN=3.0; % # of channels covered by max PERhist
PERhists_logCHwidth=log10(highBF/lowBF)/(PERhistsYCHANS-1);  % log10 channel width
PERhist_YMIN=lowBF;
PERhist_YMAX=10^((PERhistsYCHANS-1+PERhistGAIN)*PERhists_logCHwidth)*lowBF;    % sets an extra (GAIN-1) log channel widths
YLIMITS=[PERhist_YMIN PERhist_YMAX];  % Used for all plots
%% This  is ALL needed to get the right LOG Yticks!!
YLIMunit=10^floor(log10(lowBF));
YLIMS=floor(lowBF/YLIMunit)*YLIMunit*[1 100]; % Do two decades to be sure we get all the ticks
YTICKS=[YLIMS(1):YLIMunit:YLIMS(1)*10 YLIMS(1)*20:YLIMunit*10:YLIMS(end)];
BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave

% Setup parameters for title
if isempty(unit.Info.Threshold_dBSPL)
   unit.Info.Threshold_dBSPL=NaN;
end
if isempty(unit.Info.SR_sps)
   unit.Info.SR_sps=NaN;
end
if isempty(unit.Info.Q10)
   unit.Info.Q10=NaN;
end

for LEVEL=levels_dBSPL
% beep
% disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
% for LEVEL=levels_dBSPL(end)
   figure(round(LEVEL)); clf
   set(gcf,'pos',[420     4   977   976])
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Plots
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      
      PLOTnum=(ROWind-1)*NUMcols+1;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      for BFind=1:length(BFs_kHz{ROWind,LEVind})   
         if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
            LINEwidth=2;
         else
            LINEwidth=.5;
         end
         % This normalization plots each signal the same size on a log scale
         if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
            NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;
            %%% Show multiple periods of tone
            %             NperiodsTONE=5;
            DT=diff(PERhistXs_sec{ROWind,LEVind}{BFind}(end-1:end));
            clear xdata ydata
            xdata=PERhistXs_sec{ROWind,LEVind}{BFind};
            for i=2:NperiodsTONE
               xdata=[xdata PERhistXs_sec{ROWind,LEVind}{BFind}+xdata(end)+DT/2];
            end
            xdata=xdata*1000;
            ydata=repmat(PERhists_Smoothed{ROWind,LEVind}{BFind}*NormFact+BFs_kHz{ROWind,LEVind}(BFind),1,NperiodsTONE);
            semilogy(xdata,ydata,'LineWidth',LINEwidth)
            hold on
         end
      end
      semilogy(XLIMITS_perhist,unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Time (ms)')
      ylabel('Effective Best Frequency (kHz)')
      title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL', ...
         ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,'TONE', ...
         unit.Tone_reBF_simFF.levels_dBSPL(LEVind)),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
      %       xlim([0 max(PERhistXs_sec{ROWind,LEVind}{BFind}*1000)]) % Different xlim for TONES ??? 
      xlim(XLIMITS_perhist)
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
         text(PERhistXs_sec{ROWind,LEVind}{BFind}(end)*1000,YLIMITS(1),'1/f','units','data','HorizontalAlignment','center','VerticalAlignment','top')
      end
      hold off
      
      %%%% Rate Plot
      PLOTnum=(ROWind-1)*NUMcols+2;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]'))
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_rate)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Synch Plot
      PLOTnum=(ROWind-1)*NUMcols+3;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Synch Coef (to f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_synch)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Phase Plot
      PLOTnum=(ROWind-1)*NUMcols+4;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Phase (cycles of f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_phase)
      set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
   end   
   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Plots
   if isfield(unit,'EH_reBF_simFF')
      for FeatIND=FeatINDs
         ROWind=ROWind+1;
         for HarmonicsIND=1:1
            for PolarityIND=1:1
               eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
               if ~isempty(yTEMP)
                  %%%% EH_reBF plots
                  LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                  
                  %%%% Spatio-Temporal Plots
                  PLOTnum=(ROWind-1)*NUMcols+1;
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  for BFind=1:length(BFs_kHz{ROWind,LEVind})
                     if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                        LINEwidth=2;
                     else
                        LINEwidth=.5;
                     end
                     % This normalization plots each signal the same size on a log scale
                     if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;
                        semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000, ...
                           PERhists_Smoothed{ROWind,LEVind}{BFind}*NormFact+BFs_kHz{ROWind,LEVind}(BFind), ...
                           'LineWidth',LINEwidth)
                        hold on
                     end
                  end
                  xlabel('Time (ms)')
                  ylabel('Effective Best Frequency (kHz)')                  
                  if ROWind==1
                     title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                        ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  else
                     title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  end
                  xlim(XLIMITS_perhist)
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                           sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                           'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                     for BFind=1:length(BFs_kHz{ROWind,LEVind})
                        if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                           if (FeatINDPlot<=FeatIND)
                              if (FeatINDPlot>1)
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                              else
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                              end
                           end
                        end
                     end
                  end
                  hold off
                  
                  
                  %%%% Rate Plot
                  PLOTnum=(ROWind-1)*NUMcols+2;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
                  semilogy(ALSRs{ROWind,LEVind},unit.Info.BF_kHz,'go','MarkerSize',6)
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'),'FontSize',6)
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_rate)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Synch Plot
                  PLOTnum=(ROWind-1)*NUMcols+3;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_synch)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Phase Plot
                  PLOTnum=(ROWind-1)*NUMcols+4;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_phase)
                  set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_phase,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end
   
   
   Xcorner=0.05;
   Xwidth1=.5;
   Xshift1=0.05;
   Xwidth2=.1;
   Xshift2=0.03;
   
   Ycorner=0.05;
   Yshift=0.07;
   Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
   
   TICKlength=0.02;
   
   if NUMrows>4
      set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>3
      set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>2
      set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>1
      set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   
   orient landscape
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2/13/05 TODO
% *1) plot 5 periods tone
% *2) Add all-levels PERhist plot for tones
% *3) Take out DFT for tones
% *4) Update SCCs - which channels
% *5) Add SCCs for tones
% *5.25) Add extra figure with all level SCCs
% 5.5) Figure out how to show CNL is a benefit, by accounting for passive TW delay
% 6) Deal with passive TW delay, to show NSCCpeak goes up with level
% 7) Add figure to show CD(L) for tones


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADD EXTRA PERhist plot at the end with all levels on top of one anothe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1000); clf
set(gcf,'pos',[420     4   977   976])
ROWind=0;
ALLlevelsTriFiltTONE=9;
ALLlevelsTriFilt=3;

LEVEL=max(levels_dBSPL);
%%%%%%%%%%%%%%%%%%%% Tone Plots
if isfield(unit,'Tone_reBF_simFF')
   ROWind=ROWind+1;
   
   %%%% Tone_reBF plots
   LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
   
   PLOTnum=(ROWind-1)*NUMcols+1;   
   eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
   LEGtext='';
   for BFind=1:length(BFs_kHz{ROWind,LEVind})   
      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
         LINEwidth=.5;
      else
         LINEwidth=.125;
      end
      % This normalization plots each signal the same size on a log scale
      if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
         NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;
         %%% Show multiple periods of tone
         DT=diff(PERhistXs_sec{ROWind,LEVind}{BFind}(end-1:end));
         clear xdata ydata
         xdata=PERhistXs_sec{ROWind,LEVind}{BFind};
         for i=2:NperiodsTONE
            xdata=[xdata PERhistXs_sec{ROWind,LEVind}{BFind}+xdata(end)+DT/2];
         end
         xdata=xdata*1000;
         ydata=repmat(PERhists_Smoothed{ROWind,LEVind}{BFind}*NormFact+BFs_kHz{ROWind,LEVind}(BFind),1,NperiodsTONE);
         semilogy(xdata,trifilt(ydata,ALLlevelsTriFiltTONE),'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
         hold on
      end
      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
         LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
      end
      for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
         if ~isempty(PERhistXs_sec{ROWind,LEVind2}{BFind})                             
            %%% Show multiple periods of tone
            DT=diff(PERhistXs_sec{ROWind,LEVind2}{BFind}(end-1:end));
            clear xdata ydata
            xdata=PERhistXs_sec{ROWind,LEVind2}{BFind};
            for i=2:NperiodsTONE
               xdata=[xdata PERhistXs_sec{ROWind,LEVind2}{BFind}+xdata(end)+DT/2];
            end
            xdata=xdata*1000;
            ydata=repmat(PERhists_Smoothed{ROWind,LEVind2}{BFind}*NormFact+BFs_kHz{ROWind,LEVind2}(BFind),1,NperiodsTONE);
            semilogy(xdata,trifilt(ydata,ALLlevelsTriFiltTONE),'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
            hold on
            if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
               LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind2));
            end
         end
      end
      %          if strcmp(FeaturesText{FeatIND},'F1')
      hleg1000a=legend(LEGtext,1);
      set(hleg1000a,'FontSize',8)
      set(hleg1000a,'pos',[0.4451    0.8913    0.0942    0.0473])
      %          end
   end
   semilogy(XLIMITS_perhist,unit.Info.BF_kHz*[1 1],'k:')
   xlabel('Time (ms)')
   ylabel('Effective Best Frequency (kHz)')
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s', ...
      ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,'TONE'),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
   %       xlim([0 max(PERhistXs_sec{ROWind,LEVind}{BFind}*1000)]) % Different xlim for TONES ??? 
   xlim(XLIMITS_perhist)
   PLOThand=eval(['h' num2str(PLOTnum)]);
   set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
   ylim(YLIMITS)  % Same Ylimits for all plots
   text(XLIMITS_perhist(2),YLIMITS(1),'1/f','units','data','HorizontalAlignment','center','VerticalAlignment','top')
   hold off
   
   %%%% Rate Plot
   PLOTnum=(ROWind-1)*NUMcols+2;   
   eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
   semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
   hold on
   for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
      semilogy(Rates{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
      %                   semilogy(Nsps{ROWind,LEVind2}/10,BFs_kHz{ROWind,LEVind2},'m+','MarkerSize',4,'Color',LEVELcolors{LEVind2})
   end
   %    semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
   semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
   xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]'))
   PLOThand=eval(['h' num2str(PLOTnum)]);
   xlim(XLIMITS_rate)
   set(PLOThand,'XDir','reverse')
   set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
   ylim(YLIMITS)  % Same Ylimits for all plots
   hold off
   
   %%%% Synch Plot
   PLOTnum=(ROWind-1)*NUMcols+3;   
   eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
   semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
   hold on
   for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
      semilogy(Synchs{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
   end
   semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
   xlabel('Synch Coef (to f)')
   PLOThand=eval(['h' num2str(PLOTnum)]);
   xlim(XLIMITS_synch)
   set(PLOThand,'XDir','reverse')
   set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
   set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
   ylim(YLIMITS)  % Same Ylimits for all plots
   hold off
   
   %%%% Phase Plot
   PLOTnum=(ROWind-1)*NUMcols+4;   
   eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
   semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
   hold on
   for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
      semilogy(Phases{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
   end
   semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
   xlabel('Phase (cycles of f)')
   PLOThand=eval(['h' num2str(PLOTnum)]);
   xlim(XLIMITS_phase)
   set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
   set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
   ylim(YLIMITS)  % Same Ylimits for all plots
   hold off
end   


%%%%%%%%%%%%%%%%%%%% EH_reBF Plots
if isfield(unit,'EH_reBF_simFF')
   for FeatIND=FeatINDs
      ROWind=ROWind+1;
      for HarmonicsIND=1:1
         for PolarityIND=1:1
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               %%%% EH_reBF plots
               LEVind=find(yTEMP.levels_dBSPL==LEVEL);
               
               %%%% Spatio-Temporal Plots
               PLOTnum=(ROWind-1)*NUMcols+1;
               eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
               LEGtext='';
               for BFind=1:length(BFs_kHz{ROWind,LEVind})
                  if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                     LINEwidth=2;
                  else
                     LINEwidth=.5;
                  end
                  % This normalization plots each signal the same size on a log scale
                  if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                     NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;
                     semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000, ...
                        trifilt(PERhists_Smoothed{ROWind,LEVind}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,LEVind}(BFind), ...
                        'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
                     hold on
                     if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                        LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
                     end
                     for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                        if ~isempty(PERhistXs_sec{ROWind,LEVind2}{BFind})                             
                           semilogy(PERhistXs_sec{ROWind,LEVind2}{BFind}*1000, ...
                              trifilt(PERhists_Smoothed{ROWind,LEVind2}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,LEVind2}(BFind), ...
                              'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
                           if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                              LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind2));
                           end
                        end
                     end
                     if strcmp(FeaturesText{FeatIND},'F1')
                        hleg1000=legend(LEGtext,1);
                        set(hleg1000,'FontSize',8)
                        set(hleg1000,'pos',[0.4451    0.8913    0.0942    0.0473])
                     end
                  end
               end
               xlabel('Time (ms)')
               ylabel('Effective Best Frequency (kHz)')                  
               if ROWind==1
                  title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s,   Harm: %d, Polarity: %d', ...
                     ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                     HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
               else
                  title(sprintf('%s,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                     HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
               end
               xlim(XLIMITS_perhist)
               PLOThand=eval(['h' num2str(PLOTnum)]);
               set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
               ylim(YLIMITS)  % Same Ylimits for all plots
               %%%%%%%%%%%%%%%%%%%%%
               % Plot lines at all features
               for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                  if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                     semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                        sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                        'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  end
                  for BFind=1:length(BFs_kHz{ROWind,LEVind})
                     if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        if (FeatINDPlot<=FeatIND)
                           if (FeatINDPlot>1)
                              text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                           else
                              text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                           end
                        end
                     end
                  end
               end
               hold off
               
               
               %%%% Rate Plot
               PLOTnum=(ROWind-1)*NUMcols+2;   
               eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
               semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
               hold on
               %                semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4,'Color',LEVELcolors{LEVind})
               semilogy(ALSRs{ROWind,LEVind},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',LEVELcolors{LEVind})
               for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                  semilogy(Rates{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
                  %                   semilogy(Nsps{ROWind,LEVind2}/10,BFs_kHz{ROWind,LEVind2},'m+','MarkerSize',4,'Color',LEVELcolors{LEVind2})
                  semilogy(ALSRs{ROWind,LEVind2},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',LEVELcolors{LEVind2})
               end
               semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
               %                xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'),'FontSize',6)               
               xlabel(sprintf('Rate (sp/sec)\nO: ALSR'),'FontSize',8)
               PLOThand=eval(['h' num2str(PLOTnum)]);
               xlim(XLIMITS_rate)
               set(PLOThand,'XDir','reverse')
               set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
               ylim(YLIMITS)  % Same Ylimits for all plots
               %%%%%%%%%%%%%%%%%%%%%
               % Plot lines at all features
               for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                  if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                     semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  end
               end
               hold off
               
               %%%% Synch Plot
               PLOTnum=(ROWind-1)*NUMcols+3;   
               eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
               semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
               hold on
               for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                  semilogy(Synchs{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
               end
               semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
               xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
               PLOThand=eval(['h' num2str(PLOTnum)]);
               xlim(XLIMITS_synch)
               set(PLOThand,'XDir','reverse')
               set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
               set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
               ylim(YLIMITS)  % Same Ylimits for all plots
               %%%%%%%%%%%%%%%%%%%%%
               % Plot lines at all features
               for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                  if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                     semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  end
               end
               hold off
               
               %%%% Phase Plot
               PLOTnum=(ROWind-1)*NUMcols+4;   
               eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
               semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-','Color',LEVELcolors{LEVind})
               hold on
               for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                  semilogy(Phases{ROWind,LEVind2},BFs_kHz{ROWind,LEVind2},'*-','Color',LEVELcolors{LEVind2})
               end
               semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
               xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
               PLOThand=eval(['h' num2str(PLOTnum)]);
               xlim(XLIMITS_phase)
               set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
               set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
               ylim(YLIMITS)  % Same Ylimits for all plots
               %%%%%%%%%%%%%%%%%%%%%
               % Plot lines at all features
               for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                  if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                     semilogy(XLIMITS_phase,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  end
               end
               hold off
               
            end %End if data for this condition, plot
         end % End Polarity
      end % End Harmonics
   end % End Feature
end


Xcorner=0.05;
Xwidth1=.5;
Xshift1=0.05;
Xwidth2=.1;
Xshift2=0.03;

Ycorner=0.05;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

TICKlength=0.02;

if NUMrows>4
   set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>3
   set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>2
   set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>1
   set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

orient landscape

if isfield(unit,'EH_reBF_simFF')

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%% DO ALL DFT PLOTS
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   XLIMITS_dft=[0 10];
   
   for LEVEL=levels_dBSPL
      % beep
      % disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
      % for LEVEL=levels_dBSPL(end)
      figure(round(LEVEL)+1); clf
      set(gcf,'pos',[420     4   977   976])
      ROWind=0;
      
      %%%%%%%%%%%%%%%%%%%% Tone Plots
      %
      %%%% SKIP TONE DFTs!!!!!!!!!!!!!!!!  
      %
      if isfield(unit,'Tone_reBF_simFF')
         ROWind=ROWind+1;
         
         %%%% Tone_reBF plots
         LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
         
         PLOTnum=(ROWind-1)*NUMcols+1;   
         eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
         title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL', ...
            ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,'TONE', ...
            unit.Tone_reBF_simFF.levels_dBSPL(LEVind)),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
         xlim(XLIMITS_dft)
         PLOThand=eval(['h' num2str(PLOTnum)]);
         set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
         ylim(YLIMITS)  % Same Ylimits for all plots
         hold off
         
         %%%% Rate Plot
         PLOTnum=(ROWind-1)*NUMcols+2;   
         eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
         
         %%%% Synch Plot
         PLOTnum=(ROWind-1)*NUMcols+3;   
         eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
         
         %%%% Phase Plot
         PLOTnum=(ROWind-1)*NUMcols+4;   
         eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      end   
      
      
      %%%%%%%%%%%%%%%%%%%% EH_reBF Plots
      if isfield(unit,'EH_reBF_simFF')
         for FeatIND=FeatINDs
            ROWind=ROWind+1;
            for HarmonicsIND=1:1
               for PolarityIND=1:1
                  eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
                  if ~isempty(yTEMP)
                     %%%% EH_reBF plots
                     LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                     
                     %%%% Spatio-Temporal Plots
                     PLOTnum=(ROWind-1)*NUMcols+1;
                     eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                     for BFind=1:length(BFs_kHz{ROWind,LEVind})
                        if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                           LINEwidth=2;
                        else
                           LINEwidth=.5;
                        end
                        % This normalization plots each signal the same size on a log scale
                        if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                           NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/DFTsMAX;
                           %                         plot(DFTfreqs_Hz{ROWind,LEVind}{BFind}/1000, ...
                           %                            abs(DFTs{ROWind,LEVind}{BFind})*NormFact+BFs_kHz{ROWind,LEVind}(BFind),'b-x', ...
                           %                            'LineWidth',LINEwidth)
                           semilogy(DFTfreqs_Hz{ROWind,LEVind}{BFind}/1000, ...
                              abs(DFTs{ROWind,LEVind}{BFind})*NormFact+BFs_kHz{ROWind,LEVind}(BFind),'b-x', ...
                              'LineWidth',LINEwidth)
                           hold on
                        end
                     end
                     %                   plot([1e-6 1e6],[1e-6 1e6],'k')
                     semilogy(BFs_kHz{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'k')
                     xlabel('Stimulus Frequency (kHz)')
                     ylabel('Effective Best Frequency (kHz)')                  
                     if ROWind==1
                        title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                           ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                           yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                     else
                        title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                           yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                     end
                     xlim(XLIMITS_dft)
                     PLOThand=eval(['h' num2str(PLOTnum)]);
                     set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                     ylim(YLIMITS)  % Same Ylimits for all plots
                     %%%%%%%%%%%%%%%%%%%%%
                     % Plot lines at all features
                     for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                        if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                           semilogy(XLIMITS_dft,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                           text(XLIMITS_dft(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                              sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                              'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        end
                        semilogy(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],YLIMITS,':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        text(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000,YLIMITS(1)*1.0, ...
                           sprintf('%s',FeaturesText{FeatINDPlot}),'units','data','HorizontalAlignment','center','VerticalAlignment','top', ...
                           'Color',FeatureColors{-rem(FeatINDPlot,2)+2},'FontSize',6)
                        %                      for BFind=1:length(BFs_kHz{ROWind,LEVind})
                        %                         if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        %                            if (FeatINDPlot<=FeatIND)
                        %                               if (FeatINDPlot>1)
                        %                                  text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                        %                                     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        %                               else
                        %                                  text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                        %                                     'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                        %                               end
                        %                            end
                        %                         end
                        %                      end
                     end
                     hold off
                     
                     
                     %%%% Rate Plot
                     PLOTnum=(ROWind-1)*NUMcols+2;   
                     eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                     semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                     hold on
                     semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
                     semilogy(ALSRs{ROWind,LEVind},unit.Info.BF_kHz,'go','MarkerSize',6)
                     semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                     xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'))
                     PLOThand=eval(['h' num2str(PLOTnum)]);
                     xlim(XLIMITS_rate)
                     set(PLOThand,'XDir','reverse')
                     set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                     ylim(YLIMITS)  % Same Ylimits for all plots
                     %%%%%%%%%%%%%%%%%%%%%
                     % Plot lines at all features
                     for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                        if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                           semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        end
                     end
                     hold off
                     
                     %%%% Synch Plot
                     PLOTnum=(ROWind-1)*NUMcols+3;   
                     eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                     semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                     hold on
                     semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                     xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
                     PLOThand=eval(['h' num2str(PLOTnum)]);
                     xlim(XLIMITS_synch)
                     set(PLOThand,'XDir','reverse')
                     set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                     set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
                     ylim(YLIMITS)  % Same Ylimits for all plots
                     %%%%%%%%%%%%%%%%%%%%%
                     % Plot lines at all features
                     for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                        if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                           semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        end
                     end
                     hold off
                     
                     %%%% Phase Plot
                     PLOTnum=(ROWind-1)*NUMcols+4;   
                     eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                     semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                     hold on
                     semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                     xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
                     PLOThand=eval(['h' num2str(PLOTnum)]);
                     xlim(XLIMITS_phase)
                     set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
                     set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                     ylim(YLIMITS)  % Same Ylimits for all plots
                     %%%%%%%%%%%%%%%%%%%%%
                     % Plot lines at all features
                     for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                        if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                           semilogy(XLIMITS_phase,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        end
                     end
                     hold off
                     
                  end %End if data for this condition, plot
               end % End Polarity
            end % End Harmonics
         end % End Feature
      end
      
      
      Xcorner=0.05;
      Xwidth1=.5;
      Xshift1=0.05;
      Xwidth2=.1;
      Xshift2=0.03;
      
      Ycorner=0.05;
      Yshift=0.07;
      Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
      
      TICKlength=0.02;
      
      if NUMrows>4
         set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
         set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      end
      
      if NUMrows>3
         set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
         set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      end
      
      if NUMrows>2
         set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
         set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      end
      
      if NUMrows>1
         set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
         set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
         set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      end
      
      set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      
      orient landscape
   end
   
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%% ADD EXTRA DFT PLOTS with ALL levels
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   LEVEL=max(levels_dBSPL);
   figure(1001); clf
   set(gcf,'pos',[420     4   977   976])
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Plots
   %
   %%%% SKIP TONE DFTs!!!!!!!!!!!!!!!!  
   %
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
   end   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Plots
   if isfield(unit,'EH_reBF_simFF')
      for FeatIND=FeatINDs
         ROWind=ROWind+1;
         for HarmonicsIND=1:1
            for PolarityIND=1:1
               eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
               if ~isempty(yTEMP)
                  %%%% EH_reBF plots
                  LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                  
                  %%%% Spatio-Temporal Plots
                  PLOTnum=(ROWind-1)*NUMcols+1;
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  LEGtext='';
                  for BFind=1:length(BFs_kHz{ROWind,LEVind})
                     if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                        LINEwidth=2;
                     else
                        LINEwidth=.5;
                     end
                     % This normalization plots each signal the same size on a log scale
                     if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/DFTsMAX;
                        %                      plot(DFTfreqs_Hz{ROWind,LEVind}{BFind}/1000, ...
                        %                         abs(DFTs{ROWind,LEVind}{BFind})*NormFact+BFs_kHz{ROWind,LEVind}(BFind),'-x', ...
                        %                         'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
                        semilogy(DFTfreqs_Hz{ROWind,LEVind}{BFind}/1000, ...
                           abs(DFTs{ROWind,LEVind}{BFind})*NormFact+BFs_kHz{ROWind,LEVind}(BFind),'-x', ...
                           'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
                        hold on
                        if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                           LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
                        end
                        for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                           if ~isempty(PERhistXs_sec{ROWind,LEVind2}{BFind})                             
                              %                            plot(DFTfreqs_Hz{ROWind,LEVind2}{BFind}/1000, ...
                              %                               abs(DFTs{ROWind,LEVind2}{BFind})*NormFact+BFs_kHz{ROWind,LEVind2}(BFind),'-x', ...
                              %                               'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
                              semilogy(DFTfreqs_Hz{ROWind,LEVind2}{BFind}/1000, ...
                                 abs(DFTs{ROWind,LEVind2}{BFind})*NormFact+BFs_kHz{ROWind,LEVind2}(BFind),'-x', ...
                                 'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
                              if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                                 LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind2));
                              end
                           end
                        end
                        if strcmp(FeaturesText{FeatIND},'F1')
                           hleg1001=legend(LEGtext,1);
                           set(hleg1001,'FontSize',8)
                           set(hleg1001,'pos',[0.4451    0.8913    0.0942    0.0473])
                        end
                     end
                  end
                  %                plot([1e-6 1e6],[1e-6 1e6],'k')
                  semilogy(BFs_kHz{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'k')
                  xlabel('Stimulus Frequency (kHz)')
                  ylabel('Effective Best Frequency (kHz)')                  
                  if ROWind==1
                     title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                        ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  else
                     title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  end
                  xlim(XLIMITS_dft)
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_dft,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        text(XLIMITS_dft(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                           sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                           'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                     semilogy(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],YLIMITS,':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     text(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000,YLIMITS(1)*1.0, ...
                        sprintf('%s',FeaturesText{FeatINDPlot}),'units','data','HorizontalAlignment','center','VerticalAlignment','top', ...
                        'Color',FeatureColors{-rem(FeatINDPlot,2)+2},'FontSize',6)
                  end
                  hold off
                  
                  
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end
   
   
   Xcorner=0.05;
   Xwidth1=.5;
   Xshift1=0.05;
   Xwidth2=.1;
   Xshift2=0.03;
   
   Ycorner=0.05;
   Yshift=0.07;
   Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
   
   TICKlength=0.02;
   
   if NUMrows>4
      set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>3
      set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>2
      set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>1
      set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      %       set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   %    set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   %    set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   %    set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   
   orient landscape
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL SCC PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XLIMITS_scc=1*[-PERhist_XMAX PERhist_XMAX];
sccGAIN=5.0; % # of channels covered by max PERhist
scc_logCHwidth=PERhists_logCHwidth;  % log10 channel width

for LEVEL=levels_dBSPL
   % beep
   % disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
   % for LEVEL=levels_dBSPL(end)
   figure(round(LEVEL)+2); clf
   set(gcf,'pos',[420     4   977   976])
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Plots
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      
      PLOTnum=(ROWind-1)*NUMcols+1;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      
      for SCCind=1:length(NSCCs{ROWind,LEVind})
         %                      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
         %                         LINEwidth=2;
         %                      else
         LINEwidth=.5;
         %                      end
         % This normalization plots each signal the same size on a log scale
         if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
            NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind})/SCCsMAX;
            semilogy(NSCC_delays_usec{ROWind,LEVind}{SCCind}/1000, ...
               NSCCs{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
               'LineWidth',LINEwidth)
            hold on
            % Plot 0 and 1 lines
            semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
               'k-','LineWidth',LINEwidth/2)
            semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
               'k-','LineWidth',LINEwidth/2)
            semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
            % Plot values at 0 delay and CD
            semilogy(0,NSCC_0delay{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
               'sk')
            semilogy(NSCC_CDs_usec{ROWind,LEVind}{SCCind}/1000,NSCC_peaks{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
               'or')
            
         end
      end
      xlabel('Delay (ms)')
      ylabel('Effective Best Frequency (kHz)')
      title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL', ...
         ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,'TONE', ...
         unit.Tone_reBF_simFF.levels_dBSPL(LEVind)),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
      xlim(XLIMITS_scc)
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      %%%%%%%%%%%%%%%%%%%%%
      % Plot 1 Period of tone
      text(1000/(unit.Info.BF_kHz*1000),YLIMITS(1),sprintf('1/f'),'units','data', ...
         'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','b')
      hold off
      
            
      %%%% Rate Plot
      PLOTnum=(ROWind-1)*NUMcols+2;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]'))
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_rate)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Synch Plot
      PLOTnum=(ROWind-1)*NUMcols+3;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Synch Coef (to f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_synch)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Phase Plot
      PLOTnum=(ROWind-1)*NUMcols+4;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Phase (cycles of f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_phase)
      set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
   end   
   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Plots
   if isfield(unit,'EH_reBF_simFF')
      for FeatIND=FeatINDs
         ROWind=ROWind+1;
         for HarmonicsIND=1:1
            for PolarityIND=1:1
               eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
               if ~isempty(yTEMP)
                  %%%% EH_reBF plots
                  LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                  
                  %%%% Spatio-Temporal Plots
                  PLOTnum=(ROWind-1)*NUMcols+1;
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  for SCCind=1:length(NSCCs{ROWind,LEVind})
                     %                      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                     %                         LINEwidth=2;
                     %                      else
                     LINEwidth=.5;
                     %                      end
                     % This normalization plots each signal the same size on a log scale
                     if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
                        NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind})/SCCsMAX;
                        semilogy(NSCC_delays_usec{ROWind,LEVind}{SCCind}/1000, ...
                           NSCCs{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                           'LineWidth',LINEwidth)
                        hold on
                        % Plot 0 and 1 lines
                        semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                           'k-','LineWidth',LINEwidth/2)
                        semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                           'k-','LineWidth',LINEwidth/2)
                        semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
                        % Plot values at 0 delay and CD
                        semilogy(0,NSCC_0delay{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                           'sk')
                        semilogy(NSCC_CDs_usec{ROWind,LEVind}{SCCind}/1000,NSCC_peaks{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                           'or')
                        
                     end
                  end
                  xlabel('Delay (ms)')
                  ylabel('Effective Best Frequency (kHz)')
                  if ROWind==1
                     title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                        ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  else
                     title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  end
                  xlim(XLIMITS_scc)
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     %                      if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                     %                         semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     %                         text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                     %                            sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                     %                            'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     %                      end
                     for SCCind=1:length(NSCCs{ROWind,LEVind})
                        if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
                           if (FeatINDPlot<=FeatIND)
                              if (FeatINDPlot>1)
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                              else
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                              end
                           end
                        end
                     end
                  end
                  hold off
                  
                  
                  %%%% Rate Plot
                  PLOTnum=(ROWind-1)*NUMcols+2;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
                  semilogy(ALSRs{ROWind,LEVind},unit.Info.BF_kHz,'go','MarkerSize',6)
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]\nO: ALSR'))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_rate)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Synch Plot
                  PLOTnum=(ROWind-1)*NUMcols+3;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_synch)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Phase Plot
                  PLOTnum=(ROWind-1)*NUMcols+4;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_phase)
                  set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_phase,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end
   
   
   Xcorner=0.05;
   Xwidth1=.5;
   Xshift1=0.05;
   Xwidth2=.1;
   Xshift2=0.03;
   
   Ycorner=0.05;
   Yshift=0.07;
   Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
   
   TICKlength=0.02;
   
   if NUMrows>4
      set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>3
      set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>2
      set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>1
      set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   
   orient landscape
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% EXTRA SCC PLOT WITH ALL LEVElS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XLIMITS_scc=4*[-1 1];
LEVEL=max(levels_dBSPL);
figure(1002); clf
set(gcf,'pos',[420     4   977   976])
ROWind=0;
   
%%%%%%%%%%%%%%%%%%%% Tone Plots
if isfield(unit,'Tone_reBF_simFF')
   ROWind=ROWind+1;
   
   %%%% Tone_reBF plots
   LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
   
   PLOTnum=(ROWind-1)*NUMcols+1;   
   eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
   LEGtext='';
   for SCCind=1:length(NSCCs{ROWind,LEVind})
      %                      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
      %                         LINEwidth=2;
      %                      else
      LINEwidth=.5;
      %                      end
      % This normalization plots each signal the same size on a log scale
      if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
         NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind})/SCCsMAX;
         semilogy(NSCC_delays_usec{ROWind,LEVind}{SCCind}/1000, ...
            NSCCs{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
            'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
         hold on
         if SCCind==1
            LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
         end
      end
      for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
         % This normalization plots each signal the same size on a log scale
         if ~isempty(NSCCs{ROWind,LEVind2}{SCCind})
            NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind})/SCCsMAX;
            semilogy(NSCC_delays_usec{ROWind,LEVind2}{SCCind}/1000, ...
               NSCCs{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
               'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
            hold on
            if SCCind==1
               LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind2));
            end
         end
      end
      if SCCind==1
         hleg1001a=legend(LEGtext,1);
         set(hleg1001a,'FontSize',8)
         set(hleg1001a,'pos',[0.5700    0.8944    0.0942    0.0473])
      end
      if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
         % Plot 0 and 1 lines
         semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
            'k-','LineWidth',LINEwidth/2)
         semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
            'k-','LineWidth',LINEwidth/2)
         semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
         %          % Plot values at 0 delay and CD
         semilogy(0,NSCC_0delay{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
            's','Color',LEVELcolors{LEVind})
         semilogy(NSCC_CDs_usec{ROWind,LEVind}{SCCind}/1000,NSCC_peaks{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
            'o','Color',LEVELcolors{LEVind})
      end
      for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
         % Plot 0 and 1 lines
         semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
            'k-','LineWidth',LINEwidth/2)
         semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
            'k-','LineWidth',LINEwidth/2)
         semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
         %             % Plot values at 0 delay and CD
         semilogy(0,NSCC_0delay{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
            's','Color',LEVELcolors{LEVind2})
         semilogy(NSCC_CDs_usec{ROWind,LEVind2}{SCCind}/1000,NSCC_peaks{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
            'o','Color',LEVELcolors{LEVind2})
      end
   end
   xlabel('Delay (ms)')
   ylabel('Effective Best Frequency (kHz)')
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s', ...
      ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,'TONE'), ...
      'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
   xlim(XLIMITS_scc)
   PLOThand=eval(['h' num2str(PLOTnum)]);
   set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
   ylim(YLIMITS)  % Same Ylimits for all plots
   %%%%%%%%%%%%%%%%%%%%%
   % Plot 1 Period of tone
   text(1000/(unit.Info.BF_kHz*1000),YLIMITS(1),sprintf('1/f'),'units','data', ...
      'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color','b')
   hold off
   
end   


%%%%%%%%%%%%%%%%%%%% EH_reBF Plots
if isfield(unit,'EH_reBF_simFF')
   for FeatIND=FeatINDs
      ROWind=ROWind+1;
      for HarmonicsIND=1:1
         for PolarityIND=1:1
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               %%%% EH_reBF plots
               LEVind=find(yTEMP.levels_dBSPL==LEVEL);
               
               %%%% Spatio-Temporal Plots
               PLOTnum=(ROWind-1)*NUMcols+1;
               eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
               LEGtext='';
               for SCCind=1:length(NSCCs{ROWind,LEVind})
                  %                      if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                  %                         LINEwidth=2;
                  %                      else
                  LINEwidth=.5;
                  %                      end
                  % This normalization plots each signal the same size on a log scale
                  if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
                     NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind})/SCCsMAX;
                     semilogy(NSCC_delays_usec{ROWind,LEVind}{SCCind}/1000, ...
                        NSCCs{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                        'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind})
                     hold on
                     if SCCind==1
                        LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
                     end
                  end
                  for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                     % This normalization plots each signal the same size on a log scale
                     if ~isempty(NSCCs{ROWind,LEVind2}{SCCind})
                        NormFact=(10^(sccGAIN*scc_logCHwidth)-1)*geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind})/SCCsMAX;
                        semilogy(NSCC_delays_usec{ROWind,LEVind2}{SCCind}/1000, ...
                           NSCCs{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
                           'LineWidth',LINEwidth,'Color',LEVELcolors{LEVind2})
                        hold on
                        if SCCind==1
                           LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',levels_dBSPL(LEVind2));
                        end
                     end
                  end
                  if SCCind==1
                     hleg1001=legend(LEGtext,1);
                     set(hleg1001,'FontSize',8)
                     set(hleg1001,'pos',[0.5700    0.8944    0.0942    0.0473])
                  end
                  if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
                     % Plot 0 and 1 lines
                     semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                        'k-','LineWidth',LINEwidth/2)
                     semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                        'k-','LineWidth',LINEwidth/2)
                     semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
                     %                      % Plot values at 0 delay and CD
                     semilogy(0,NSCC_0delay{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                        's','Color',LEVELcolors{LEVind})
                     semilogy(NSCC_CDs_usec{ROWind,LEVind}{SCCind}/1000,NSCC_peaks{ROWind,LEVind}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind}{SCCind}), ...
                        'o','Color',LEVELcolors{LEVind})                     
                  end
                  for LEVind2=fliplr(find(levels_dBSPL~=max(levels_dBSPL)))
                     % Plot 0 and 1 lines
                     semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
                        'k-','LineWidth',LINEwidth/2)
                     semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
                        'k-','LineWidth',LINEwidth/2)
                     semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
                     %             % Plot values at 0 delay and CD
                     semilogy(0,NSCC_0delay{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
                        's','Color',LEVELcolors{LEVind2})
                     semilogy(NSCC_CDs_usec{ROWind,LEVind2}{SCCind}/1000,NSCC_peaks{ROWind,LEVind2}{SCCind}*NormFact+geomean(NSCC_BFskHz{ROWind,LEVind2}{SCCind}), ...
                        'o','Color',LEVELcolors{LEVind2})
                  end
                  
               end
               xlabel('Delay (ms)')
               ylabel('Effective Best Frequency (kHz)')
               if ROWind==1
                  title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s,   Harm: %d, Polarity: %d', ...
                     ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                     HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
               else
                  title(sprintf('%s,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                     HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
               end
               xlim(XLIMITS_scc)
               PLOThand=eval(['h' num2str(PLOTnum)]);
               set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
               ylim(YLIMITS)  % Same Ylimits for all plots
               %%%%%%%%%%%%%%%%%%%%%
               % Plot lines at all features
               for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
                  %                      if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                  %                         semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  %                         text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                  %                            sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                  %                            'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                  %                      end
                  for SCCind=1:length(NSCCs{ROWind,LEVind})
                     if ~isempty(NSCCs{ROWind,LEVind}{SCCind})
                        if (FeatINDPlot<=FeatIND)
                           if (FeatINDPlot>1)
                              text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                           else
                              text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                           end
                        end
                     end
                  end
               end
               hold off
               
            end %End if data for this condition, plot
         end % End Polarity
      end % End Harmonics
   end % End Feature
end


Xcorner=0.05;
Xwidth1=.5;
Xshift1=0.05;
Xwidth2=.1;
Xshift2=0.03;

Ycorner=0.05;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

TICKlength=0.02;

if NUMrows>4
   set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>3
   set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>2
   set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>1
   set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
%    set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
% set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

orient landscape




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO SMP PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SMPlegFontSize=7;
if isfield(unit,'EH_reBF_simFF')
   FIGnum=1;
   figure(FIGnum); clf
   set(gcf,'pos',[753     4   643   976])

   
   %%%%%%%%%%%%%%%%%%%%% CALCULATE SLOPES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Rate vs. Feature Level
   SensitivitySlopes_rate=NaN+ones(size(levels_dBSPL));
   SensitivityIntercepts_rate=NaN+ones(size(levels_dBSPL));
   for LEVind=1:length(levels_dBSPL)
      if sum(~isnan(SMP_rate{LEVind}))>1
         x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_rate{LEVind})));
         y=SMP_rate{LEVind}(find(~isnan(SMP_rate{LEVind})));
         [Cfit,MSE,fit]=fit1slope(x,y);
         SensitivitySlopes_rate(LEVind)=Cfit(1);
         SensitivityIntercepts_rate(LEVind)=Cfit(2);
      end
   end

   % ALSR vs. Feature Level
   SensitivitySlopes_alsr=NaN+ones(size(levels_dBSPL));
   SensitivityIntercepts_alsr=NaN+ones(size(levels_dBSPL));
   for LEVind=1:length(levels_dBSPL)
      if sum(~isnan(SMP_alsr{LEVind}))>1
         x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_alsr{LEVind})));
         y=SMP_alsr{LEVind}(find(~isnan(SMP_alsr{LEVind})));
         [Cfit,MSE,fit]=fit1slope(x,y);
         SensitivitySlopes_alsr(LEVind)=Cfit(1);
         SensitivityIntercepts_alsr(LEVind)=Cfit(2);
      end
   end

   % NSCC_CD(1) vs. Feature Level
   SensitivitySlopes_nsccCD{1}=NaN+ones(size(levels_dBSPL));
   SensitivityIntercepts_nsccCD{1}=NaN+ones(size(levels_dBSPL));
   for LEVind=1:length(levels_dBSPL)
      if sum(~isnan(SMP_NSCC_CD{1,LEVind}))>1
         x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_NSCC_CD{1,LEVind})));
         y=SMP_NSCC_CD{1,LEVind}(find(~isnan(SMP_NSCC_CD{1,LEVind})));
         [Cfit,MSE,fit]=fit1slope(x,y);
         SensitivitySlopes_nsccCD{1}(LEVind)=Cfit(1);
         SensitivityIntercepts_nsccCD{1}(LEVind)=Cfit(2);
      end
   end

   % NSCC_0delay(1) vs. Feature Level
   SensitivitySlopes_nscc0{1}=NaN+ones(size(levels_dBSPL));
   SensitivityIntercepts_nscc0{1}=NaN+ones(size(levels_dBSPL));
   for LEVind=1:length(levels_dBSPL)
      if sum(~isnan(SMP_NSCC_0delay{1,LEVind}))>1
         x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_NSCC_0delay{1,LEVind})));
         y=SMP_NSCC_0delay{1,LEVind}(find(~isnan(SMP_NSCC_0delay{1,LEVind})));
         [Cfit,MSE,fit]=fit1slope(x,y);
         SensitivitySlopes_nscc0{1}(LEVind)=Cfit(1);
         SensitivityIntercepts_nscc0{1}(LEVind)=Cfit(2);
      end
   end

   % NSCC_peak(1) vs. Feature Level
   SensitivitySlopes_nsccPEAK{1}=NaN+ones(size(levels_dBSPL));
   SensitivityIntercepts_nsccPEAK{1}=NaN+ones(size(levels_dBSPL));
   for LEVind=1:length(levels_dBSPL)
      if sum(~isnan(SMP_NSCC_peak{1,LEVind}))>1
         x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,LEVind})));
         y=SMP_NSCC_peak{1,LEVind}(find(~isnan(SMP_NSCC_peak{1,LEVind})));
         [Cfit,MSE,fit]=fit1slope(x,y);
         SensitivitySlopes_nsccPEAK{1}(LEVind)=Cfit(1);
         SensitivityIntercepts_nsccPEAK{1}(LEVind)=Cfit(2);
      end
   end


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%% Rate-FeatureLevels for each OAL
   subplot(421)
   LEGtext='';
   for LEVind=1:length(levels_dBSPL)
      %% Plot fitted lines
      if ~isnan(SensitivitySlopes_rate(LEVind))
         xdata=[min(FeatureLevels_dB(find(~isnan(SMP_rate{LEVind})))+levels_dBSPL(LEVind)) max(FeatureLevels_dB(find(~isnan(SMP_rate{LEVind})))+levels_dBSPL(LEVind))];
         ydata=SensitivitySlopes_rate(LEVind)*xdata+SensitivityIntercepts_rate(LEVind);
         plot(xdata,ydata,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
         hold on
         LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.2f sp/sec/dB)',levels_dBSPL(LEVind),SensitivitySlopes_rate(LEVind));
      end
   end
   for LEVind=1:length(levels_dBSPL)
      for FeatIND=FeatINDs
         if ~isnan(SMP_rate{LEVind}(FeatIND))
            plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_rate{LEVind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
            hold on
         end
      end
   end
   YLIMITS_SMPrate=[0 300];
   XLIMITS_SMPrate=[0 100];
   ylim(YLIMITS_SMPrate)  % Fixed ordinate for all plots
   xlim(XLIMITS_SMPrate)
   %%%%%%%%%%%%%%%%%%%%%
   % Home-made Feature Symbol Legend
   %%%%%%%%%%%%%%%%%%%%%
   LEGXleft=.8; LEGYbottom=0.15; LEGYstep=0.1; LEGXsymbOFFset=0.05;
   FeatNum=0;
   for FeatIND=fliplr(FeatINDs)
      FeatNum=FeatNum+1;
      plot(XLIMITS_SMPrate(1)+diff(XLIMITS_SMPrate)*(LEGXleft),YLIMITS_SMPrate(1)+diff(YLIMITS_SMPrate)*(LEGYbottom+(FeatNum-1)*LEGYstep),FEATmarkers{FeatIND}, ...
         'Color','k','MarkerSize',6)
      text(LEGXleft+LEGXsymbOFFset,LEGYbottom+(FeatNum-1)*LEGYstep,FeaturesText{FeatIND},'Units','norm','FontSize',10)
   end
   ylabel('Rate (sp/sec)')
   xlabel('Feature Level (dB SPL)')
   hleg=legend(LEGtext,1);
   set(hleg,'FontSize',SMPlegFontSize)
   hold off
   set(gca,'FontSize',FIG.FontSize)
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
      ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10), ...
      'units','norm')
   

   %%%% ALSR-FeatureLevels for each OAL
   subplot(422)
   LEGtext='';
   for LEVind=1:length(levels_dBSPL)
      %% Plot fitted lines
      if ~isnan(SensitivitySlopes_alsr(LEVind))
         xdata=[min(FeatureLevels_dB(find(~isnan(SMP_alsr{LEVind})))+levels_dBSPL(LEVind)) max(FeatureLevels_dB(find(~isnan(SMP_alsr{LEVind})))+levels_dBSPL(LEVind))];
         ydata=SensitivitySlopes_alsr(LEVind)*xdata+SensitivityIntercepts_alsr(LEVind);
         plot(xdata,ydata,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
         hold on
         LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.2f sp/sec/dB)',levels_dBSPL(LEVind),SensitivitySlopes_alsr(LEVind));
      end
   end
   for LEVind=1:length(levels_dBSPL)
      for FeatIND=FeatINDs
         if ~isnan(SMP_alsr{LEVind}(FeatIND))
            plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_alsr{LEVind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
            hold on
         end
      end
   end
   YLIMITS_SMPalsr=[0 300];
   XLIMITS_SMPalsr=XLIMITS_SMPrate;
   ylim(YLIMITS_SMPalsr)  % Fixed ordinate for all plots
   xlim(XLIMITS_SMPalsr)
   ylabel('ALSR (sp/sec)')
   xlabel('Feature Level (dB SPL)')
   hleg=legend(LEGtext,1);
   set(hleg,'FontSize',SMPlegFontSize)
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% NSCC_0delay-FeatureLevels for each OAL
   subplot(424)
   LEGtext='';
   for LEVind=1:length(levels_dBSPL)
      %% Plot fitted lines
      if ~isnan(SensitivitySlopes_nscc0{1}(LEVind))
         xdata=[min(FeatureLevels_dB(find(~isnan(SMP_NSCC_0delay{1,LEVind})))+levels_dBSPL(LEVind)) max(FeatureLevels_dB(find(~isnan(SMP_NSCC_0delay{1,LEVind})))+levels_dBSPL(LEVind))];
         ydata=SensitivitySlopes_nscc0{1}(LEVind)*xdata+SensitivityIntercepts_nscc0{1}(LEVind);
         plot(xdata,ydata,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
         hold on
         LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.3f 1/dB)',levels_dBSPL(LEVind),SensitivitySlopes_nscc0{1}(LEVind));
      end
   end
   for LEVind=1:length(levels_dBSPL)
      for FeatIND=FeatINDs
         if ~isnan(SMP_NSCC_0delay{1,LEVind}(FeatIND))
            plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_NSCC_0delay{1,LEVind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
            hold on
         end
      end
   end
   YLIMITS_SMPnscc0=[0 5];
   XLIMITS_SMPnscc0=XLIMITS_SMPrate;
   ylim(YLIMITS_SMPnscc0)  % Fixed ordinate for all plots
   xlim(XLIMITS_SMPnscc0)
   ylabel('NSCC (at 0 delay)')
   title(sprintf('NSCC1: BF+%.2f octaves re BF-%.2f octaves',SCC_octOFFSET),'color','red')
   xlabel('Feature Level (dB SPL)')
   hleg=legend(LEGtext,1);
   set(hleg,'FontSize',SMPlegFontSize)
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   
   %%%% NSCC_peak-FeatureLevels for each OAL
   subplot(423)
   LEGtext='';
   for LEVind=1:length(levels_dBSPL)
      %% Plot fitted lines
      if ~isnan(SensitivitySlopes_nsccPEAK{1}(LEVind))
         xdata=[min(FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,LEVind})))+levels_dBSPL(LEVind)) max(FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,LEVind})))+levels_dBSPL(LEVind))];
         ydata=SensitivitySlopes_nsccPEAK{1}(LEVind)*xdata+SensitivityIntercepts_nsccPEAK{1}(LEVind);
         plot(xdata,ydata,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
         hold on
         LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.3f 1/dB)',levels_dBSPL(LEVind),SensitivitySlopes_nsccPEAK{1}(LEVind));
      end
   end
   for LEVind=1:length(levels_dBSPL)
      for FeatIND=FeatINDs
         if ~isnan(SMP_NSCC_peak{1,LEVind}(FeatIND))
            plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_NSCC_peak{1,LEVind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
            hold on
         end
      end
   end
   YLIMITS_SMPnsccPEAK=[0 5];
   XLIMITS_SMPnsccPEAK=XLIMITS_SMPrate;
   ylim(YLIMITS_SMPnsccPEAK)  % Fixed ordinate for all plots
   xlim(XLIMITS_SMPnsccPEAK)
   ylabel('Peak NSCC (at CD)')
   xlabel('Feature Level (dB SPL)')
   hleg=legend(LEGtext,1);
   set(hleg,'FontSize',SMPlegFontSize)
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   
   %%%% NSCC_CD-FeatureLevels for each OAL
   subplot(427)
   LEGtext='';
   for LEVind=1:length(levels_dBSPL)
      %% Plot fitted lines
      if ~isnan(SensitivitySlopes_nsccCD{1}(LEVind))
         xdata=[min(FeatureLevels_dB(find(~isnan(SMP_NSCC_CD{1,LEVind})))+levels_dBSPL(LEVind))  ...
               max(FeatureLevels_dB(find(~isnan(SMP_NSCC_CD{1,LEVind})))+levels_dBSPL(LEVind))];
         ydata=SensitivitySlopes_nsccCD{1}(LEVind)*xdata+SensitivityIntercepts_nsccCD{1}(LEVind);
         plot(xdata,ydata/1000,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
         hold on
         %          LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.2f usec/dB)',levels_dBSPL(LEVind),SensitivitySlopes_nsccCD{1}(LEVind));
      end
   end
   for LEVind=1:length(levels_dBSPL)
      for FeatIND=FeatINDs
         if ~isnan(SMP_NSCC_CD{1,LEVind}(FeatIND))
            plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_NSCC_CD{1,LEVind}(FeatIND)/1000, ...
               'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
            hold on
         end
      end
   end
   YLIMITS_SMPnsccCD=1*[-1 1];
   XLIMITS_SMPnsccCD=XLIMITS_SMPrate;
   ylim(YLIMITS_SMPnsccCD)  % Fixed ordinate for all plots
   xlim(XLIMITS_SMPnsccCD)
   plot(XLIMITS_SMPnsccCD,zeros(1,2),'k:')
   ylabel('Characteristic Dealy of NSCC (msec)')
   xlabel('Feature Level (dB SPL)')
   %    hleg=legend(LEGtext,1);
   %    set(hleg,'FontSize',SMPlegFontSize)
   hold off
   set(gca,'FontSize',FIG.FontSize)
   

   
   %%%%%%%%%%%%%% Interactive Section to pick delays to try SMP on
   
   INTERACTyes=1;
   ARBdelay_msec=NaN;
   
   while INTERACTyes
      
      %%%%%%%%%%%%%%% ADD INTERACTIVE STUFF HERE
      
      temp=input('Enter arbitrary delay to compute SMP (in msec) [Return to quit]: ');
      if isempty(temp)
         INTERACTyes=0;
      else
         ARBdelay_msec=temp;   
         
         
         for LEVEL=levels_dBSPL
            ROWind=0;
            %%%%%%%%%%%%%%%%%%%% Tone Calcs
            if isfield(unit,'Tone_reBF_simFF')
               ROWind=ROWind+1;
            end   % if Tone data
            
            %%%%%%%%%%%%%%%%%%%% EH_reBF Calcs
            if isfield(unit,'EH_reBF_simFF')
               LEVind=find(levels_dBSPL==LEVEL);
               for FeatIND=FeatINDs
                  ROWind=ROWind+1;
                  
                  for SCCind=1:length(NSCC_BFinds)  % index of SCC to calculate
                     if ~isempty(NSCCs{ROWind,LEVind})
                        [yy,ABRdelayIND]=min(abs(NSCC_delays_usec{ROWind,LEVind}{SCCind}-ARBdelay_msec*1000));
                        NSCC_ARBdelay{ROWind,LEVind}{SCCind}=NSCCs{ROWind,LEVind}{SCCind}(ABRdelayIND);
                     else  % No Data to compute SCCs
                        NSCC_ARBdelay{ROWind,LEVind}{SCCind}=NaN;
                     end
                     
                     %%%%%%%%%%%%%%%%
                     % Store SMP data
                     %%%%%%%%%%%%%%%%
                     SMP_NSCC_ARBdelay{SCCind,LEVind}(FeatIND)=NSCC_ARBdelay{ROWind,LEVind}{i};
                     
                  end %End if data for this condition, plot
               end % End Feature
            end % If EHrBF data
         end % Levels
         
         % NSCC_ARBdelay(1) vs. Feature Level
         SensitivitySlopes_nsccARB{1}=NaN+ones(size(levels_dBSPL));
         SensitivityIntercepts_nsccARB{1}=NaN+ones(size(levels_dBSPL));
         for LEVind=1:length(levels_dBSPL)
            if sum(~isnan(SMP_NSCC_ARBdelay{1,LEVind}))>1
               x=levels_dBSPL(LEVind)+FeatureLevels_dB(find(~isnan(SMP_NSCC_ARBdelay{1,LEVind})));
               y=SMP_NSCC_ARBdelay{1,LEVind}(find(~isnan(SMP_NSCC_ARBdelay{1,LEVind})));
               [Cfit,MSE,fit]=fit1slope(x,y);
               SensitivitySlopes_nsccARB{1}(LEVind)=Cfit(1);
               SensitivityIntercepts_nsccARB{1}(LEVind)=Cfit(2);
            end
         end
         
         %%%% NSCC_ARBdelay-FeatureLevels for each OAL
         figure(FIGnum)
         subplot(425)
         LEGtext='';
         for LEVind=1:length(levels_dBSPL)
            %% Plot fitted lines
            if ~isnan(SensitivitySlopes_nsccARB{1}(LEVind))
               xdata=[min(FeatureLevels_dB(find(~isnan(SMP_NSCC_ARBdelay{1,LEVind})))+levels_dBSPL(LEVind)) max(FeatureLevels_dB(find(~isnan(SMP_NSCC_ARBdelay{1,LEVind})))+levels_dBSPL(LEVind))];
               ydata=SensitivitySlopes_nsccARB{1}(LEVind)*xdata+SensitivityIntercepts_nsccARB{1}(LEVind);
               plot(xdata,ydata,'Marker','none','Color',LEVELcolors{LEVind},'LineStyle','-')
               hold on
               LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (%.3f 1/dB)',levels_dBSPL(LEVind),SensitivitySlopes_nsccARB{1}(LEVind));
            end
         end
         for LEVind=1:length(levels_dBSPL)
            for FeatIND=FeatINDs
               if ~isnan(SMP_NSCC_ARBdelay{1,LEVind}(FeatIND))
                  plot(levels_dBSPL(LEVind)+FeatureLevels_dB(FeatIND),SMP_NSCC_ARBdelay{1,LEVind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',LEVELcolors{LEVind},'LineStyle','none')
                  hold on
               end
            end
         end
         YLIMITS_SMPnsccARB=YLIMITS_SMPnscc0;
         XLIMITS_SMPnsccARB=XLIMITS_SMPnscc0;
         ylim(YLIMITS_SMPnsccARB)  % Fixed ordinate for all plots
         xlim(XLIMITS_SMPnsccARB)
         ylabel(sprintf('NSCC (at delay = %.2f msec)',ARBdelay_msec))
         title(sprintf('NSCC1: BF+%.2f octaves re BF-%.2f octaves',SCC_octOFFSET,SCC_octOFFSET),'color','red')
         xlabel('Feature Level (dB SPL)')
         hleg=legend(LEGtext,1);
         set(hleg,'FontSize',SMPlegFontSize)
         hold off
         set(gca,'FontSize',FIG.FontSize)
      end
      
      
   end % While
      
      
end %End if data for this condition, plot


orient tall



% Turn off saved PICS feature
SavedPICS=[]; SavedPICnums=[];
SavedPICSuse=0;

% % For ARO2005: ANmodel figure
% save SMP_FIG

return;

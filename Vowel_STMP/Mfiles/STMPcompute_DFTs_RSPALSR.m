function STMPcompute_DFTs_RSPALSR(ExpDate,UnitName,STIMtype,DATAtype)
% M. Heinz Jun 10, 2007
% FROM: STMPcompute_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% FROM: STMPconvert_STs_Nchans1stim(ExpDate,UnitName,STIMtype)
%
% Computes DFTs (from PERhistograms) for a STMP POPULATION, and then
% computes: Rate, Synch, Phase, and ALSR for a STMP population.  Assumes
% STMP format for  data storage, but this can be either: DATAtype= 
%     - Nchans1stim (PHs.Info.STMPshift=1)
%     - 1unitNstim (PHs.Info.STMPshift=0)
%
% NOTE: Synchronized Rates (sps) = abs(DFT);
%
% Setup to be general for STMP data.
%
% Calls basic function: DFT_PERhist, ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
global FeaturesText

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=5;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03), 
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
   elseif TESTunitNUM==5
      ExpDate='041805'; UnitName='2.04'; STIMtype='TrBFi'; DATAtype='Nchans1stim';  % TONE re BFI, with Tone RLV
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
RAWdata_dir=fullfile(STMP_dir,'ExpData',ExpName);
eval(['cd ''' RAWdata_dir ''''])
UNITinfo_dir=fullfile(STMP_dir,'ExpData',ExpName,'UNITinfo');   % For general unit info (BF, SR, bad lines, ...)
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify PHs s exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHs_filename=sprintf('PHs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
PHs_dataname=sprintf('PHs_%s',DATAtype);
eval(['load ''' fullfile(STMPanal_dir,PHs_filename) ''''])
eval(['PHs = ' PHs_dataname ';']);
eval(['clear ' PHs_dataname]);

disp(sprintf('... STMP: Computing ''%s'' DFTs for: Experiment: ''%s''; Unit: %d.%02d',DATAtype,ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% See if PHs_filename exists already
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFTs_filename=sprintf('DFTs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
DFTs_dataname=sprintf('DFTs_%s',DATAtype);
eval(['ddd=dir(''' fullfile(STMPanal_dir,DFTs_filename) ''');'])
if ~isempty(ddd)
   beep
   TEMP = input(sprintf('File: ''%s'' already exists!!\n  ***** Do you want to re-compute DFTs, or leave as is?  [0]: LEAVE AS IS; RERUN: 1;  ',DFTs_filename));
   if isempty(TEMP)
      TEMP=0;
   end
   if TEMP~=1
      beep
      disp(sprintf(' FILE NOT ALTERED\n'));
      return;
   else
      disp(' ... Re-Running PERhists - SAVING NEW FILE!');
   end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup PERhists data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFTs=PHs; 
DFTs=rmfield(DFTs,'PERhists');

%              Info: [1x1 struct]             % UNCHANGED
%          STIMtype: 'EHrBFi'                 % UNCHANGED
%       ChannelInfo: [1x1 struct]             % UNCHANGED
%     ParameterInfo: [1x1 struct]             % UNCHANGED
%     ConditionInfo: [1x1 struct]             % UNCHANGED
%          StimInfo: [1x1 struct]             % UNCHANGED
%          PERhists: {2x1 cell}               %  USED AND DISCARDED
%              DFTs: {2x1 cell}               % COMPUTED 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this unit, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(PHs.PERhists.PERhists{1,1,1},1);
NumP=size(PHs.PERhists.PERhists{1,1,1},2);
NumFEATURES=size(PHs.PERhists.PERhists,1);
NumPOL=size(PHs.PERhists.PERhists,2);
NumHARMS=size(PHs.PERhists.PERhists,3);

DFTs.DFTs.DFTs=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.DFTs.DFT_freqbinwidth_Hz=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.DFTs.PH_orig_NumDrivenSpikes = PHs.PERhists.PH_orig_NumDrivenSpikes;

DFTs.Feature_RSP.rate_sps=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.Feature_RSP.synch=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.Feature_RSP.phase_cyc=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.Feature_RSP.Rayleigh_P=cell(NumFEATURES,NumPOL,NumHARMS);

DFTs.F0_RSP.rate_sps=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.F0_RSP.synch=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.F0_RSP.phase_cyc=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.F0_RSP.Rayleigh_P=cell(NumFEATURES,NumPOL,NumHARMS);

DFTs.ALSR.ALSRs_sps=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.ALSR.INDsUsed=cell(NumFEATURES,NumPOL,NumHARMS);
DFTs.ALSR.oct_range=cell(NumFEATURES,NumPOL,NumHARMS);

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         DFTs.DFTs.DFTs{FeatIND,PolIND,HarmIND}=cell(NumCH,NumP);
         DFTs.DFTs.DFT_freqbinwidth_Hz{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);

         DFTs.Feature_RSP.rate_sps{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.Feature_RSP.synch{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.Feature_RSP.phase_cyc{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.Feature_RSP.Rayleigh_P{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);

         DFTs.F0_RSP.rate_sps{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.F0_RSP.synch{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.F0_RSP.phase_cyc{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.F0_RSP.Rayleigh_P{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);

         DFTs.ALSR.ALSRs_sps{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
         DFTs.ALSR.INDsUsed{FeatIND,PolIND,HarmIND}=cell(NumCH,NumP);
         DFTs.ALSR.oct_range{FeatIND,PolIND,HarmIND}=NaN+zeros(NumCH,NumP);
       
         for ParamIND=1:NumP
            for ChanIND=1:NumCH
               
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%% Compute DFT from PERhists to obtain Synchronized Rate functions               
               PERhist=PHs.PERhists.PERhists{FeatIND,PolIND,HarmIND}{ChanIND,ParamIND};
               PH_binwidth_sec=PHs.PERhists.PH_binwidth_sec{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND);
               
               [SynchRate_sps,DFT_freqbinwidth_Hz]=SynchRate_PERhist(PERhist,PH_binwidth_sec);

               DFTs.DFTs.DFTs{FeatIND,PolIND,HarmIND}{ChanIND,ParamIND}=SynchRate_sps;
               DFTs.DFTs.DFT_freqbinwidth_Hz{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=DFT_freqbinwidth_Hz;

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%%%% Compute rate, synch, and phase at a given frequency from Synchronized Rate functions
               %% NOTE: non-significant synch and phase returned as NaN
               %
               %% First for FEATURE frequency of interest                
               if strcmp(STIMtype,'TrBFi')
                  freq = PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{ChanIND};
               else
                  REALfeatIND=find(strcmp(FeaturesText,PHs.ConditionInfo.features{FeatIND}));
                  freq = PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{ChanIND}(REALfeatIND);
               end
               % Need to take ORIGINAL numspikes for STMP-shifted
               % conditions for Rayleigh stats
               NumSpikes = PHs.PERhists.PH_orig_NumDrivenSpikes{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND);

               [rate_sps,synch,phase_cyc,Rayleigh_P]=RSP_SynchRate(SynchRate_sps,freq,DFT_freqbinwidth_Hz,NumSpikes);
               
%                if FeatIND==1 & ChanIND==6
%                   input('Enter to cont')
%                end
%                
               DFTs.Feature_RSP.rate_sps{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=rate_sps;
               DFTs.Feature_RSP.synch{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=synch;
               DFTs.Feature_RSP.phase_cyc{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=phase_cyc;
               DFTs.Feature_RSP.Rayleigh_P{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=Rayleigh_P;

               %% Second for fundamental frequency (F0) 
               freq = PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{ChanIND}(find(strcmp(FeaturesText,'F0')));
               % Need to take ORIGINAL numspikes for STMP-shifted
               % conditions for Rayleigh stats
               NumSpikes = PHs.PERhists.PH_orig_NumDrivenSpikes{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND);

               [rate_sps,synch,phase_cyc,Rayleigh_P]=RSP_SynchRate(SynchRate_sps,freq,DFT_freqbinwidth_Hz,NumSpikes);
               
               DFTs.F0_RSP.rate_sps{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=rate_sps;
               DFTs.F0_RSP.synch{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=synch;
               DFTs.F0_RSP.phase_cyc{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=phase_cyc;
               DFTs.F0_RSP.Rayleigh_P{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=Rayleigh_P;

            end  % end Chan

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Compute ALSRs if this is POP data
            % ALSR is the average SynchRates_sps over a specified +- OCT_range
            if strcmp(DATAtype,'Nchans1stim')
               CenterFreq=PHs.Info.BF_kHz;
               chan_BFs=PHs.ChannelInfo.channel_values{FeatIND,PolIND,HarmIND};
               SynchRates_sps=DFTs.Feature_RSP.rate_sps{FeatIND,PolIND,HarmIND}(:,ParamIND).*DFTs.Feature_RSP.synch{FeatIND,PolIND,HarmIND}(:,ParamIND);
               [yyy,BFind]=min(abs(log2(chan_BFs/PHs.Info.BF_kHz)));
               
               [ALSR_sps,INDsUsed,oct_range]=ALSRcalc(SynchRates_sps,chan_BFs,CenterFreq);
               
               DFTs.ALSR.ALSRs_sps{FeatIND,PolIND,HarmIND}(BFind,ParamIND)=ALSR_sps;
               DFTs.ALSR.INDsUsed{FeatIND,PolIND,HarmIND}{BFind,ParamIND}=INDsUsed;
               DFTs.ALSR.oct_range{FeatIND,PolIND,HarmIND}(BFind,ParamIND)=oct_range;
            end

         end  % end Param
         

        
         
      end % End FormsatHarms
   end % End Invert Polarity
end % End: FeatIND

%% SAVE PERhists
eval([DFTs_dataname ' =DFTs;'])
clear DFTs

disp(['   *** Saving new file: "' DFTs_filename '" ...'])
eval(['save ''' fullfile(STMPanal_dir,DFTs_filename) '''' ' ' DFTs_dataname])

return;

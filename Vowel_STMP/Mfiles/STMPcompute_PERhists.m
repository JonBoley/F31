function STMPcompute_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% FROM: STMPconvert_STs_Nchans1stim(ExpDate,UnitName,STIMtype)
% M. Heinz Jun 09, 2007
% Computes PERhistograms for a STMP POPULATION.  Assumes STMP format for
% data storage, but this can be either: DATAtype=
%     - 1unitNstim (PHs.Info.STMPshift=0)
%     - Nchans1stim (PHs.Info.STMPshift=1)
%
% Setup to be general for STMP data.
%
% Calls basic function: PERhist
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList

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
RAWdata_dir=fullfile(STMP_dir,'ExpData',ExpName);
eval(['cd ''' RAWdata_dir ''''])
UNITinfo_dir=fullfile(STMP_dir,'ExpData',ExpName,'UNITinfo');   % For general unit info (BF, SR, bad lines, ...)
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify unitINFO and STs exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitINFO_filename=sprintf('unitINFO.%d.%02d.mat',TrackNum,UnitNum);
eval(['load ''' fullfile(UNITinfo_dir,unitINFO_filename) ''''])
STs_filename=sprintf('STs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
STs_dataname=sprintf('STs_%s',DATAtype);
eval(['load ''' fullfile(STMPanal_dir,STs_filename) ''''])
eval(['STs = ' STs_dataname ';']);
eval(['clear ' STs_dataname]);

disp(sprintf('... STMP: Computing ''%s'' PERhists  for:  Experiment: ''%s''; Unit: %d.%02d',DATAtype,ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% See if PHs_filename exists already
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHs_filename=sprintf('PHs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
PHs_dataname=sprintf('PHs_%s',DATAtype);
eval(['ddd=dir(''' fullfile(STMPanal_dir,PHs_filename) ''');'])
if ~isempty(ddd)
   beep
   TEMP = input(sprintf('File: ''%s'' already exists!!\n  ***** Do you want to re-compute PERhists, or leave as is?  [0]: LEAVE AS IS; RERUN: 1;  ',PHs_filename));
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
PHs=STs; 
PHs=rmfield(PHs,'SpikeTrains');

%              Info: [1x1 struct]             % UNCHANGED
%          STIMtype: 'EHrBFi'                 % UNCHANGED
%       ChannelInfo: [1x1 struct]             % UNCHANGED
%     ParameterInfo: [1x1 struct]             % UNCHANGED
%     ConditionInfo: [1x1 struct]             % UNCHANGED
%          StimInfo: [1x1 struct]             % UNCHANGED
%       SpikeTrains: {2x1 cell}               % USED AND DISCARDED
%          PERhists: {2x1 cell}               % COMPUTED 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this unit, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(STs.SpikeTrains{1,1,1},1);
NumP=size(STs.SpikeTrains{1,1,1},2);
NumFEATURES=size(STs.SpikeTrains,1);
NumPOL=size(STs.SpikeTrains,2);
NumHARMS=size(STs.SpikeTrains,3);

PHs.PERhists.PERhists=cell(NumFEATURES,NumPOL,NumHARMS);
PHs.PERhists.PH_orig_NumDrivenSpikes=cell(NumFEATURES,NumPOL,NumHARMS);
PHs.PERhists.PH_binwidth_sec=cell(NumFEATURES,NumPOL,NumHARMS);

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         PHs.PERhists.PERhists{FeatIND,PolIND,HarmIND}=cell(NumCH,NumP);
         for ChanIND=1:NumCH
            for ParamIND=1:NumP
               
               SpikeTrains=STs.SpikeTrains{FeatIND,PolIND,HarmIND}{ChanIND,ParamIND};
               if PHs.Info.STMPshift
                  StimDur_msec=STs.ChannelInfo.StimDur_msec{FeatIND,PolIND,HarmIND}(ChanIND);
               else
                  StimDur_msec=STs.StimInfo.StimDur_msec;
               end
               StimF0_Hz=STs.ChannelInfo.channel_F0_Hz{FeatIND,PolIND,HarmIND}(ChanIND);
               
               [PERhist_sps,PH_binwidth_sec,PH_orig_NumDrivenSpikes] = PERhist(SpikeTrains,StimDur_msec,StimF0_Hz);
               
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               %%% For STMP-shifted data, the PERhist's overall rate needs
               %%% to be scaled to compensate for the scaling of time in STMP.
               if PHs.Info.STMPshift
                  TimeFact = 1/STs.ChannelInfo.FreqFact{FeatIND,PolIND,HarmIND}(ChanIND);
                  PERhist_sps = PERhist_sps * TimeFact;
               end
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

               PHs.PERhists.PERhists{FeatIND,PolIND,HarmIND}{ChanIND,ParamIND}=PERhist_sps;
               % This is the actual number of recorded spikes used to create this PERhist
               % (use this for stats, etc ...) - MAY DIFFER from
               % sum(PERhist) for STMP-shifted conditions due to
               % overall-rate scaling  
               PHs.PERhists.PH_orig_NumDrivenSpikes{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=PH_orig_NumDrivenSpikes;
               PHs.PERhists.PH_binwidth_sec{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND)=PH_binwidth_sec;

            end  % end Param
         end  % end Freq
      end % End FormsatHarms
   end % End Invert Polarity
end % End: FeatIND

%% SAVE PERhists
eval([PHs_dataname ' =PHs;'])
clear PHs

disp(['   *** Saving new file: "' PHs_filename '" ...'])
eval(['save ''' fullfile(STMPanal_dir,PHs_filename) '''' ' ' PHs_dataname])

return;

function UnitAnal_TrBF(ExpDate,UnitName,PLOTyes)
% File: UnitAnal_TrBF.m
% Date: 09Sep2004 (M. Heinz)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Performs Tone_reBF analysis for a given experiment and unit.  Loads 'DataList' file, and unit file.
% UnitPlot_TrBF.m plots results from this analysis
%

global NOHR_dir NOHR_ExpList

if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='070804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
%%% HARD CODE FOR NOW
UnitName='3.07'
   while ~ischar(UnitName)
      UnitName=input('Enter Unit Name (e.g., ''3.07''): ');
   end
end
if ~exist('PLOTyes','var')
   PLOTyes=0;
end

% Find the full Experiment Name 
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
   break
end

% Parse out the Track and Unit Number 
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(sprintf('Processing  (TreBF) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

% Load DataList file with all picture numbers
if length(dir(FileName))
   disp(['   *** Loading file: "' FileName '"'])
   eval(['load ' FileName])
else
   error(sprintf('%s: FILE NOT FOUND!!!!!',FileName))
end

% Veriy Unit Number is legitimate
[NumTracks,NumUnits]=size(DataList.Units);
if (TrackNum>NumTracks)|(UnitNum>NumUnits)
   error(sprintf('Unit Number: %d.%02d does not exist!!!!!',TrackNum,UnitNum))
else
   if isempty(DataList.Units{TrackNum,UnitNum})
      error(sprintf('Unit Number: %d.%02d does not exist!!!!!',TrackNum,UnitNum))
   elseif ~isfield(DataList.Units{TrackNum,UnitNum},'Tone_reBF')
      error(sprintf('"Tone_reBF" data does not exist for Unit Number: %d.%02d',TrackNum,UnitNum))
   end
end

%%%%%%%% Create unit structure for this unit
% LOAD UNITSdata file, if it exists, otherwise take from DataList
%%%%%%%%%%%%%%%
% ???LATER??? How to make sure DataList is not updated, re unitdata????%
% ? Can we write code to synchronize the two????
%
%%%%%%%%%%%%%%%
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
%%%%%%%%%%%%%%%%
%% LATER: load unit data if it exists,
if ~isempty(dir(fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum))))
   eval(['load ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''''])
else
   unit=DataList.Units{TrackNum,UnitNum};
   unit.IgnorePicNums=[];
end

% Update IgnorePics, from stored unit data and master list
UnitPicList=findPics('*',[TrackNum,UnitNum]);
unit.IgnorePicNums=union(intersect(UnitPicList,DataList.IgnorePicNums),unit.IgnorePicNums);


%%%%%%%%%%% Calculate rate, synch, and phase for all conditions
unit.Tone_reBF.rate=NaN+zeros(size(unit.Tone_reBF.picNums));
unit.Tone_reBF.synch=NaN+zeros(size(unit.Tone_reBF.picNums));
unit.Tone_reBF.phase=NaN+zeros(size(unit.Tone_reBF.picNums));
unit.Tone_reBF.RaySig=NaN+zeros(size(unit.Tone_reBF.picNums));

for freqIND=1:length(unit.Tone_reBF.freqs_kHz)
   for levelIND=1:length(unit.Tone_reBF.levels_dBSPL)
      if ~isempty(unit.Tone_reBF.picNums{levelIND,freqIND})
         
         % Gather spikes from relevant pictures
         PIC=concatPICS_NOHR(unit.Tone_reBF.picNums{levelIND,freqIND},unit.Tone_reBF.excludeLines{levelIND,freqIND});
         % Calc Synch and Phase from PERIOD histogram
         PIC=calcSynchRate_PERhist(PIC);

         %          unit.Tone_reBF.rate(levelIND,freqIND)=PIC.PERhist.NumDrivenSpikes/PIC.x.Stimuli.fully_presented_lines/ ...
         %             diff(PIC.PERhist.params.PERhist_window_sec);  
         unit.Tone_reBF.rate(levelIND,freqIND)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);
         unit.Tone_reBF.synch(levelIND,freqIND)=PIC.SynchRate_PERhist.FeatureSynchs;
         unit.Tone_reBF.phase(levelIND,freqIND)=PIC.SynchRate_PERhist.FeaturePhases;
         unit.Tone_reBF.RaySig(levelIND,freqIND)=PIC.SynchRate_PERhist.FeatureRaySig;
         
      end
   end
end

%%%%%%%%%%%% Save unit data at end
eval(['save ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''' unit'])


%%%%%%%%%%%% Plot Tone_reBF data
if PLOTyes
   UnitPlot_TrBF(ExpDate,UnitName)
end
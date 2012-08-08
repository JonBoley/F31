function UnitAnal_EHvNrBF(ExpDate,UnitName,PLOTyes)
% File: UnitAnal_EHvNrBF.m
% Modified Date: 27Dec2005 - from UnitAnal_EHrBF.m
% Modified Date: 07Jan2005 - consolidated with UnitAnal_EHrBFi - NOW handles both
% Date: 17Sep2004 (M. Heinz) (Modified from UnitAnal_TrBF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Performs EHvN_reBF analysis for a given experiment and unit.  Loads 'DataList' file, and unit file.
%

global NOHR_dir NOHR_ExpList FeaturesText
global SavedPICS SavedPICnums SavedPICSuse
SavedPICSuse=1; SavedPICS=[]; SavedPICnums=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ExpDate','var')
   ExpDate=0;
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
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
	error('STOPPED')
end

% Parse out the Track and Unit Number 
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData',ExpName);
unitdata_dir=fullfile(data_dir,'UNITSdata');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

cd(data_dir)
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(sprintf('Processing (EHvNreBF) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

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
   elseif ~isfield(DataList.Units{TrackNum,UnitNum},'EHvN_reBF')
       if ~isfield(DataList.Units{TrackNum,UnitNum},'EHvLTASS_reBF')
           error(sprintf('"EHvN_reBF" data does not exist for Unit Number: %d.%02d',TrackNum,UnitNum))
       end
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
if ~isempty(dir(fullfile(unitdata_dir,UnitFileName)))
   eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
else
   unit=DataList.Units{TrackNum,UnitNum};
   unit.IgnorePicNums=[];
end

% Update IgnorePics, from stored unit data and master list
if ~isfield(unit,'IgnorePicNums')
	unit.IgnorePicNums=[];
end

UnitPicList=findPics('*',[TrackNum,UnitNum]);
unit.IgnorePicNums=union(intersect(UnitPicList,DataList.IgnorePicNums),unit.IgnorePicNums);


%%%%%%%%%%%%%%%%% NEED TO FIND ALL FEATURES!!!!!!!!!!!!!!
try
    UnitFeats=fieldnames(unit.EHvN_reBF);
catch
    UnitFeats=fieldnames(unit.EHvLTASS_reBF);
end
UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
clear FeatINDs
for FeatIND=1:length(UnitFeats)
   FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
end
for FeatIND=FeatINDs  % Step through each Feature we have data for
   
	disp(sprintf('       ... processing Feature: %s (EHvN_reBF)',FeaturesText{FeatIND}))

   for HarmonicsIND=1:2
      for PolarityIND=1:2
         
         try
             eval(['yTEMP=unit.EHvN_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         catch
             eval(['yTEMP=unit.EHvLTASS_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         end
         if ~isempty(yTEMP)
            %%%%%%%%%%% Calculate rate, synch, and phase for all conditions
            yTEMP.rate=NaN+zeros(size(yTEMP.picNums));
            yTEMP.FeatureFreqs_Hz=cell(1,size(yTEMP.picNums,2));
            yTEMP.FeatureLevels_dB=[];
            yTEMP.synch=cell(size(yTEMP.picNums));
            yTEMP.phase=cell(size(yTEMP.picNums));
            yTEMP.RaySig=cell(size(yTEMP.picNums));
            
            for freqIND=1:length(yTEMP.freqs_kHz)
               for attenIND=1:length(yTEMP.Nattens_dB)
                  if ~isempty(yTEMP.picNums{attenIND,freqIND})
                     
							%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
							%%%% Update excludeLines from DataList before
							%%%% proceeding.  This updated list will be stored in
							%%%% unit and used throughout.
							for picIND = 1:length(yTEMP.picNums{attenIND,freqIND})
								masterBADlines = DataList.picBADlines{yTEMP.picNums{attenIND,freqIND}(picIND)};
								oldBADlines = yTEMP.excludeLines{attenIND,freqIND}{picIND};
								newBADlines = union(oldBADlines,masterBADlines);
								yTEMP.excludeLines{attenIND,freqIND}{picIND} = newBADlines;
							end

							% Gather spikes from relevant pictures
                     PIC=concatPICS_NOHR(yTEMP.picNums{attenIND,freqIND},yTEMP.excludeLines{attenIND,freqIND});
                     % Calc Synch and Phase from PERIOD histogram
                     if ~isfield(PIC,'FundamentalFreq_Hz')
                         F0index=4+strfind(PIC.x.Stimuli.list{1},'F0at'); 
                         PIC.FundamentalFreq_Hz = str2num(PIC.x.Stimuli.list{1}(F0index:F0index+2));
                     end
                     PIC=calcSynchRate_PERhist(PIC);
                     
                     %                      yTEMP.rate(attenIND,freqIND)=PIC.PERhist.NumDrivenSpikes/PIC.x.Stimuli.fully_presented_lines/ ...
                     %                         diff(PIC.PERhist.params.PERhist_window_sec);
                     yTEMP.rate(attenIND,freqIND)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);  % 1st DFT coefficient is Mean Rate
                     
                     yTEMP.synch{attenIND,freqIND}=PIC.SynchRate_PERhist.FeatureSynchs;
                     yTEMP.phase{attenIND,freqIND}=PIC.SynchRate_PERhist.FeaturePhases;
                     yTEMP.RaySig{attenIND,freqIND}=PIC.SynchRate_PERhist.FeatureRaySig;
                     
                     yTEMP.FeatureFreqs_Hz{freqIND}=PIC.FeatureFreqs_Hz;
                     yTEMP.FeatureLevels_dB=[NaN PIC.FeatureLevels_dB];
                  end
               end
            end
            
            try
                eval(['unit.EHvN_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND}=yTEMP;'])
            catch
                eval(['unit.EHvLTASS_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND}=yTEMP;'])
            end
            
         end
      end
   end
end

%%%%%%%%%%%% Save unit data at end
eval(['save ''' fullfile(unitdata_dir,UnitFileName) ''' unit'])
disp(['   *** Saving file: "' UnitFileName '"'])

% Turn off saved PICS feature
SavedPICS=[]; SavedPICnums=[];
SavedPICSuse=0;

return;

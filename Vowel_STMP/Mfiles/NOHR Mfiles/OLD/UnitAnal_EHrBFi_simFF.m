function UnitAnal_EHrBFi_simFF(ExpDate,UnitName,PLOTyes)
% File: UnitAnal_EHrBFi_simFF.m
% Modified Date: 122804 (M Heinz)
% Modified From: UnitAnal_EHrBF_simFF.m to use Interleaved Data
%
% File: UnitAnal_EHrBF_simFF.m
% Date: 28Sep2004 (M. Heinz) (Modified from UnitAnal_EHrBF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Performs EH_reBF_simFF analysis for a given experiment and unit.  Loads 'DataList' file, and unit file.
% UnitPlot_EHrBF_simFF.m plots results from this analysis
%

global NOHR_dir NOHR_ExpList 
global FeaturesText FormsAtHarmonicsText InvertPolarityText

%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='111804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
%%% HARD CODE FOR NOW
UnitName='1.28'
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

disp(sprintf('Processing (EHreBFi_simFF) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

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
   elseif ~isfield(DataList.Units{TrackNum,UnitNum},'EH_reBF')  %updateBFi
      error(sprintf('"EH_reBF" data does not exist for Unit Number: %d.%02d',TrackNum,UnitNum))
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
   
   % If EH_reBFi analysis is not completed, run here
   UnitFeats=fieldnames(unit.EH_reBF);  %updateBFi
   UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
   clear processedBOOL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         eval(['emptyBOOL=~isempty(unit.EH_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND});'])   %updateBFi
         if emptyBOOL
            eval(['processedBOOL=isfield(unit.EH_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND},''synch'');'])           %updateBFi
            if ~processedBOOL
               UnitAnal_EHrBFi(ExpDate,UnitName,0);   %updateBFi
               eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
            end
            break
         end
      end
      if exist('processedBOOL','var')
         break
      end
   end
   
else
   unit=DataList.Units{TrackNum,UnitNum};
   unit.IgnorePicNums=[];
end

% Update IgnorePics, from stored unit data and master list
UnitPicList=findPics('*',[TrackNum,UnitNum]);
unit.IgnorePicNums=union(intersect(UnitPicList,DataList.IgnorePicNums),unit.IgnorePicNums);


%%%%%%%%%%%%%%%%% NEED TO FIND ALL FEATURES!!!!!!!!!!!!!!
UnitFeats=fieldnames(unit.EH_reBF);   %updateBFi
UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
clear FeatINDs
for FeatIND=1:length(UnitFeats)
   FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
end
for FeatIND=FeatINDs  % Step through each Feature we have data for
   
   eval(['unit.EH_reBF_simFF.' FeaturesText{FeatIND} '=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));'])     %updateBFi
   
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         
         eval(['yTEMP0=unit.EH_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])   %updateBFi
         if ~isempty(yTEMP0)
            %%%%%%%%%%% Calculate rate, synch, and phase for all conditions

            yTEMP.BFs_kHz=NaN+zeros(size(yTEMP0.freqs_kHz));
            yTEMP.levels_dBSPL=yTEMP0.levels_dBSPL;
            yTEMP.picNums=yTEMP0.picNums;
            yTEMP.excludeLines=yTEMP0.excludeLines;
            
            yTEMP.rate=NaN+zeros(size(yTEMP.picNums));
            yTEMP.FeatureFreqs_Hz=cell(1,size(yTEMP.picNums,2));
            yTEMP.FeatureLevels_dB=[];
            yTEMP.synch=cell(size(yTEMP.picNums));
            yTEMP.phase=cell(size(yTEMP.picNums));
            yTEMP.RaySig=cell(size(yTEMP.picNums));
            
            for freqIND=1:length(yTEMP0.freqs_kHz)
               for levelIND=1:length(yTEMP.levels_dBSPL)
                  if ~isempty(yTEMP.picNums{levelIND,freqIND})
                     
                     % Gather spikes from relevant pictures
                     PIC=concatPICS_NOHR(yTEMP.picNums{levelIND,freqIND},yTEMP.excludeLines{levelIND,freqIND});
                     
                     % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                     PIC=simFF_PICshift(PIC);
                     yTEMP.BFs_kHz(freqIND)=PIC.BF_Hz/1000;
                     
                     % Calc Synch and Phase from PERIOD histogram
                     PIC=calcSynchRate_PERhist(PIC);
                     
                     %                      yTEMP.rate(levelIND,freqIND)=PIC.PERhist.NumDrivenSpikes/PIC.x.Stimuli.fully_presented_lines/ ...
                     %                         diff(PIC.PERhist.params.PERhist_window_sec);
                     yTEMP.rate(levelIND,freqIND)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);  % 1st DFT coefficient is Mean Rate
                     
                     yTEMP.synch{levelIND,freqIND}=PIC.SynchRate_PERhist.FeatureSynchs;
                     yTEMP.phase{levelIND,freqIND}=PIC.SynchRate_PERhist.FeaturePhases;
                     yTEMP.RaySig{levelIND,freqIND}=PIC.SynchRate_PERhist.FeatureRaySig;
                     
                     yTEMP.FeatureFreqs_Hz{freqIND}=PIC.FeatureFreqs_Hz;
                     yTEMP.FeatureLevels_dB=[NaN PIC.FeatureLevels_dB];
                  end
               end
            end
            
            eval(['unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND}=yTEMP;'])    %updateBFi
            
         end
      end
   end
end

%%%%%%%%%%%% Save unit data at end
eval(['save ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''' unit'])


%%%%%%%%%%%% Plot EH_reBF_simFF data
if PLOTyes
   UnitPlot_EHrBFi_simFF(ExpDate,UnitName)
end
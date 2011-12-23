function UnitAnal_EHvNrBF_simFF(ExpDate,UnitName,PLOTyes)
% File: UnitAnal_EHvNrBF_simFF.m
% Modified Date: 07Jan2005: incorporated EHrBFi (interleaving) here
% Date: 28Sep2004 (M. Heinz) (Modified from UnitAnal_EHrBF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Performs EHINvN_reBF_simFF analysis for a given experiment and unit.  Loads 'DataList' file, and unit file.
%

global NOHR_dir NOHR_ExpList 
global FeaturesText FormsAtHarmonicsText InvertPolarityText
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

cd(data_dir)
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(sprintf('Processing (EHvNreBF_simFF) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

% Load DataList file with all picture numbers
if length(dir(FileName))
   disp(['   *** Loading file: "' FileName '"'])
   eval(['load ' FileName])
else
   error(sprintf('%s: FILE NOT FOUND!!!!!',FileName))
end

% Verify Unit Number is legitimate
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
   
   % If EHvN_reBF analysis is not completed, run here
   try
       UnitFeats=fieldnames(unit.EHvLTASS_reBF);
   catch
       UnitFeats=fieldnames(unit.EHvN_reBF);
   end
   UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
	% CHECK ALL POSSIBLE Harm/Pol places to look
	clear processedBOOL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
          try
              eval(['emptyBOOL=~isempty(unit.EHvLTASS_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND});'])
          catch
              eval(['emptyBOOL=~isempty(unit.EHvN_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND});'])
          end
         if emptyBOOL
             try
                 eval(['processedBOOL=isfield(unit.EHvLTASS_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND},''synch'');'])
             catch
                 eval(['processedBOOL=isfield(unit.EHvN_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND},''synch'');'])
             end
            if ~processedBOOL
               UnitAnal_EHvNrBF(ExpDate,UnitName,0);   
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
	UnitAnal_EHvNrBF(ExpDate,UnitName,0);
	eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
end

% Update IgnorePics, from stored unit data and master list
UnitPicList=findPics('*',[TrackNum,UnitNum]);
unit.IgnorePicNums=union(intersect(UnitPicList,DataList.IgnorePicNums),unit.IgnorePicNums);


%%%%%%%%%%%%%%%%% NEED TO FIND ALL FEATURES!!!!!!!!!!!!!!
try
    UnitFeats=fieldnames(unit.EHvLTASS_reBF);
catch
    UnitFeats=fieldnames(unit.EHvN_reBF);
end
UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
clear FeatINDs
for FeatIND=1:length(UnitFeats)
   FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
end
for FeatIND=FeatINDs  % Step through each Feature we have data for
   try
      eval(['unit.EHvLTASS_reBF_simFF.' FeaturesText{FeatIND} '=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));']) 
   catch   
       eval(['unit.EHvN_reBF_simFF.' FeaturesText{FeatIND} '=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));'])   
   end
   disp(sprintf('       ... processing Feature: %s (EHvN_reBF_simFF)',FeaturesText{FeatIND}))
	
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         try
            eval(['yTEMP0=unit.EHvLTASS_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         catch 
             eval(['yTEMP0=unit.EHvN_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         end
         if ~isempty(yTEMP0)
            %%%%%%%%%%% Calculate rate, synch, and phase for all conditions

            yTEMP.BFs_kHz=NaN+zeros(size(yTEMP0.freqs_kHz));
            yTEMP.levels_dBSPL=yTEMP0.levels_dBSPL;
            yTEMP.SORTattens_indices=yTEMP0.SORTattens_indices;
            yTEMP.EqualSPL_index=yTEMP0.EqualSPL_index;
            yTEMP.EqualSL_index=yTEMP0.EqualSL_index;
            yTEMP.SNR_EqualSL=yTEMP0.SNR_EqualSL;
            yTEMP.Nattens_dB=yTEMP0.Nattens_dB;
            yTEMP.picNums=yTEMP0.picNums;
            yTEMP.excludeLines=yTEMP0.excludeLines;
            
            yTEMP.rate=NaN+zeros(size(yTEMP.picNums));
            yTEMP.FeatureFreqs_Hz=cell(1,size(yTEMP.picNums,2));
            yTEMP.FeatureLevels_dB=[];
            yTEMP.synch=cell(size(yTEMP.picNums));
            yTEMP.phase=cell(size(yTEMP.picNums));
            yTEMP.RaySig=cell(size(yTEMP.picNums));
            yTEMP.TimeFact=cell(size(yTEMP.picNums));
            
            for freqIND=1:length(yTEMP0.freqs_kHz)
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
                     
                     % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                     PIC=simFF_PICshift(PIC);
                     yTEMP.BFs_kHz(freqIND)=PIC.BF_Hz/1000;
                     
                     % Calc Synch and Phase from PERIOD histogram
                     PIC=calcSynchRate_PERhist(PIC);
                     
                     %                      yTEMP.rate(attenIND,freqIND)=PIC.PERhist.NumDrivenSpikes/PIC.x.Stimuli.fully_presented_lines/ ...
                     %                         diff(PIC.PERhist.params.PERhist_window_sec);
                     yTEMP.rate(attenIND,freqIND)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);  % 1st DFT coefficient is Mean Rate
                     
                     yTEMP.synch{attenIND,freqIND}=PIC.SynchRate_PERhist.FeatureSynchs;
                     yTEMP.phase{attenIND,freqIND}=PIC.SynchRate_PERhist.FeaturePhases;
                     yTEMP.RaySig{attenIND,freqIND}=PIC.SynchRate_PERhist.FeatureRaySig;
                     yTEMP.TimeFact{attenIND,freqIND}=PIC.simFF_PICshift.TimeFact;
                     
                     yTEMP.FeatureFreqs_Hz{freqIND}=PIC.FeatureFreqs_Hz;
                     yTEMP.FeatureLevels_dB=[NaN PIC.FeatureLevels_dB];
                  end
               end
            end
            
            try
               eval(['unit.EHvLTASS_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND}=yTEMP;']) 
            catch
                eval(['unit.EHvN_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND}=yTEMP;'])
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

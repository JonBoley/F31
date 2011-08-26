function updateNOHR_DataList_excludeLines(ExpDate,CHOOSEone)
% File: updateNOHR_DataList_excludeLines.m
% Date: 27Apr2006 (M. Heinz)
% For: R03 Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Modified From: makeNOHRdataList.m 
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
% Updates excludeLines fields in DataList, after storeBADlines has been run
% (which is run after initial makeNOHRdataList).
%
% So, order is: makeNOHRdataList.m, storeBADlines.m,
% updateNOHR_DataList_excludeLines.m. 
%
% OBVIOUSLY, this can be consolidated later to clean this up.  
% THIS is setup for 041905 experiment, but could be extended later
%
% 042706 - STOPPED SETTING THIS UP B/C TOO COMPLICATED AND NOT ENOUGH TIME
%        - just went straight to computations to include both sources of
%        exclude lines!!!!!
%        - comme back to this when more time
%

global NOHR_dir NOHR_ExpList NOHR_IMPExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

% Globals for functions
global TUlist UNITind DataList OLD_DataList

if ~exist('ExpDate','var')
   ExpDate=0;
   ExpDate='041805'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end

% Allows user to choose one unit to process, rather than ALL 
% [0: process ALL, 1: choose 1 unit to process]
if ~exist('CHOOSEone','var')
   CHOOSEone=0;
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
   error('STOPPED');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(['Processing Experiment: ''' ExpName ''' ... updating excludeLines ... '])

%%%%%% Make sure file exist already and save old version to be careful
if length(dir(FileName))
	disp(['   *** Loading old file: "' FileName '"; ... only excludeLines fields will be updated'])
   eval(['load ' FileName])
	eval(['save OLD_' FileName ' DataList'])
	disp(['   *** Saving copy of old file: "OLD_' FileName '";'])
else
	error('DataList file does not exist')
end


%%%%%%%%%%%%%%%% Get Track/Unit List for Experiment
TUlist=getTrackUnitList;
TOTALunits=size(TUlist,1);
if isempty(TUlist)
   error('******NO UNITS: File not setup yet to handle this!!!!')
end

%%%%%%%%%%% ALLOW only ONE UNIT to be done
% e.g., for DEBUGGING
if CHOOSEone
	Track=input('Enter Track (e.g., 1): ');
	Unit=input('Enter Unit (e.g., 5): ');
	ONEunit=find((TUlist(:,1)==Track)&(TUlist(:,2)==Unit));
	TUlist_ORIGINAL=TUlist;
	TUlist=TUlist(ONEunit,:);
end

for UNITind=1:size(TUlist,1)   
   disp(sprintf('   Unit: %d.%02d',TUlist(UNITind,:)))
    
   %%%%%%%%%%%%% TONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
   %%%%%%%%%%%%%
   %%% TONErlv (at BF) Info
   %%%%%%%%%%%%%
%    DataList=TONErlv_update(DataList);   
   
   %%%%%%%%%%%%%
   %%% Tone_reBF Info
   %%%%%%%%%%%%%
%    DataList=TrBF_update(DataList);   
   
   %%%%%%%%%%%%%
   %%% Tone_reBFi Info
   %%%%%%%%%%%%%
%    DataList=TrBFi_update(DataList);   
   
   %%%%%%%%%%%%%
   %%% Tone_reFF Info
   %%%%%%%%%%%%%
%    DataList=TrFF_update(DataList);
   
   %%%%%%%%%%%%% VOWELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   %%%%%%%%%%%%%
   %%% EHrlv Info
   %%%%%%%%%%%%%
%    DataList=EHrlv_update(DataList);   

   %%%%%%%%%%%%%
   %%% EHINrlv Info
   %%%%%%%%%%%%%
%    DataList=EHINrlv_update(DataList);   
   

   
   
   %%%%%%%%%%%%%
   %%% EHvNrBFi Info
   %%%%%%%%%%%%%
   DataList=EHvNrBFi_update(DataList);   
  

   %%%%%%%%%%%%%
   %%% EHvLTASSrBFi Info
   %%%%%%%%%%%%%
   DataList=EHvLTASSrBFi_update(DataList); 



   
   %%%%%%%%%%%%%
   %%% EH_reBF Info
   %%%%%%%%%%%%%
%    DataList=EHrBF_update(DataList);
   
   %%%%%%%%%%%%%
   %%% EH_reBFi Info
   %%%%%%%%%%%%%
%    DataList=EHrBFi_update(DataList);
   
   %%%%%%%%%%%%%
   %%% EH_reFF Info
   %%%%%%%%%%%%%
%    DataList=EHrFF_update(DataList);
   
   
   %%%%%%%%%%%%%
   %%% SACrlvQ Info
   %%%%%%%%%%%%%
%    DataList=SACrlvQ_update(DataList);   
   
   %%%%%%%%%%%%%
   %%% SACrlv Info
   %%%%%%%%%%%%%
%    DataList=SACrlv_update(DataList);   
   
end % end step through UNITS


disp(['   *** Saving new file: "' FileName '" ...'])
% disp(['   *** Saving new file: "' FileName '" ...'])
% eval(['save ' FileName ' DataList'])

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHrlv_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText

picList=findPics('EHrlv',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('EHrlv data Mismatch for Picture: %d',picList(ind)));
      end
      
      %%% CALC frequency from Stimuli.Condition (b/c this will be the same for all Features, even if not exactly the right frequency)
      OffsetSign=sign(find(strcmp(x.Stimuli.Condition.Offset_Direction,{'below ','above '}))-1.5);
      if isstr(x.Stimuli.Condition.FreqOffset_octs)
         Offset_octs=str2num(x.Stimuli.Condition.FreqOffset_octs);
      else
         Offset_octs=x.Stimuli.Condition.FreqOffset_octs;
      end
      Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
      %          Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
      Temp.FeatureIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.Feature),FeaturesText));
      Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
      Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
   end
   
   SORTfreqs=unique(Temp.freqs_Hz);
   %      SORTlevels=unique(Temp.levels_dBSPL);
   NumF=length(SORTfreqs);
   NumL=1;  % Keep this just to keep format the same as for other EH stimuli
   SORTfeatures=unique(Temp.FeatureIND);
   
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            TempINDs=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)&(Temp.InvertPolarityIND==InvertPolarityIND));
            
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=SORTfreqs/1000;
               %                   yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
               
               %%%%%%%%%%% STORE picNUMS
               for FreqIND=1:NumF
                  LevelIND=1;
                  TempINDs2=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)& ...
                     (Temp.InvertPolarityIND==InvertPolarityIND)&(Temp.freqs_Hz==SORTfreqs(FreqIND)));
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
               end                  
               
            end
         end % End Invert Polarity
      end % End FormsatHarms
      eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EHrlv.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
   end % End: FeatIND
end

return;  % EHrlv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHINrlv_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText

picList=findPics('EHINrlv',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('EHINrlv data Mismatch for Picture: %d',picList(ind)));
      end
      
      %%% CALC frequency from Stimuli.Condition (b/c this will be the same for all Features, even if not exactly the right frequency)
      OffsetSign=sign(find(strcmp(x.Stimuli.Condition.Offset_Direction,{'below ','above '}))-1.5);
      if isstr(x.Stimuli.Condition.FreqOffset_octs)
         Offset_octs=str2num(x.Stimuli.Condition.FreqOffset_octs);
      else
         Offset_octs=x.Stimuli.Condition.FreqOffset_octs;
      end
      Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
      %          Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
      Temp.FeatureIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.Feature),FeaturesText));
      Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
      Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
   end
   
   SORTfreqs=unique(Temp.freqs_Hz);
   SORTlevels=unique(Temp.levels_dBSPL);
   NumF=length(SORTfreqs);
   NumL=length(SORTlevels);
   SORTfeatures=unique(Temp.FeatureIND);
   
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            TempINDs=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)&(Temp.InvertPolarityIND==InvertPolarityIND));
            
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=SORTfreqs/1000;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
               
               %%%%%%%%%%% STORE picNUMS
               for FreqIND=1:NumF
                  for LevelIND=1:NumL;
                     TempINDs2=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)& ...
                        (Temp.InvertPolarityIND==InvertPolarityIND)&(Temp.freqs_Hz==SORTfreqs(FreqIND)));
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
                  end   
               end
               
            end
         end % End Invert Polarity
      end % End FormsatHarms
      eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EHINrlv.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
   end % End: FeatIND
end

return;  % EHINrlv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHrBF_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText

picList=findPics('EHrBF',TUlist(UNITind,:));
if ~isempty(picList)
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBF.interleaved='no';
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('EHrBF data Mismatch for Picture: %d',picList(ind)));
      end
      
      %%% CALC frequency from Stimuli.Condition (b/c this will be the same for all Features, even if not exactly the right frequency)
      OffsetSign=sign(find(strcmp(x.Stimuli.Condition.Offset_Direction,{'below ','above '}))-1.5);
      if isstr(x.Stimuli.Condition.FreqOffset_octs)
         Offset_octs=str2num(x.Stimuli.Condition.FreqOffset_octs);
      else
         Offset_octs=x.Stimuli.Condition.FreqOffset_octs;
      end
      Temp.octshifts(ind)=OffsetSign*Offset_octs;
      Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
      Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
      Temp.FeatureIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.Feature),FeaturesText));
      Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
      Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
   end
   
   SORTocts=unique(Temp.octshifts);
   SORTfreqs=unique(Temp.freqs_Hz);
   SORTlevels=unique(Temp.levels_dBSPL);
   NumF=length(SORTfreqs);
   NumL=length(SORTlevels);
   SORTfeatures=unique(Temp.FeatureIND);
   
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            TempINDs=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)&(Temp.InvertPolarityIND==InvertPolarityIND));
            
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.octshifts=SORTocts;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=SORTfreqs/1000;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
               
               %%%%%%%%%%% STORE picNUMS
               for FreqIND=1:NumF
                  for LevelIND=1:NumL
                     TempINDs2=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)& ...
                        (Temp.InvertPolarityIND==InvertPolarityIND)&(Temp.freqs_Hz==SORTfreqs(FreqIND))& ...
                        (Temp.levels_dBSPL==SORTlevels(LevelIND)));
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
                  end
               end                  
               
            end
         end % End Invert Polarity
      end % End FormsatHarms
      eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
   end % End: FeatIND
end 

return; %EHrBF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHrBFi_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText
% Modified 11Apr2005 M. Heinz to include TEH data files   

picList=findPics('EHrBFi',TUlist(UNITind,:));
if ~isempty(picList)
   if isfield(DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)},'EH_reBF')
      error('This unit has both EH_reBF and EH_reBFi data, makeNOHRdataList IS NOT SET UP TO COMBINE THESE!!!! (NEED TO DO IT NOW)')
   end
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBF.interleaved='yes';
   
   % FIRST, find all octshifts, Features, freqs 
   clear Temp xx
   Temp.octshifts=cell(1,length(picList));
   Temp.FeatureINDs=cell(1,length(picList));
   Temp.freqs_Hz=cell(1,length(picList));
   Temp.FeaturesList=cell(1,length(picList));
   Temp.levels_dBSPL=cell(1,length(picList));
   Temp.FirstBADline=Inf+ones(size(picList));
   for ind=1:length(picList)
      disp(sprintf('      ... Loading picture: %d',picList(ind)))
      x=loadPic(picList(ind));
      xx{ind}=x;
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('EHrBFi data Mismatch for Picture: %d',picList(ind)));
      end
      
      %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
      if ~isempty(x.Stimuli.bad_lines)
         Temp.FirstBADline(ind)=x.Stimuli.bad_lines(1);
         disp(sprintf('EHrBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(ind),Temp.FirstBADline(ind)));
      end
      
      Temp.octshifts{ind}=x.Stimuli.Used.OctShifts_List;
      Temp.FeatureINDs{ind}=find(strcmp(deblank(x.Stimuli.Condition.Features{1}),FeaturesText));
      for i=2:length(x.Stimuli.Condition.Features)
         Temp.FeatureINDs{ind}=[Temp.FeatureINDs{ind} find(strcmp(deblank(x.Stimuli.Condition.Features{i}),FeaturesText))];
      end
      Temp.freqs_Hz{ind}=x.Stimuli.Used.FeatureTarget_Hz_List;
      Temp.FeaturesList{ind}=x.Stimuli.Used.Features_List;
      %          Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
      if isfield(x.Stimuli.Used,'Levels_dBSPL_List')
         Temp.levels_dBSPL{ind}=x.Stimuli.Used.Levels_dBSPL_List;
      else
         Temp.levels_dBSPL{ind}=x.Stimuli.Condition.Level_dBSPL;
      end
      Temp.picNums(ind)=picList(ind);
      Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
      Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
   end
   
   TempFeatureINDs=[];
   Tempoctshifts=[];
   Templevels_dBSPL=[];
   for i=1:length(picList)
      Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
      TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
      Templevels_dBSPL=[Templevels_dBSPL Temp.levels_dBSPL{i}];
   end
   SORTocts=unique(Tempoctshifts);
   SORTfeatures=unique(TempFeatureINDs);
   SORTlevels=unique(Templevels_dBSPL);
   NumF=length(SORTocts);
   NumL=length(SORTlevels);
   
   TONEfeatureIND=find(strcmp(FeaturesText,'TN'));
   %% Do all NON-tone features first
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            
            % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
            clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
            for i=1:length(picList)
               TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
            end
            TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
            TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
            TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity
            
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.octshifts=SORTocts;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=[];
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
               
               %%%%%%%%%%% STORE picNUMS for each Freq and Level
               for FreqIND=1:NumF
                  
                  % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
                  CONDindFULL=[];
                  for PICind=TempINDs
                     CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                  end
                  CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
                  if length(CONDfreq_kHz)~=1
                     error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
                  end
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz(FreqIND)=CONDfreq_kHz;
                  
                  for LevelIND=1:NumL
                     % Verify unique conditions
                     CONDindFULL=[];
                     for PICind=TempINDs
                        CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.levels_dBSPL{PICind}==SORTlevels(LevelIND))& ...
                              (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                     end
                     CONDind=unique(CONDindFULL);
                     if length(CONDind)>1
                        error(sprintf('Non-unique CONDind in TEHeBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND))
                     end
                     
                     % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
                     clear TempOctINDs TempLevINDs
                     for i=1:length(picList)
                        TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
                        TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL{i},SORTlevels(LevelIND)))>0;
                     end
                     TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Freq,Level
                     
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
                     % Store excludeLines for each picture, based on CONDind and TOTALconds
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
                     ii=0;
                     for i=TempINDs2
                        ii=ii+1;
                        TOTALconds=length(xx{i}.Stimuli.Used.FeatureTarget_Hz_List);
                        GOODlines=CONDind:TOTALconds:xx{i}.Stimuli.fully_presented_lines;
                        %% Take out any lines beyond badlines
                        GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
                        IGNORElines=setdiff(1:xx{i}.Stimuli.fully_presented_lines,GOODlines);
                        yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}{ii}=IGNORElines;
                     end
                  end
               end
            end
         end % End Invert Polarity
      end % End FormsatHarms
      % Need to save 'TN' to a different place
      if FeatIND==TONEfeatureIND  %Tone_reBFi
         if isempty(yTEMP{1,1})
            error(sprintf('Unit %d.%02d has TEH_reBFi data for TONE (TN) not in [Harm,Pol]=[1,1], need to fix code!!!',TUlist(UNITind,1),TUlist(UNITind,2)))
         else
            DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF=yTEMP{1,1};
            DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.interleaved='yes';
         end
      else % EHreBF
         eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
      end
   end % End: FeatIND
end 

return; %EHrBFi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHvNrBFi_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CHECK THIS LATER
disp('CHECK THIS LATER')
beep
beep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

picList=DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.picNums;
if ~isempty(picList)
	for ind=1:length(picList)
		picNUM=picList(ind);
		oldLIST = DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind};
		newLIST = union(oldLIST,DataList.picBADlines{picNUM});
		DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind}=newLIST;
	end
end



      
  
   NumF=length(SORTocts);
   NumL=length(SORTlevels);
   NumA=length(SORTattens);

   if NumL~=1
      error('There is more than one vowel level in this EHvNreBFI data, makeNOHRdataList NOT SETUP FOR THIS!')
   end
   
   %    TONEfeatureIND=find(strcmp(FeaturesText,'TN'));
   %% Do all NON-tone features first
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            
            % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
            clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
            for i=1:length(picList)
               TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
            end
            TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
            TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
            TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity
            
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.octshifts=SORTocts;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=[];
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.Nattens_dB=SORTattens;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumA,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumA,NumF);
               
               %%%%%%%%%%% STORE picNUMS for each Freq and Level
               for FreqIND=1:NumF
                  
                  % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
                  CONDindFULL=[];
                  for PICind=TempINDs
                     CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                  end
                  CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
                  if length(CONDfreq_kHz)~=1
                     error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
                  end
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz(FreqIND)=CONDfreq_kHz;
                  
                  for AttenIND=1:NumA
                     % Verify unique conditions
                     CONDindFULL=[];
                     for PICind=TempINDs
                        CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.Nattens_dB{PICind}==SORTattens(AttenIND))& ...
                              (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                     end
                     CONDind=unique(CONDindFULL);
                     if length(CONDind)>1
                        error(sprintf('Non-unique CONDind in EHvNeBFi, when looking at FreqIND=%d and AttenIND=%d',FreqIND,AttenIND))
                     end
                     
                     % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
                     clear TempOctINDs TempAttINDs
                     for i=1:length(picList)
                        TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
                        TempAttINDs(i)=sum(ismember(Temp.Nattens_dB{i},SORTattens(AttenIND)))>0;
                     end
                     TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempAttINDs);  % Mask for Feature,Harm,Polarity,Freq,Atten
                     
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{AttenIND,FreqIND}=Temp.picNums(TempINDs2);
                     % Store excludeLines for each picture, based on CONDind and TOTALconds
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{AttenIND,FreqIND}=cell(size(TempINDs2));
                     ii=0;
                     for i=TempINDs2
                        ii=ii+1;
                        TOTALconds=length(xx{i}.Stimuli.Used.FeatureTarget_Hz_List);
                        GOODlines=CONDind:TOTALconds:xx{i}.Stimuli.fully_presented_lines;
                        %% Take out any lines beyond badlines
                        GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
                        IGNORElines=setdiff(1:xx{i}.Stimuli.fully_presented_lines,GOODlines);
                        yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{AttenIND,FreqIND}{ii}=IGNORElines;
                     end
                  end
               end
            end
         end % End Invert Polarity
      end % End FormsatHarms
      % HERE, 'TN' is saved in the same place
      eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EHvN_reBF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
   end % End: FeatIND
% end 

return; %EHvNrBFi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHvLTASSrBFi_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText

picList=DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.picNums;
if ~isempty(picList)
    for ind=1:length(picList)
        picNUM=picList(ind);
        oldLIST = DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind};
        newLIST = union(oldLIST,DataList.picBADlines{picNUM});
        DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind}=newLIST;
    end
end

NumF=length(SORTocts);
NumL=length(SORTlevels);
NumA=length(SORTattens);

if NumL~=1
    error('There is more than one vowel level in this EHvLTASSreBFI data, makeNOHRdataList NOT SETUP FOR THIS!')
end

%    TONEfeatureIND=find(strcmp(FeaturesText,'TN'));
%% Do all NON-tone features first
for FeatIND=SORTfeatures
    yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
    for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
        for InvertPolarityIND=1:length(InvertPolarityText)
            
            % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
            clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
            for i=1:length(picList)
                TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
            end
            TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
            TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
            TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity
            
            if ~isempty(TempINDs)
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.octshifts=SORTocts;
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=[];
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.Nattens_dB=SORTattens;
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumA,NumF);
                yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumA,NumF);
                
                %%%%%%%%%%% STORE picNUMS for each Freq and Level
                for FreqIND=1:NumF
                    
                    % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
                    CONDindFULL=[];
                    for PICind=TempINDs
                        CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                    end
                    CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
                    if length(CONDfreq_kHz)~=1
                        error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
                    end
                    yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz(FreqIND)=CONDfreq_kHz;
                    
                    for AttenIND=1:NumA
                        % Verify unique conditions
                        CONDindFULL=[];
                        for PICind=TempINDs
                            CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.Nattens_dB{PICind}==SORTattens(AttenIND))& ...
                                (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                        end
                        CONDind=unique(CONDindFULL);
                        if length(CONDind)>1
                            error(sprintf('Non-unique CONDind in EHvLTASSeBFi, when looking at FreqIND=%d and AttenIND=%d',FreqIND,AttenIND))
                        end
                        
                        % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
                        clear TempOctINDs TempAttINDs
                        for i=1:length(picList)
                            TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
                            TempAttINDs(i)=sum(ismember(Temp.Nattens_dB{i},SORTattens(AttenIND)))>0;
                        end
                        TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempAttINDs);  % Mask for Feature,Harm,Polarity,Freq,Atten
                        
                        yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{AttenIND,FreqIND}=Temp.picNums(TempINDs2);
                        % Store excludeLines for each picture, based on CONDind and TOTALconds
                        yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{AttenIND,FreqIND}=cell(size(TempINDs2));
                        ii=0;
                        for i=TempINDs2
                            ii=ii+1;
                            TOTALconds=length(xx{i}.Stimuli.Used.FeatureTarget_Hz_List);
                            GOODlines=CONDind:TOTALconds:xx{i}.Stimuli.fully_presented_lines;
                            %% Take out any lines beyond badlines
                            GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
                            IGNORElines=setdiff(1:xx{i}.Stimuli.fully_presented_lines,GOODlines);
                            yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{AttenIND,FreqIND}{ii}=IGNORElines;
                        end
                    end
                end
            end
        end % End Invert Polarity
    end % End FormsatHarms
    % HERE, 'TN' is saved in the same place
    eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EHvLTASS_reBF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
end % End: FeatIND
return; %EHvLTASSrBFi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=EHrFF_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText FormsAtHarmonicsText InvertPolarityText

picList=findPics('EHrFF',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('EHrFF data Mismatch for Picture: %d',picList(ind)));
      end
      
      %%% frequency comes straight from Stimuli.Condition 
      Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*1000;
      Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
      Temp.FeatureIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.Feature),FeaturesText));
      Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
      Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
   end
   
   SORTfreqs=unique(Temp.freqs_Hz);
   SORTlevels=unique(Temp.levels_dBSPL);
   NumF=length(SORTfreqs);
   NumL=length(SORTlevels);
   SORTfeatures=unique(Temp.FeatureIND);
   
   for FeatIND=SORTfeatures
      yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
      for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
         for InvertPolarityIND=1:length(InvertPolarityText)
            TempINDs=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)&(Temp.InvertPolarityIND==InvertPolarityIND));
            if ~isempty(TempINDs)
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=SORTfreqs/1000;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
               yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
               
               %%%%%%%%%%% STORE picNUMS
               for FreqIND=1:NumF
                  for LevelIND=1:NumL
                     TempINDs2=find((Temp.FeatureIND==FeatIND)&(Temp.FormsAtHarmsIND==FormsAtHarmsIND)& ...
                        (Temp.InvertPolarityIND==InvertPolarityIND)&(Temp.freqs_Hz==SORTfreqs(FreqIND))& ...
                        (Temp.levels_dBSPL==SORTlevels(LevelIND)));
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
                  end
               end                  
               
            end
         end % End Invert Polarity
      end % End FormsatHarms
      eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reFF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
   end % End: FeatIND
end

return;  % EHrFF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=TrBFi_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global FeaturesText

picList=findPics('TrBFi',TUlist(UNITind,:));
if ~isempty(picList)
   if isfield(DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)},'Tone_reBF')
      error('This unit has both Tone_reBF and Tone_reBFi data, makeNOHRdataList IS NOT SET UP TO COMBINE THESE!!!! (NEED TO DO IT NOW)')
   end
   yTEMP.interleaved='yes';
   
   % FIRST, find all octshifts, Features, freqs 
   clear Temp xx
   Temp.octshifts=cell(1,length(picList));
   Temp.FeatureINDs=cell(1,length(picList));
   Temp.freqs_Hz=cell(1,length(picList));
   Temp.FeaturesList=cell(1,length(picList));
   Temp.levels_dBSPL=cell(1,length(picList));
   Temp.FirstBADline=Inf+ones(size(picList));
   for ind=1:length(picList)
      disp(sprintf('      ... Loading picture: %d',picList(ind)))
      x=loadPic(picList(ind));
      xx{ind}=x;
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('TonerBFi data Mismatch for Picture: %d',picList(ind)));
      end
      
      %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
      if ~isempty(x.Stimuli.bad_lines)
         Temp.FirstBADline(ind)=x.Stimuli.bad_lines(1);
         disp(sprintf('Tone_reBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(ind),Temp.FirstBADline(ind)));
      end
      
      Temp.octshifts{ind}=x.Stimuli.Used.OctShifts_List;
      %% PROBABLY NOT NEEDED!!
      Temp.FeatureINDs{ind}=find(strcmp(deblank(x.Stimuli.Used.Features_List{1}),FeaturesText));
      for i=2:length(x.Stimuli.Used.Features_List)
         Temp.FeatureINDs{ind}=[Temp.FeatureINDs{ind} find(strcmp(deblank(x.Stimuli.Used.Features_List{i}),FeaturesText))];
      end
      Temp.freqs_Hz{ind}=x.Stimuli.Used.FeatureTarget_Hz_List;
      Temp.FeaturesList{ind}=x.Stimuli.Used.Features_List;
      %          Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
      Temp.levels_dBSPL{ind}=x.Stimuli.Used.Levels_dBSPL_List;
      Temp.picNums(ind)=picList(ind);
   end
   
   TempFeatureINDs=[];
   Tempoctshifts=[];
   Templevels_dBSPL=[];
   for i=1:length(picList)
      Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
      TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
      Templevels_dBSPL=[Templevels_dBSPL Temp.levels_dBSPL{i}];
   end
   SORTocts=unique(Tempoctshifts);
   SORTfeatures=unique(TempFeatureINDs);
   SORTlevels=unique(Templevels_dBSPL);
   NumF=length(SORTocts);
   NumL=length(SORTlevels);
   
   % Check if TONE only
   if ~isempty(find(SORTfeatures~=find(strcmp('TN',FeaturesText))))
	   beep
	   disp(sprintf('Tone_reBFi picture %d: THIS PICTURE HAS MORE FEATURES THAN JUST TONE (TN), IGNORING OTHER FEATURES',picList(ind)));
   end      
   FeatIND=find(strcmp('TN',FeaturesText));  % LIMIT TO ONLY 'TN'
   
   % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
   clear TempFeatINDs TempINDs
   for i=1:length(picList)
      TempFeatINDs(i)=(sum(ismember(Temp.FeatureINDs{i},FeatIND))>0);
   end
   TempINDs=TempFeatINDs;  % 1 for pictures with current: Feature, Harm, Polarity (left over from EHrBFi)
   
   if ~isempty(TempINDs)
      yTEMP.octshifts=SORTocts;
      yTEMP.freqs_kHz=[];
      yTEMP.levels_dBSPL=SORTlevels;
      yTEMP.picNums=cell(NumL,NumF);
      yTEMP.excludeLines=cell(NumL,NumF);
      
      %%%%%%%%%%% STORE picNUMS for each Freq and Level
      for FreqIND=1:NumF
         
         % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
         CONDindFULL=[];
         for PICind=TempINDs
            CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
         end
         CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
         if length(CONDfreq_kHz)~=1
            error('Non-unique CONDfreq_kHz in Tone_reBFi, when looking at FreqIND=%d',FreqIND)
         end
         yTEMP.freqs_kHz(FreqIND)=CONDfreq_kHz;
         
         for LevelIND=1:NumL
            % Verify unique conditions
            CONDindFULL=[];
            for PICind=TempINDs
               CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.levels_dBSPL{PICind}==SORTlevels(LevelIND))& ...
                     (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
            end
            CONDind=unique(CONDindFULL);
            if length(CONDind)~=1
               error('Non-unique CONDind in TONEreBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND)
            end
            
            % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
            clear TempOctINDs TempLevINDs
            for i=1:length(picList)
               TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
               TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL{i},SORTlevels(LevelIND)))>0;
            end
            TempINDs2=find(TempFeatINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Freq,Level
            
            yTEMP.picNums{LevelIND,FreqIND}=Temp.picNums(TempINDs2);
            % Store excludeLines for each picture, based on CONDind and TOTALconds
            yTEMP.excludeLines{LevelIND,FreqIND}=cell(size(TempINDs2));
            ii=0;
            for i=TempINDs2
               ii=ii+1;
               TOTALconds=length(xx{i}.Stimuli.Used.FeatureTarget_Hz_List);
               GOODlines=CONDind:TOTALconds:xx{i}.Stimuli.fully_presented_lines;
               %% Take out any lines beyond badlines
               GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
               IGNORElines=setdiff(1:xx{i}.Stimuli.fully_presented_lines,GOODlines);
               yTEMP.excludeLines{LevelIND,FreqIND}{ii}=IGNORElines;
            end
         end
      end                  
   end
   eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF=yTEMP;'])
end  

return;  % TONE_reBFi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=TrBF_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind

picList=findPics('TrBF',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('TrBF data Mismatch for Picture: %d',picList(ind)));
      end
      Temp.freqs_Hz(ind)=x.Stimuli.main.tone.freq;
      Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
   end
   SORTfreqs=unique(Temp.freqs_Hz);
   SORTlevels=unique(Temp.levels_dBSPL);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.freqs_kHz=SORTfreqs/1000;
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.levels_dBSPL=SORTlevels;
   NumF=length(SORTfreqs);
   NumL=length(SORTlevels);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.picNums=cell(NumL,NumF);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.excludeLines=cell(NumL,NumF);
   for FreqIND=1:NumF
      for LevelIND=1:NumL
         DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.picNums{LevelIND,FreqIND}= ...
            Temp.picNums(find((Temp.freqs_Hz==SORTfreqs(FreqIND))&(Temp.levels_dBSPL==SORTlevels(LevelIND))));
         DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.excludeLines{LevelIND,FreqIND}= ...
            cell(size(DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reBF.picNums{LevelIND,FreqIND}));
      end
   end
end

return;  % TrBF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=TrFF_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind

picList=findPics('TrFF',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('TrFF data Mismatch for Picture: %d',picList(ind)));
      end
      Temp.freqs_Hz(ind)=x.Stimuli.main.tone.freq;
      Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
   end
   SORTfreqs=unique(Temp.freqs_Hz);
   SORTlevels=unique(Temp.levels_dBSPL);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.freqs_kHz=SORTfreqs/1000;
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.levels_dBSPL=SORTlevels;
   NumF=length(SORTfreqs);
   NumL=length(SORTlevels);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.picNums=cell(NumL,NumF);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.excludeLines=cell(NumL,NumF);
   for FreqIND=1:NumF
      for LevelIND=1:NumL
         DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.picNums{LevelIND,FreqIND}= ...
            Temp.picNums(find((Temp.freqs_Hz==SORTfreqs(FreqIND))&(Temp.levels_dBSPL==SORTlevels(LevelIND))));
         DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.excludeLines{LevelIND,FreqIND}= ...
            cell(size(DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_reFF.picNums{LevelIND,FreqIND}));
      end
   end
end  % Tone_reFF

return;  % TrFF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=TONErlv_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind

picList=findPics('TB',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   for ind=1:length(picList)
      x=loadPic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('TONErlv data Mismatch for Picture: %d',picList(ind)));
      end
      Temp.freqs_Hz(ind)=x.Stimuli.main.tone.freq;
      %          Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
      Temp.picNums(ind)=picList(ind);
   end
   SORTfreqs=unique(Temp.freqs_Hz);
   %       SORTlevels=unique(Temp.levels_dBSPL);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.freqs_kHz=SORTfreqs/1000;
   NumF=length(SORTfreqs);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.picNums=cell(1,NumF);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.excludeLines=cell(1,NumF);
   for FreqIND=1:NumF
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.picNums{FreqIND}= ...
         Temp.picNums(find(Temp.freqs_Hz==SORTfreqs(FreqIND)));
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.excludeLines{FreqIND}= ...
         cell(size(DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Tone_rlv.picNums{FreqIND}));
   end
end 

return;  %Tone_rlv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=SACrlvQ_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CHECK THIS LATER
disp('CHECK THIS LATER')
beep
beep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

picList=DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.picNums;
if ~isempty(picList)
	for ind=1:length(picList)
		picNUM=picList(ind);
		oldLIST = DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind};
		newLIST = union(oldLIST,DataList.picBADlines{picNUM});
		DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.quick.excludeLines{ind}=newLIST;
	end
end

return;  % SACrlvQ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataList=SACrlv_update(DataList) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global TUlist UNITind
global InvertPolarityText

picList=findPics('SACrlv',TUlist(UNITind,:));
if ~isempty(picList)
   clear Temp
   Temp.Nattens_dB=cell(1,length(picList));
   Temp.FirstBADline=Inf+ones(size(picList));
   Temp.picNums=NaN+ones(size(picList));
   Temp.InvertPolarityIND=NaN+ones(size(picList));
   for ind=1:length(picList)
      disp(sprintf('      ... Gathering conditions from picture: %d',picList(ind)))
      x=loadPic(picList(ind));
      xx{ind}=x;
      
      if strcmp(deblank(x.Stimuli.short_description),'SACrlv')
         if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
            error(sprintf('SACrlv data Mismatch for Picture: %d',picList(ind)));
         end
         
         %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
         if ~isempty(x.Stimuli.bad_lines)
            Temp.FirstBADline(ind)=x.Stimuli.bad_lines(1);
            disp(sprintf('SACrlv picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(ind),Temp.FirstBADline(ind)));
         end
         
         Temp.Nattens_dB{ind}=x.Stimuli.attens;
         Temp.picNums(ind)=picList(ind);
         Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
      end
   end
   
   TempNattens_dB=[];
   for i=1:length(picList)
      TempNattens_dB=[TempNattens_dB Temp.Nattens_dB{i}];
   end
   SORTattens=unique(TempNattens_dB);
   NumA=length(SORTattens);
   SORTpolarity=unique(Temp.InvertPolarityIND);
   
   if NumA>1
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.interleaved='yes';
   end
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.Nattens_dB=SORTattens;
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.InvertPolarityText=InvertPolarityText;
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.picNums=cell(length(InvertPolarityText),NumA);
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.excludeLines=cell(length(InvertPolarityText),NumA);
   
   for InvertPolarityIND=1:length(InvertPolarityText)
      for AttenIND=1:NumA
         
         % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
         clear TempPolINDs TempINDs TempAttINDs
         TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
         
         % Find relevant picture indices for each Atten condition (1: yes, 0:no)
         for i=1:length(picList)
            TempAttINDs(i)=sum(ismember(Temp.Nattens_dB{i},SORTattens(AttenIND)))>0;
         end
         TempINDs=find(TempPolINDs&TempAttINDs);  
         
         if ~isempty(TempINDs)
            
            % Find Condition Index
            CONDindFULL=[];
            for PICind=TempINDs
               CONDindFULL=[CONDindFULL find(Temp.Nattens_dB{PICind}==SORTattens(AttenIND))];
            end
            CONDind=unique(CONDindFULL);
            
            DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.picNums{InvertPolarityIND,AttenIND}=Temp.picNums(TempINDs);
            % Store excludeLines for each picture, based on CONDind and TOTALconds
            DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.excludeLines{InvertPolarityIND,AttenIND}=cell(size(TempINDs));
            ii=0;
            for i=TempINDs
               ii=ii+1;
               TOTALconds=length(xx{i}.Stimuli.attens);
               GOODlines=CONDind:TOTALconds:xx{i}.Stimuli.fully_presented_lines;
               %% Take out any lines beyond badlines
               GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
               IGNORElines=setdiff(1:xx{i}.Stimuli.fully_presented_lines,GOODlines);
               DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.SAC.rlv.excludeLines{InvertPolarityIND,AttenIND}{ii}=IGNORElines;
            end
         end
      end % END AttenIND
   end % End InvertPolarityIND
end

return;  % SACrlv

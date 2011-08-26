function makeNOHRdataList(ExpDate)
% File: makeNOHRdataList.m
% Date: 08Sep2004 (M. Heinz)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Modified From: makeCNdataList.m 
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
% Builds basic info data list for each experiment
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO 3/19/04
% *1) Setup saving file at end
% *2) check all stimuli and be sure storing the right thing, and everything we want saved
% *3) Setup re-run, to save fields that need to be kept on re-run
% *4) How to mark which CALIB file goes with which pictures?
%     -CalibFile_ToUse(index=1:MAXpicNum)=Calib_Pic_Num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO 9/8/04
% *1) Setup general utility files into Matlab_Mfiles, from CNexps
% *2) Setup basic Tone_re_BF dataBase
% *3) Setup basic Tone_re_BF Analysis & Plotting
%      *- make general files for analysis
%  4) Setup EHreBF *database, analysis, and plotting
%  5) Setup TrrFF databse, analysis and plotting
%  6) Setup EHreFF databse, analysis and plotting

global NOHR_dir NOHR_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

if ~exist('ExpDate','var')
   ExpDate=0;
ExpDate='111804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(['Processing Experiment: ''' ExpName ''' ... '])

%%%%%% IF File exist already, save old version to recover specific fields at END
if length(dir(FileName))
   disp(['   *** Loading old file: "' FileName '"; specified fields will be saved as is ...'])
   eval(['load ' FileName])
   OLD_DataList=DataList;
   clear DataList
end


%%% General info about Experiment
%%%%%%%% SAVED on re-load %%%%%%%%%%%%%%
DataList.General.date=ExpDateText;
% Status_Text={'Normal','Impaired'};
DataList.General.Status='Normal';
if strcmp(DataList.General.Status,'Impaired')
   DataList.General.ExposureLevel_dBSPL=[];
end
DataList.General.TOTALpics=length(findPics('p','notTAG'));

%%% Ignore Picture List
%%%%%%%% SAVED on re-load %%%%%%%%%%%%%%
DataList.IgnorePicNums=[];

%%% Calibration info
DataList.CALIB.picNUMs=findPics('calib');
%%%%%%%% SAVED on re-load %%%%%%%%%%%%%%
%%%%%%%%%%%%% NEED way to label which CALIB pic is to be used with which pictures
DataList.CALIB.CalibFile_ToUse=NaN*ones(1,DataList.General.TOTALpics);  % TO BE FILLED IN LATER

%%% CAP info
DataList.CAP.picNUM=findPics('CAP1f');
for picIND=DataList.CAP.picNUM
   x=loadPic(picIND);
   DataList.CAP.freq_Hz(find(DataList.CAP.picNUM==picIND))=x.Stimuli.freq_hz;
end
%%% TO SAVE TIME FOR NOW!!!
DataList.CAP.freq_Hz=NaN*ones(size(DataList.CAP.picNUM));

%%%%%%%%%%%%%%%% Get Track/Unit List for Experiment
TUlist=getTrackUnitList;
TOTALunits=size(TUlist,1);
DataList.General.TOTALunits=TOTALunits;
if isempty(TUlist)
   error('******NO UNITS: File not setup yet to handle this!!!!')
end

DataList.Units=cell(max(TUlist));
for UNITind=1:size(TUlist,1)
   
% for UNITind=28
%    beep
%    disp(sprintf('HARD CODED: UNITind: %d',UNITind))
   
   disp(sprintf('   Unit: %d.%02d',TUlist(UNITind,:)))
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% General Unit Info (postAnal)
   %%%%%%%% SAVED on re-load %%%%%%%%%%%%%%
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.ExpName=ExpName; 
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.Unit=sprintf('%d.%02d',TUlist(UNITind,:)); 
   
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.BF_kHz=[];  % For now, replace later with TC data
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.Threshold_dBatten=[];  % For now, replace later with TC data
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.Comment=[];
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.SR_sps=[];
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.Q10=[];
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.TCindToUse=1;
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.QUALITY=[];
  
   %%%%%%%%%%%%%
   %%% Tuning-Curve Info
   %%%%%%%%%%%%%
   picList=findPics('tc',TUlist(UNITind,:));
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.TC=cell(1,length(picList));
   for ind=1:length(picList)
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.TC{ind}.picNUM=picList(ind);
      x=loadpic(picList(ind));
      if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
         error(sprintf('TC data Mismatch for Picture: %d',picList(ind)));
      end
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.TC{ind}.book_BF_kHz=x.Thresh.BF;
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.TC{ind}.book_Thr_dBatt=x.Thresh.thresh;
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.BF_kHz=x.Thresh.BF;  % Use TC data, for cases were UnitInfo is BAD
      DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.Threshold_dBatten=x.Thresh.thresh;
   end
   
   
   %%%%%%%%%%%%%
   %%% Tone_reBF Info
   %%%%%%%%%%%%%
   picList=findPics('TrBF',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      for ind=1:length(picList)
         x=loadpic(picList(ind));
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
   end  % Tone_reBF
   
   %%%%%%%%%%%%%
   %%% EH_reBF Info
   %%%%%%%%%%%%%
   picList=findPics('EHrBF',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      for ind=1:length(picList)
         x=loadpic(picList(ind));
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
         Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
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
         eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBF.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
      end % End: FeatIND
   end  % EH_reBF


   %%%%%%%%%%%%%
   % Added 121304 (M. Heinz)

   %%%%%%%%%%%%%
   %%% EHrlv Info
   %%%%%%%%%%%%%
   picList=findPics('EHrlv',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      for ind=1:length(picList)
         x=loadpic(picList(ind));
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
   end  % EHrlv


   %%%%%%%%%%%%%
   %%% EH_reBFi Info
   %%%%%%%%%%%%%
   picList=findPics('EHrBFi',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      
      Temp.octshifts={};
      Temp.FeatureINDs={};
      Temp.freqs_Hz={};
      Temp.FeaturesList={};
      Temp.FirstBADline=Inf+ones(size(picList));
      for ind=1:length(picList)
         disp(sprintf('      ... Loading picture: %d',picList(ind)))
         x=loadpic(picList(ind));
         xx{ind}=x;
         if (x.General.picture_number~=picList(ind))|(x.General.track~=TUlist(UNITind,1))|(x.General.unit~=TUlist(UNITind,2))
            error(sprintf('EHrBFi data Mismatch for Picture: %d',picList(ind)));
         end

         %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
         if ~isempty(x.Stimuli.bad_lines)
            Temp.FirstBADline(ind)=x.Stimuli.bad_lines;
            disp(sprintf('EHrBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picList(ind),Temp.FirstBADline(ind)));
         end
         
         Temp.octshifts{ind}=x.Stimuli.Used.OctShifts_List;
         Temp.FeatureINDs{length(Temp.FeatureINDs)+1}=find(strcmp(deblank(x.Stimuli.Condition.Features{1}),FeaturesText));
         for i=2:length(x.Stimuli.Condition.Features)
            Temp.FeatureINDs{length(Temp.FeatureINDs)}=[Temp.FeatureINDs{length(Temp.FeatureINDs)} find(strcmp(deblank(x.Stimuli.Condition.Features{i}),FeaturesText))];
         end
         Temp.freqs_Hz{ind}=x.Stimuli.Used.FeatureTarget_Hz_List;
         Temp.FeaturesList{ind}=x.Stimuli.Used.Features_List;
         %          Temp.freqs_Hz(ind)=x.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
         Temp.levels_dBSPL(ind)=x.Stimuli.Condition.Level_dBSPL;
         Temp.picNums(ind)=picList(ind);
         Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
         Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x.Stimuli.Condition.InvertPolarity),InvertPolarityText));
      end
      
      SORTlevels=unique(Temp.levels_dBSPL);
      TempFeatureINDs=[];
      Tempoctshifts=[];
      for i=1:length(picList)
         Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
         TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
      end
      SORTocts=unique(Tempoctshifts);
      SORTfeatures=unique(TempFeatureINDs);
      NumF=length(SORTocts);
      NumL=length(SORTlevels);
      
      for FeatIND=SORTfeatures
         yTEMP=cell(length(FormsAtHarmonicsText),length(InvertPolarityText));
         for FormsAtHarmsIND=1:length(FormsAtHarmonicsText)
            for InvertPolarityIND=1:length(InvertPolarityText)
               
               % Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picList)
               for i=1:length(picList)
                  TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{1},FeatIND));
               end
               TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
               TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
               
               TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % For Featurs, Harm, Polarity
               
               if ~isempty(TempINDs)
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.octshifts=SORTocts;
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz=[];
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.levels_dBSPL=SORTlevels;
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.picNums=cell(NumL,NumF);
                  yTEMP{FormsAtHarmsIND,InvertPolarityIND}.excludeLines=cell(NumL,NumF);
                                    
                  %%%%%%%%%%% STORE picNUMS
                  for FreqIND=1:NumF
                     
                     % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
                     CONDind=[];
                     for PICind=TempINDs
                        CONDind=[CONDind find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
                     end
                     
                     if length(unique(CONDind))~=1
                        error('Non-unique CONDind in EHreBFi, when looking at FreqIND=%d',FreqIND)
                     end
                     yTEMP{FormsAtHarmsIND,InvertPolarityIND}.freqs_kHz(FreqIND)=Temp.freqs_Hz{PICind}(CONDind)/1000;
                     
                     for LevelIND=1:NumL
                        % Find relevant picture indices for each condition (1: yes, 0:no)
                        for i=1:length(picList)
                           TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
                           TempLevINDs(i)=Temp.levels_dBSPL(i)==SORTlevels(LevelIND);
                        end
                        TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);

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
         eval(['DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.EH_reBFi.' FeaturesText{FeatIND}(1:2) '=yTEMP;'])
      end % End: FeatIND
   end  % EH_reBF


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End addition 121304   
   
   %%%%%%%%%%%%%
   %%% Tone_reFF Info
   %%%%%%%%%%%%%
   picList=findPics('TrFF',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      for ind=1:length(picList)
         x=loadpic(picList(ind));
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

   
   %%%%%%%%%%%%%
   %%% EH_reFF Info
   %%%%%%%%%%%%%
   picList=findPics('EHrFF',TUlist(UNITind,:));
   if ~isempty(picList)
      clear Temp
      for ind=1:length(picList)
         x=loadpic(picList(ind));
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
   end  % EH_reBF
   
end

%%%%%% IF File existed already, recover specific fields
if exist('OLD_DataList','var')
   disp(['   *** Saving specified fields as is ...'])
   DataList.CALIB.CalibFile_ToUse=OLD_DataList.CALIB.CalibFile_ToUse;
   %    DataList.Info=OLD_DataList.Info;
   %    DataList.General=OLD_DataList.General;
   DataList.IgnorePicNums=OLD_DataList.IgnorePicNums;
   [Ntracks,Nunits]=size(DataList.Units);
   %    for trackIND=1:Ntracks
   %       for unitIND=1:Nunits
   %          if isfield(DataList.Units{trackIND,unitIND},'Tone_reBF')
   %             DataList.Units{trackIND,unitIND}.Tone_reBF.excludeLines=OLD_DataList.Units{trackIND,unitIND}.Tone_reBF.excludeLines;
   %          end
   %       end
   %    end  
end

disp(['   *** Saving new file: "' FileName '" ...'])
eval(['save ' FileName ' DataList'])
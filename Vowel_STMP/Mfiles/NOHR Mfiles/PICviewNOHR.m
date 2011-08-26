function PICview_NOHR(picNUMs)
% File: PSTview_PickParams.m
% M. Heinz May 19, 2005
% From: makeNOHRdataList.m
%
% pass a set of picNums with the same TAG, program will catalog
% a list of existing conditions and then user chooses which condition to show.
% This is a more user-friendly version of PSTview, which avoids the user having to calculate 'exclude-lines'
%
% Can use to look through whole units at a time: e.g.,
% Usage: PICviewNOHR(findpics('TrBF',[1,12]))
%
%

global FeaturesText FormsAtHarmonicsText InvertPolarityText

if isempty(picNUMs)
   error('NO PICTURES IN LIST!!!')
end
for i=1:length(picNUMs)
   x{i}=loadPic(picNUMs(i));
   TAGs{i}=getTAG(getFileName(picNUMs(i)));
end
% Verify TAGs match across PICs, but then just get all conditions
for i=2:length(picNUMs)
   if ~strcmp(TAGs{i},TAGs{1})
      % Allow TB and EHrlv to be plotted together
      if strcmp(TAGs{1},'EHrlv')&strcmp(TAGs{i},'TB')
      elseif strcmp(TAGs{1},'TB')&strcmp(TAGs{i},'EHrlv')
      else
         error('TAGs don''t match')
      end   
   end
end

if strcmp(TAGs{1},'EHrlv')
   EHrlv_parse(picNUMs,x)
elseif strcmp(TAGs{1},'EHINrlv')
   EHINrlv_parse(picNUMs,x)
elseif strcmp(TAGs{1},'EHrBFi')
   EHrBFi_parse(picNUMs,x)
elseif strcmp(TAGs{1},'EHvNrBFi')
   EHvNrBFi_parse(picNUMs,x)
elseif strcmp(TAGs{1},'EHrBF')
   EHrBF_parse(picNUMs,x)
elseif strcmp(TAGs{1},'EHrFF')
   EHrFF_parse(picNUMs,x)
elseif strcmp(TAGs{1},'TrBF')
   TrBF_parse(picNUMs,x)
elseif strcmp(TAGs{1},'TrBFi')
   TrBFi_parse(picNUMs,x)
elseif strcmp(TAGs{1},'TrFF')
   TrFF_parse(picNUMs,x)
elseif strcmp(TAGs{1},'TB')
   EHrlv_parse(picNUMs,x)

%%%%%%%%%%%%%% TODO 
elseif strcmp(TAGs{1},'SACrlvQ')
   EHrlv_parse(picNUMs,x)
elseif strcmp(TAGs{1},'SACrlv')
   SACrlv_parse(picNUMs,x)
   
   
   
   
else
   beep
   sprintf('TAG: %s, NOT SET UP FOR PICviewNOHR!!!',TAGs{1})
end

%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHrlv_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% Try to get CalibFile form DataList
DataFile=dir('DataList*');
if ~isempty(DataFile)
   eval(['load ' DataFile.name])
   %       TrackUnit=getTrackUnit(getFileName(picNUMs(1)));
   if length(unique(DataList.CALIB.CalibFile_ToUse(picNUMs)))~=1
      error('Different CALIB files for these picNUMs')
   end
   CALIBpic=DataList.CALIB.CalibFile_ToUse(picNUMs(1));
else
   CALIBpicNUMs=findPics('calib');
   if length(CALIBpicNUMs)==1
      CALIBpic=CALIBpicNUMs; 
   else
      CALIBpic=-1;
      while isempty(find(CALIBpicNUMs==CALIBpic))
         CALIBpic=input(sprintf('Pick CALIBpic: [%s]: ',mat2str(CALIBpicNUMs)));
      end
   end
end

quick_EHrlfs(picNUMs,CALIBpic);   
return;  % END EHrlv


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHINrlv_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
% Try to get CalibFile form DataList
DataFile=dir('DataList*');
if ~isempty(DataFile)
   eval(['load ' DataFile.name])
   %       TrackUnit=getTrackUnit(getFileName(picNUMs(1)));
   if length(unique(DataList.CALIB.CalibFile_ToUse(picNUMs)))~=1
      error('Different CALIB files for these picNUMs')
   end
   CALIBpic=DataList.CALIB.CalibFile_ToUse(picNUMs(1));
else
   CALIBpicNUMs=findPics('calib');
   if length(CALIBpicNUMs)==1
      CALIBpic=CALIBpicNUMs; 
   else
      CALIBpic=-1;
      while isempty(find(CALIBpicNUMs==CALIBpic))
         CALIBpic=input(sprintf('Pick CALIBpic: [%s]: ',mat2str(CALIBpicNUMs)));
      end
   end
end

quick_EHINrlfs(picNUMs,CALIBpic);

%% Also run EHrlv
TU=getTrackUnit(getFileName(picNUMs(1)));
PICviewNOHR(findPics('EHrlv',[TU(1),TU(2)]))

return;  % END EHrlv



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHrBFi_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FeaturesText FormsAtHarmonicsText InvertPolarityText

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.octshifts=cell(1,length(picNUMs));
Temp.FeatureINDs=cell(1,length(picNUMs));
Temp.freqs_Hz=cell(1,length(picNUMs));
Temp.FeaturesList=cell(1,length(picNUMs));
Temp.levels_dBSPL=cell(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.FormsAtHarmsIND=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBFi Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
   if ~isempty(x{ind}.Stimuli.bad_lines)
      Temp.FirstBADline(ind)=x{ind}.Stimuli.bad_lines(1);
      disp(sprintf('EHrBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picNUMs(ind),Temp.FirstBADline(ind)));
   end
   
   Temp.octshifts{ind}=x{ind}.Stimuli.Used.OctShifts_List;
   Temp.FeatureINDs{ind}=find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{1}),FeaturesText));
   for i=2:length(x{ind}.Stimuli.Condition.Features)
      Temp.FeatureINDs{ind}=[Temp.FeatureINDs{ind} find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{i}),FeaturesText))];
   end
   Temp.freqs_Hz{ind}=x{ind}.Stimuli.Used.FeatureTarget_Hz_List;
   Temp.FeaturesList{ind}=x{ind}.Stimuli.Used.Features_List;
   %          Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
   if isfield(x{ind}.Stimuli.Used,'Levels_dBSPL_List')
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Used.Levels_dBSPL_List;
   else
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Condition.Level_dBSPL;
   end
   Temp.picNums(ind)=picNUMs(ind);
   Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

TempFeatureINDs=[];
Tempoctshifts=[];
Templevels_dBSPL=[];
for i=1:length(picNUMs)
   Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
   TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
   Templevels_dBSPL=[Templevels_dBSPL Temp.levels_dBSPL{i}];
end
SORTocts=unique(Tempoctshifts);
SORTfeatures=unique(TempFeatureINDs);
SORTlevels=unique(Templevels_dBSPL);
NumF=length(SORTocts);
NumL=length(SORTlevels);
SORTpolarity=unique(Temp.InvertPolarityIND);
SORTharmonics=unique(Temp.FormsAtHarmsIND);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Feature
FeatIND=SORTfeatures(1);
if length(SORTfeatures)>1
   Ftext='';
   for i=1:length(SORTfeatures)
      Ftext=sprintf('%s %s[%d],',Ftext,FeaturesText{SORTfeatures(i)},SORTfeatures(i));
   end
   Ftext=Ftext(2:end-1);
   FeatIND=-1;
   while isempty(find(SORTfeatures==FeatIND))
      FeatIND=input(sprintf('Pick feature:  %s: ',Ftext));
   end
end

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Octave Shift
FreqIND=1;
if length(SORTocts)>1
   Otext='';
   for i=1:length(SORTocts)
      Otext=sprintf('%s %.2f[%d],',Otext,SORTocts(i),i);
   end
   Otext=Otext(2:end-1);
   FreqIND=-1;
   while isempty(find((1:length(SORTocts))==FreqIND))
      FreqIND=input(sprintf('Pick octave shift:  %s: ',Otext));
   end
end

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

% Formants at Harmonics
FormsAtHarmsIND=SORTharmonics(1);
if length(SORTharmonics)>1
   Htext='';
   for i=1:length(SORTharmonics)
      Htext=sprintf('%s [1], %s [2]',FormsAtHarmonicsText{1},FormsAtHarmonicsText{2});
   end
   FormsAtHarmsIND=-1;
   while isempty(find((1:length(SORTharmonics))==FormsAtHarmsIND))
      FormsAtHarmsIND=input(sprintf('Pick Formants-at-Harmonics:   %s: ',Htext));
   end
end

disp(sprintf('   *** Analyzing Condition: %s at BF+ %.2f octaves, @ %.f dB SPL, Polarity Inversion: %s, Formants at Harmonics: %s', ...
   FeaturesText{FeatIND},SORTocts(FreqIND),SORTlevels(LevelIND),upper(InvertPolarityText{InvertPolarityIND}),upper(FormsAtHarmonicsText{FormsAtHarmsIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
for i=1:length(picNUMs)
   TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
end
TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
   end
   CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
   if length(CONDfreq_kHz)~=1
      error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
   end
   
   % Verify unique conditions
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.levels_dBSPL{PICind}==SORTlevels(LevelIND))& ...
            (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
   end
   CONDind=unique(CONDindFULL);
   if length(CONDind)~=1
      error('Non-unique CONDind in TEHeBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND)
   end
   
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempOctINDs TempLevINDs
   for i=1:length(picNUMs)
      TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL{i},SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Freq,Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      TOTALconds=length(x{i}.Stimuli.Used.FeatureTarget_Hz_List);
      GOODlines=CONDind:TOTALconds:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview_simFF(picNums,excludeLines,2);
   
else
   error(upper('No picture with these params!'))
end
return; % END EHreBFi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHvNrBFi_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FeaturesText FormsAtHarmonicsText InvertPolarityText

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.octshifts=cell(1,length(picNUMs));
Temp.FeatureINDs=cell(1,length(picNUMs));
Temp.freqs_Hz=cell(1,length(picNUMs));
Temp.FeaturesList=cell(1,length(picNUMs));
Temp.levels_dBSPL=cell(1,length(picNUMs));
Temp.Nattens_dB=cell(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.FormsAtHarmsIND=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBFi Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
   if ~isempty(x{ind}.Stimuli.bad_lines)
      Temp.FirstBADline(ind)=x{ind}.Stimuli.bad_lines(1);
      disp(sprintf('EHrBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picNUMs(ind),Temp.FirstBADline(ind)));
   end
   
   Temp.octshifts{ind}=x{ind}.Stimuli.Used.OctShifts_List;
   Temp.FeatureINDs{ind}=find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{1}),FeaturesText));
   for i=2:length(x{ind}.Stimuli.Condition.Features)
      Temp.FeatureINDs{ind}=[Temp.FeatureINDs{ind} find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{i}),FeaturesText))];
   end
   Temp.freqs_Hz{ind}=x{ind}.Stimuli.Used.FeatureTarget_Hz_List;
   Temp.FeaturesList{ind}=x{ind}.Stimuli.Used.Features_List;
   %          Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
   if isfield(x{ind}.Stimuli.Used,'Levels_dBSPL_List')
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Used.Levels_dBSPL_List;
   else
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Condition.Level_dBSPL;
   end
   Temp.Nattens_dB{ind}=x{ind}.Stimuli.Used.NoiseAttens_dB_List;
   Temp.picNums(ind)=picNUMs(ind);
   Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

TempFeatureINDs=[];
Tempoctshifts=[];
Templevels_dBSPL=[];
TempNattens_dB=[];
for i=1:length(picNUMs)
   Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
   TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
   Templevels_dBSPL=[Templevels_dBSPL Temp.levels_dBSPL{i}];
   TempNattens_dB=[TempNattens_dB Temp.Nattens_dB{i}];
end
SORTocts=unique(Tempoctshifts);
SORTfeatures=unique(TempFeatureINDs);
SORTlevels=unique(Templevels_dBSPL);
SORTattens=unique(TempNattens_dB);
NumF=length(SORTocts);
NumL=length(SORTlevels);
NumA=length(SORTattens);
SORTpolarity=unique(Temp.InvertPolarityIND);
SORTharmonics=unique(Temp.FormsAtHarmsIND);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Feature
FeatIND=SORTfeatures(1);
if length(SORTfeatures)>1
   Ftext='';
   for i=1:length(SORTfeatures)
      Ftext=sprintf('%s %s[%d],',Ftext,FeaturesText{SORTfeatures(i)},SORTfeatures(i));
   end
   Ftext=Ftext(2:end-1);
   FeatIND=-1;
   while isempty(find(SORTfeatures==FeatIND))
      FeatIND=input(sprintf('Pick feature:  %s: ',Ftext));
   end
end

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Noise Attens
AttenIND=1;
if length(SORTattens)>1
   Atext='';
   for i=1:length(SORTattens)
      Atext=sprintf('%s %d dB [%d],',Atext,SORTattens(i),i);
   end
   Atext=Atext(2:end-1);
   AttenIND=-1;
   while isempty(find(find(SORTattens)==AttenIND))
      AttenIND=input(sprintf('Pick noise atten:  %s: ',Atext));
   end
end

% Octave Shift
FreqIND=1;
if length(SORTocts)>1
   Otext='';
   for i=1:length(SORTocts)
      Otext=sprintf('%s %.2f[%d],',Otext,SORTocts(i),i);
   end
   Otext=Otext(2:end-1);
   FreqIND=-1;
   while isempty(find((1:length(SORTocts))==FreqIND))
      FreqIND=input(sprintf('Pick octave shift:  %s: ',Otext));
   end
end

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

% Formants at Harmonics
FormsAtHarmsIND=SORTharmonics(1);
if length(SORTharmonics)>1
   Htext='';
   for i=1:length(SORTharmonics)
      Htext=sprintf('%s [1], %s [2]',FormsAtHarmonicsText{1},FormsAtHarmonicsText{2});
   end
   FormsAtHarmsIND=-1;
   while isempty(find((1:length(SORTharmonics))==FormsAtHarmsIND))
      FormsAtHarmsIND=input(sprintf('Pick Formants-at-Harmonics:   %s: ',Htext));
   end
end

disp(sprintf('   *** Analyzing Condition: %s at BF+ %.2f octaves, vowel @ %.f dB SPL, noise @ %.f dB atten, Polarity Inversion: %s, Formants at Harmonics: %s', ...
   FeaturesText{FeatIND},SORTocts(FreqIND),SORTlevels(LevelIND),SORTattens(AttenIND),upper(InvertPolarityText{InvertPolarityIND}),upper(FormsAtHarmonicsText{FormsAtHarmsIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
for i=1:length(picNUMs)
   TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
end
TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
   end
   CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
   if length(CONDfreq_kHz)~=1
      error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
   end
   
   % Verify unique conditions
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.levels_dBSPL{PICind}==SORTlevels(LevelIND))& ...
            (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND}))&(Temp.Nattens_dB{PICind}==SORTattens(AttenIND)))];
   end
   CONDind=unique(CONDindFULL);
   if length(CONDind)~=1
      error('Non-unique CONDind in TEHeBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND)
   end
   
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempOctINDs TempLevINDs TempAttINDs
   for i=1:length(picNUMs)
      TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL{i},SORTlevels(LevelIND)))>0;
      TempAttINDs(i)=sum(ismember(Temp.Nattens_dB{i},SORTattens(AttenIND)))>0;
   end
   TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs&TempAttINDs);  % Mask for Feature,Harm,Polarity,Freq,Level,Natten
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      TOTALconds=length(x{i}.Stimuli.Used.FeatureTarget_Hz_List);
      GOODlines=CONDind:TOTALconds:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview_simFF(picNums,excludeLines,2);
   
else
   error(upper('No picture with these params!'))
end
return; % END EHvNreBFi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHrBF_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FeaturesText FormsAtHarmonicsText InvertPolarityText

%% Only reason to include all lists, is for lists of picNUMs, may want to do this sometime.

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.octshifts=zeros(1,length(picNUMs));
Temp.FeatureINDs=zeros(1,length(picNUMs));
Temp.freqs_Hz=zeros(1,length(picNUMs));
Temp.FeaturesList=cell(1,length(picNUMs));
Temp.levels_dBSPL=zeros(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.FormsAtHarmsIND=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBF Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %%% CALC frequency from Stimuli.Condition (b/c this will be the same for all Features, even if not exactly the right frequency)
   OffsetSign=sign(find(strcmp(x{ind}.Stimuli.Condition.Offset_Direction,{'below ','above '}))-1.5);
   if isstr(x{ind}.Stimuli.Condition.FreqOffset_octs)
      Offset_octs=str2num(x{ind}.Stimuli.Condition.FreqOffset_octs);
   else
      Offset_octs=x{ind}.Stimuli.Condition.FreqOffset_octs;
   end
   Temp.octshifts(ind)=OffsetSign*Offset_octs;
   Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
   Temp.levels_dBSPL(ind)=x{ind}.Stimuli.Condition.Level_dBSPL;
   Temp.picNums(ind)=picNUMs(ind);
   Temp.FeatureINDs(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.Feature),FeaturesText));
   Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

SORTocts=unique(Temp.octshifts);
SORTfreqs=unique(Temp.freqs_Hz);
SORTlevels=unique(Temp.levels_dBSPL);
SORTpolarity=unique(Temp.InvertPolarityIND);
SORTharmonics=unique(Temp.FormsAtHarmsIND);
NumF=length(SORTfreqs);
NumL=length(SORTlevels);
SORTfeatures=unique(Temp.FeatureINDs);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Feature
FeatIND=SORTfeatures(1);
if length(SORTfeatures)>1
   Ftext='';
   for i=1:length(SORTfeatures)
      Ftext=sprintf('%s %s[%d],',Ftext,FeaturesText{SORTfeatures(i)},SORTfeatures(i));
   end
   Ftext=Ftext(2:end-1);
   FeatIND=-1;
   while isempty(find(SORTfeatures==FeatIND))
      FeatIND=input(sprintf('Pick feature:  %s: ',Ftext));
   end
end

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Octave Shift
FreqIND=1;
if length(SORTocts)>1
   Otext='';
   for i=1:length(SORTocts)
      Otext=sprintf('%s %.2f[%d],',Otext,SORTocts(i),i);
   end
   Otext=Otext(2:end-1);
   FreqIND=-1;
   while isempty(find((1:length(SORTocts))==FreqIND))
      FreqIND=input(sprintf('Pick octave shift:  %s: ',Otext));
   end
end

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

% Formants at Harmonics
FormsAtHarmsIND=SORTharmonics(1);
if length(SORTharmonics)>1
   Htext='';
   for i=1:length(SORTharmonics)
      Htext=sprintf('%s [1], %s [2]',FormsAtHarmonicsText{1},FormsAtHarmonicsText{2});
   end
   FormsAtHarmsIND=-1;
   while isempty(find((1:length(SORTharmonics))==FormsAtHarmsIND))
      FormsAtHarmsIND=input(sprintf('Pick Formants-at-Harmonics:   %s: ',Htext));
   end
end

disp(sprintf('   *** Analyzing Condition: %s at BF+ %.2f octaves, @ %.f dB SPL, Polarity Inversion: %s, Formants at Harmonics: %s', ...
   FeaturesText{FeatIND},SORTocts(FreqIND),SORTlevels(LevelIND),upper(InvertPolarityText{InvertPolarityIND}),upper(FormsAtHarmonicsText{FormsAtHarmsIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
for i=1:length(picNUMs)
   TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs(i),FeatIND));
end
TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempOctINDs TempLevINDs
   for i=1:length(picNUMs)
      TempOctINDs(i)=sum(ismember(Temp.octshifts(i),SORTocts(FreqIND)))>0;
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL(i),SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Freq,Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      GOODlines=1:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview_simFF(picNums,excludeLines,2);
   
else
   error(upper('No picture with these params!'))
end
return;  % END EHreBF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EHrFF_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FeaturesText FormsAtHarmonicsText InvertPolarityText

%% Only reason to include all lists, is for lists of picNUMs, may want to do this sometime.

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.FeatureINDs=zeros(1,length(picNUMs));
Temp.freqs_Hz=zeros(1,length(picNUMs));
Temp.FeaturesList=cell(1,length(picNUMs));
Temp.levels_dBSPL=zeros(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.FormsAtHarmsIND=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBF Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*1000;
   Temp.levels_dBSPL(ind)=x{ind}.Stimuli.Condition.Level_dBSPL;
   Temp.picNums(ind)=picNUMs(ind);
   Temp.FeatureINDs(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.Feature),FeaturesText));
   Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

SORTfreqs=unique(Temp.freqs_Hz);
SORTlevels=unique(Temp.levels_dBSPL);
SORTpolarity=unique(Temp.InvertPolarityIND);
SORTharmonics=unique(Temp.FormsAtHarmsIND);
NumF=length(SORTfreqs);
NumL=length(SORTlevels);
SORTfeatures=unique(Temp.FeatureINDs);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Feature
FeatIND=SORTfeatures(1);
if length(SORTfeatures)>1
   Ftext='';
   for i=1:length(SORTfeatures)
      Ftext=sprintf('%s %s[%d],',Ftext,FeaturesText{SORTfeatures(i)},SORTfeatures(i));
   end
   Ftext=Ftext(2:end-1);
   FeatIND=-1;
   while isempty(find(SORTfeatures==FeatIND))
      FeatIND=input(sprintf('Pick feature:  %s: ',Ftext));
   end
end

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Frequencies
if NumF~=1
   error('These pictures have more than one Fixed-Frequency for EH_reFF!!!')
end
FreqIND=1;

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

% Formants at Harmonics
FormsAtHarmsIND=SORTharmonics(1);
if length(SORTharmonics)>1
   Htext='';
   for i=1:length(SORTharmonics)
      Htext=sprintf('%s [1], %s [2]',FormsAtHarmonicsText{1},FormsAtHarmonicsText{2});
   end
   FormsAtHarmsIND=-1;
   while isempty(find((1:length(SORTharmonics))==FormsAtHarmsIND))
      FormsAtHarmsIND=input(sprintf('Pick Formants-at-Harmonics:   %s: ',Htext));
   end
end

disp(sprintf('   *** Analyzing Condition: %s at FF=%.2f Hz, @ %.f dB SPL, Polarity Inversion: %s, Formants at Harmonics: %s', ...
   FeaturesText{FeatIND},SORTfreqs(FreqIND),SORTlevels(LevelIND),upper(InvertPolarityText{InvertPolarityIND}),upper(FormsAtHarmonicsText{FormsAtHarmsIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
for i=1:length(picNUMs)
   TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs(i),FeatIND));
end
TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity

%%%%%%% HERE,
% take our SORTocts to SORTfreqs


if ~isempty(TempINDs)
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempLevINDs
   for i=1:length(picNUMs)
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL(i),SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      GOODlines=1:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview(picNums,excludeLines,1);
   
else
   error(upper('No picture with these params!'))
end
return;  % END EHreFF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrBFi_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FeaturesText FormsAtHarmonicsText InvertPolarityText

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.octshifts=cell(1,length(picNUMs));
Temp.FeatureINDs=cell(1,length(picNUMs));
Temp.freqs_Hz=cell(1,length(picNUMs));
Temp.FeaturesList=cell(1,length(picNUMs));
Temp.levels_dBSPL=cell(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.FormsAtHarmsIND=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBFi Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
   if ~isempty(x{ind}.Stimuli.bad_lines)
      Temp.FirstBADline(ind)=x{ind}.Stimuli.bad_lines(1);
      disp(sprintf('EHrBFi picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picNUMs(ind),Temp.FirstBADline(ind)));
   end
   
   Temp.octshifts{ind}=x{ind}.Stimuli.Used.OctShifts_List;
   Temp.FeatureINDs{ind}=find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{1}),FeaturesText));
   for i=2:length(x{ind}.Stimuli.Condition.Features)
      Temp.FeatureINDs{ind}=[Temp.FeatureINDs{ind} find(strcmp(deblank(x{ind}.Stimuli.Condition.Features{i}),FeaturesText))];
   end
   Temp.freqs_Hz{ind}=x{ind}.Stimuli.Used.FeatureTarget_Hz_List;
   Temp.FeaturesList{ind}=x{ind}.Stimuli.Used.Features_List;
   %          Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
   if isfield(x{ind}.Stimuli.Used,'Levels_dBSPL_List')
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Used.Levels_dBSPL_List;
   else
      Temp.levels_dBSPL{ind}=x{ind}.Stimuli.Condition.Level_dBSPL;
   end
   Temp.picNums(ind)=picNUMs(ind);
   Temp.FormsAtHarmsIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.FormsAtHarmonics),FormsAtHarmonicsText));
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

TempFeatureINDs=[];
Tempoctshifts=[];
Templevels_dBSPL=[];
for i=1:length(picNUMs)
   Tempoctshifts=[Tempoctshifts Temp.octshifts{i}];
   TempFeatureINDs=[TempFeatureINDs Temp.FeatureINDs{i}];
   Templevels_dBSPL=[Templevels_dBSPL Temp.levels_dBSPL{i}];
end
SORTocts=unique(Tempoctshifts);
SORTfeatures=unique(TempFeatureINDs);
SORTlevels=unique(Templevels_dBSPL);
NumF=length(SORTocts);
NumL=length(SORTlevels);
SORTpolarity=unique(Temp.InvertPolarityIND);
SORTharmonics=unique(Temp.FormsAtHarmsIND);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Feature
FeatIND=SORTfeatures(1);
if length(SORTfeatures)>1
   Ftext='';
   for i=1:length(SORTfeatures)
      Ftext=sprintf('%s %s[%d],',Ftext,FeaturesText{SORTfeatures(i)},SORTfeatures(i));
   end
   Ftext=Ftext(2:end-1);
   FeatIND=-1;
   while isempty(find(SORTfeatures==FeatIND))
      FeatIND=input(sprintf('Pick feature:  %s: ',Ftext));
   end
end

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Octave Shift
FreqIND=1;
if length(SORTocts)>1
   Otext='';
   for i=1:length(SORTocts)
      Otext=sprintf('%s %.2f[%d],',Otext,SORTocts(i),i);
   end
   Otext=Otext(2:end-1);
   FreqIND=-1;
   while isempty(find((1:length(SORTocts))==FreqIND))
      FreqIND=input(sprintf('Pick octave shift:  %s: ',Otext));
   end
end

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

% Formants at Harmonics
FormsAtHarmsIND=SORTharmonics(1);
if length(SORTharmonics)>1
   Htext='';
   for i=1:length(SORTharmonics)
      Htext=sprintf('%s [1], %s [2]',FormsAtHarmonicsText{1},FormsAtHarmonicsText{2});
   end
   FormsAtHarmsIND=-1;
   while isempty(find((1:length(SORTharmonics))==FormsAtHarmsIND))
      FormsAtHarmsIND=input(sprintf('Pick Formants-at-Harmonics:   %s: ',Htext));
   end
end

disp(sprintf('   *** Analyzing Condition: %s at BF+ %.2f octaves, @ %.f dB SPL, Polarity Inversion: %s', ...
   FeaturesText{FeatIND},SORTocts(FreqIND),SORTlevels(LevelIND),upper(InvertPolarityText{InvertPolarityIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempFeatINDs TempHarmINDs TempPolINDs TempINDs
for i=1:length(picNUMs)
   TempFeatINDs(i)=sum(ismember(Temp.FeatureINDs{i},FeatIND));
end
TempHarmINDs=(Temp.FormsAtHarmsIND==FormsAtHarmsIND);
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);
TempINDs=find(TempFeatINDs&TempHarmINDs&TempPolINDs);  % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Store exact freqs for ALL conditions for this Feature (will vary from feature to feature)
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
   end
   CONDfreq_kHz=unique(Temp.freqs_Hz{PICind}(CONDindFULL)/1000);
   if length(CONDfreq_kHz)~=1
      error('Non-unique CONDfreq_kHz in TEH_reBFi, when looking at FreqIND=%d',FreqIND)
   end
   
   % Verify unique conditions
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.octshifts{PICind}==SORTocts(FreqIND))&(Temp.levels_dBSPL{PICind}==SORTlevels(LevelIND))& ...
            (strcmp(deblank(Temp.FeaturesList{PICind}),FeaturesText{FeatIND})))];
   end
   CONDind=unique(CONDindFULL);
   if length(CONDind)~=1
      error('Non-unique CONDind in TEHeBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND)
   end
   
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempOctINDs TempLevINDs
   for i=1:length(picNUMs)
      TempOctINDs(i)=sum(ismember(Temp.octshifts{i},SORTocts(FreqIND)))>0;
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL{i},SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempFeatINDs&TempHarmINDs&TempPolINDs&TempOctINDs&TempLevINDs);  % Mask for Feature,Harm,Polarity,Freq,Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      TOTALconds=length(x{i}.Stimuli.Used.FeatureTarget_Hz_List);
      GOODlines=CONDind:TOTALconds:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview_simFF(picNums,excludeLines,2);
   
else
   error(upper('No picture with these params!'))
end
% END TreBFi
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrBF_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FeaturesText FormsAtHarmonicsText InvertPolarityText

%% Only reason to include all lists, is for lists of picNUMs, may want to do this sometime.

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.octshifts=zeros(1,length(picNUMs));
Temp.freqs_Hz=zeros(1,length(picNUMs));
Temp.levels_dBSPL=zeros(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('TrBF Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %%% CALC frequency from Stimuli.Condition (b/c this will be the same for all Features, even if not exactly the right frequency)
   OffsetSign=sign(find(strcmp(x{ind}.Stimuli.Condition.Offset_Direction,{'below ','above '}))-1.5);
   if isstr(x{ind}.Stimuli.Condition.FreqOffset_octs)
      Offset_octs=str2num(x{ind}.Stimuli.Condition.FreqOffset_octs);
   else
      Offset_octs=x{ind}.Stimuli.Condition.FreqOffset_octs;
   end
   Temp.octshifts(ind)=OffsetSign*Offset_octs;
   Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*2^(OffsetSign*Offset_octs)*1000;
   Temp.levels_dBSPL(ind)=x{ind}.Stimuli.Condition.Level_dBSPL;
   Temp.picNums(ind)=picNUMs(ind);
end

SORTocts=unique(Temp.octshifts);
SORTfreqs=unique(Temp.freqs_Hz);
SORTlevels=unique(Temp.levels_dBSPL);
NumF=length(SORTfreqs);
NumL=length(SORTlevels);

%%%%%%%%% ASK WHICH CONDITION TO SHOW

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Octave Shift
FreqIND=1;
if length(SORTocts)>1
   Otext='';
   for i=1:length(SORTocts)
      Otext=sprintf('%s %.2f[%d],',Otext,SORTocts(i),i);
   end
   Otext=Otext(2:end-1);
   FreqIND=-1;
   while isempty(find((1:length(SORTocts))==FreqIND))
      FreqIND=input(sprintf('Pick octave shift:  %s: ',Otext));
   end
end

disp(sprintf('   *** Analyzing Condition: Tone at BF+ %.2f octaves, @ %.f dB SPL',SORTocts(FreqIND),SORTlevels(LevelIND)))


% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempINDs
TempINDs=ones(size(picNUMs));    % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempOctINDs TempLevINDs
   for i=1:length(picNUMs)
      TempOctINDs(i)=sum(ismember(Temp.octshifts(i),SORTocts(FreqIND)))>0;
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL(i),SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempINDs&TempOctINDs&TempLevINDs);  % Mask for Freq,Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      GOODlines=1:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview_simFF(picNums,excludeLines,2);
   
else
   error(upper('No picture with these params!'))
end
return;  % END TreBF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TrFF_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Only reason to include all lists, is for lists of picNUMs, may want to do this sometime.

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.freqs_Hz=zeros(1,length(picNUMs));
Temp.levels_dBSPL=zeros(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('TrFF Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   Temp.freqs_Hz(ind)=x{ind}.Stimuli.Condition.BaseFrequency_kHz*1000;
   Temp.levels_dBSPL(ind)=x{ind}.Stimuli.Condition.Level_dBSPL;
   Temp.picNums(ind)=picNUMs(ind);
end

SORTfreqs=unique(Temp.freqs_Hz);
SORTlevels=unique(Temp.levels_dBSPL);
NumF=length(SORTfreqs);
NumL=length(SORTlevels);

%%%%%%%%% ASK WHICH CONDITION TO SHOW

% Level
LevelIND=1;
if length(SORTlevels)>1
   Ltext='';
   for i=1:length(SORTlevels)
      Ltext=sprintf('%s %d dB SPL[%d],',Ltext,SORTlevels(i),i);
   end
   Ltext=Ltext(2:end-1);
   LevelIND=-1;
   while isempty(find(find(SORTlevels)==LevelIND))
      LevelIND=input(sprintf('Pick level:  %s: ',Ltext));
   end
end

% Frequencies
if NumF~=1
   error('These pictures have more than one Fixed-Frequency for Tone_reFF!!!')
end
FreqIND=1;

disp(sprintf('   *** Analyzing Condition: Tone at FF = %.2f Hz, @ %.f dB SPL',SORTfreqs(FreqIND),SORTlevels(LevelIND)))


% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempINDs
TempINDs=ones(size(picNUMs));    % 1 for pictures with current: Feature, Harm, Polarity

if ~isempty(TempINDs)
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempLevINDs
   for i=1:length(picNUMs)
      TempLevINDs(i)=sum(ismember(Temp.levels_dBSPL(i),SORTlevels(LevelIND)))>0;
   end
   TempINDs2=find(TempINDs&TempLevINDs);  % Mask for Level
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      GOODlines=1:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   PSTview(picNums,excludeLines,1);
   
else
   error(upper('No picture with these params!'))
end
return;  % END TreFF




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SACrlv_parse(picNUMs,x) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FeaturesText FormsAtHarmonicsText InvertPolarityText

% FIRST, find all octshifts, Features, freqs 
clear Temp
Temp.Nattens_dB=cell(1,length(picNUMs));
Temp.FirstBADline=Inf+ones(size(picNUMs));
Temp.picNums=NaN+ones(size(picNUMs));
Temp.InvertPolarityIND=NaN+ones(size(picNUMs));
for ind=1:length(picNUMs)
   disp(sprintf('      ... Gathering conditions from picture: %d',picNUMs(ind)))
   if (x{ind}.General.picture_number~=picNUMs(ind))
      error(sprintf('EHrBFi Mismatch for Picture: %d',picNUMs(ind)));
   end
   
   %% Need to verify no bad-line errors, ow, counting will be off - SO FOR NOW, DISCARD ALL LINES AFTER THIS
   if ~isempty(x{ind}.Stimuli.bad_lines)
      Temp.FirstBADline(ind)=x{ind}.Stimuli.bad_lines(1);
      disp(sprintf('SACrlv picture %d: LINE ERRORS, IGNORING ALL LINES STARTING at line #: %d',picNUMs(ind),Temp.FirstBADline(ind)));
   end
   
   Temp.Nattens_dB{ind}=x{ind}.Stimuli.attens;
   Temp.picNums(ind)=picNUMs(ind);
   Temp.InvertPolarityIND(ind)=find(strcmp(deblank(x{ind}.Stimuli.Condition.InvertPolarity),InvertPolarityText));
end

TempNattens_dB=[];
for i=1:length(picNUMs)
   TempNattens_dB=[TempNattens_dB Temp.Nattens_dB{i}];
end
SORTattens=unique(TempNattens_dB);
NumA=length(SORTattens);
SORTpolarity=unique(Temp.InvertPolarityIND);

%%%%%%%%% ASK WHICH CONDITION TO SHOW  
% Noise Attens
AttenIND=1;
if length(SORTattens)>1
   Atext='';
   for i=1:length(SORTattens)
      Atext=sprintf('%s %d dB [%d],',Atext,SORTattens(i),i);
   end
   Atext=Atext(2:end-1);
   AttenIND=-1;
   while isempty(find(find(SORTattens)==AttenIND))
      AttenIND=input(sprintf('Pick noise atten:  %s: ',Atext));
   end
end

% Invert Polarity
InvertPolarityIND=SORTpolarity(1);
if length(SORTpolarity)>1
   Ptext='';
   for i=1:length(SORTpolarity)
      Ptext=sprintf('%s [1], %s [2]',InvertPolarityText{1},InvertPolarityText{2});
   end
   InvertPolarityIND=-1;
   while isempty(find((1:length(SORTpolarity))==InvertPolarityIND))
      InvertPolarityIND=input(sprintf('Pick Polarity Inversion:   %s: ',Ptext));
   end
end

disp(sprintf('   *** Analyzing Condition: @ %.f dB, Polarity Inversion: %s',SORTattens(AttenIND),upper(InvertPolarityText{InvertPolarityIND})))

% Find relevant picture indices for each condition (These are Masks [1: yes, 0:no] for each picture in picNUMs)
clear TempPolINDs TempINDs
TempPolINDs=(Temp.InvertPolarityIND==InvertPolarityIND);   % 1 for pictures with current: Feature, Harm, Polarity
TempINDs=find(TempPolINDs);

if ~isempty(TempINDs)
   
   % Verify unique conditions
   CONDindFULL=[];
   for PICind=TempINDs
      CONDindFULL=[CONDindFULL find((Temp.Nattens_dB{PICind}==SORTattens(AttenIND)))];
   end
   CONDind=unique(CONDindFULL);
   if length(CONDind)~=1
      error('Non-unique CONDind in TEHeBFi, when looking at FreqIND=%d and LevelIND=%d',FreqIND,LevelIND)
   end
   
   % Find relevant picture indices for each Freq/Level condition (1: yes, 0:no)
   clear TempAttINDs
   for i=1:length(picNUMs)
      TempAttINDs(i)=sum(ismember(Temp.Nattens_dB{i},SORTattens(AttenIND)))>0;
   end
   TempINDs2=find(TempPolINDs&TempAttINDs);  % Mask for Polarity,Natten
   
   picNums=Temp.picNums(TempINDs2);
   excludeLines=cell(size(TempINDs2));
   ii=0;
   for i=TempINDs2
      ii=ii+1;
      TOTALconds=length(x{i}.Stimuli.attens);
      GOODlines=CONDind:TOTALconds:x{i}.Stimuli.fully_presented_lines;
      %% Take out any lines beyond badlines
      GOODlines=GOODlines(find(GOODlines<Temp.FirstBADline(i)));
      IGNORElines=setdiff(1:x{i}.Stimuli.fully_presented_lines,GOODlines);
      excludeLines{ii}=IGNORElines;
   end
   
   SACview(picNums,excludeLines,1);
   
else
   error(upper('No picture with these params!'))
end
return; % END SACrlv



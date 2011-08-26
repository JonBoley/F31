function SR_sps=UnitCalc_SR(ExpDate,UnitName)
% File: UnitCalc_SR.m
% Date: 04Jan2005 (M. Heinz) 
% For: NOHR Experiments
%
% Estimates SR for a given unit. Assumes DataList_xxx.m exists for this experiment (in case UNITsdata/unit does not exist).
% Creates unit.Tn.Um.mat file if it does not exist, otherwise loads and adds to it.
% Returns SRestimate to allow saving to DataList.
%
% ***LATER TO DO***: Need a more accurate way to estimate SR than taking all PICs, because most PST pics are at HLs and therefore
% have reduced SRs.  For now, store all data, and take AVG, but later select by hand for each unit which pictures to average.
%

global NOHR_dir NOHR_ExpList FeaturesText


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

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(sprintf('Processing (calc SR) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

% Verify Unit Number is legitimate
UnitPicList=findPics('*',[TrackNum,UnitNum]);
if isempty(UnitPicList)
   error(sprintf('Unit Number: %d.%02d does not exist!!!!!',TrackNum,UnitNum))
end

%%%%%%%% Create unit structure for this unit
% LOAD UNITSdata file, if it exists, otherwise take from DataList
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
%%%%%%%%%%%%%%%%
%% load unit data if it exists,
if ~isempty(dir(fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum))))
   eval(['load ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''''])
else
   disp(['   *** Loading file: "' FileName '" because UNITSdata/unit does not exist yet'])
   eval(['load ' FileName])
   unit=DataList.Units{TrackNum,UnitNum};
   unit.IgnorePicNums=intersect(UnitPicList,DataList.IgnorePicNums);
end


% See if SR already calculated
CALCsr=1;
if isfield(unit,'SR')
   beep
   temp=input(sprintf('--- SR field has already been calculated for this unit, do you want to recalculate (''y''/[''n''])? '));
   if isempty(temp)
      CALCsr=0;
   elseif ~strcmp(upper(temp(1)),'Y')
      CALCsr=0;
   end
end

if CALCsr
   SRpics=[];
   SRests_sps=[];
   SRdata={};

   TCpicList=findPics('tc',[TrackNum,UnitNum]);
   for PICnum=setdiff(setdiff(UnitPicList,TCpicList),unit.IgnorePicNums)  % All pics except TC and any in the IGNORE list
      %       [SR_sps,lineSRs_sps]=calcSRpic(PICnum,[],[],[.8 1]);  % No excludelines for Now
      [SR_sps,lineSRs_sps]=calcSRpic(PICnum);  % No excludelines for Now
      
      if ~isnan(SR_sps)
         SRpics=[SRpics PICnum];
         SRests_sps=[SRests_sps SR_sps];
         SRdata{length(SRdata)+1}=lineSRs_sps;
      end
      
   end

   unit.SR.SR_sps=mean(SRests_sps);
   unit.SR.SRpics=SRpics;
   unit.SR.SRests=SRests_sps;
   unit.SR.SRdata=SRdata;
   
   unit.Info.SR_sps=unit.SR.SR_sps;

   %    SRests_sps
   
   %%%%%%%%%%%% Save unit data only if updated
   disp(['   *** Saving unit file: "' fullfile(ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) '" ...']) 
   %    beep
   %    disp(' SAVING TURNED OFF')
   eval(['save ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''' unit'])
   
end

SR_sps=unit.SR.SR_sps;


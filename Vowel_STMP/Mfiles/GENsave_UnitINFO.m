function GENsave_UnitINFO(TrackNum,UnitNum,skipVERIFYdir,modelPARAMS)
% File: GENsave_UnitINFO(TrackNum,UnitNum,skipVERIFYdir)
% From: STMPsave_UnitINFO(TrackNum,UnitNum,skipVERIFYdir,modelPARAMS)
% - should be able to be used for any type of data/experiment, thus
% GENERAL!
% FROM: makeNOHRdataList(ExpDate,CHOOSEone)
% Jun 5 2007 M. Heinz
%
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
% Builds BASIC info for a UNIT.  ASSUMES you are in the rawdata directory -
% will verify unless skipped!!
%

if ~exist('TrackNum','var')|~exist('UnitNum','var')
   error('No Track and/or Unit specified')
end
if ~exist('skipVERIFYdir','var')
   skipVERIFYdir=0;
end
if nargin>3 %modelPARAMS supplied
    importINFO=1;
else
    importINFO=0;
end

%% Verify you are in desired directory (with RAW data) 
[pathstr,ExpName]=fileparts(pwd);
if ~skipVERIFYdir
   beep
   TEMP = input(sprintf('\n** Verify current directory (should have RAW data): ''%s'':  0:NO; [1]: YES;  ',ExpName));
   if isempty(TEMP)
      TEMP=1;
   end
   if TEMP~=1
      beep
      disp('***** WRONG CURRENT DIRECTORY - change, and re-run *****');
      return;
   end
end

disp(sprintf('\n... STMP: Processing unitINFO for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum));
UNITinfo_dir='UNITinfo';
if ~exist('UNITinfo','dir')
   mkdir('UNITinfo');
end

unitINFO_filename=sprintf('unitINFO.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile(UNITinfo_dir,unitINFO_filename) ''');'])
if ~isempty(ddd)
   beep
   TEMP = input(sprintf('File: ''%s'' already exists!!\n  ***** Do you want to re-run unitINFO, or leave it as is?  [0]: LEAVE AS IS; RERUN: 1;  ',unitINFO_filename));
   if isempty(TEMP)
      TEMP=0;
   end
   if TEMP~=1
      beep
      disp(sprintf(' FILE NOT ALTERED\n'));
      return;
   else
      disp(' ... Re-Running unitINFO - SAVING NEW FILE!');
   end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General Unit Info (postAnal)
%%%%%%%% SAVED on re-load %%%%%%%%%%%%%%

if importINFO
    unitINFO.Info = modelPARAMS.unit.Info;
else
    unitINFO.Info.BF_kHz=[];  % For now, replace later with TC data
    unitINFO.Info.SR_sps=[];
end

unitINFO.Info.ExpName=ExpName;
unitINFO.Info.Unit=sprintf('%d.%02d',TrackNum,UnitNum);
unitINFO.Info.Threshold_dBSPL=[];  % For now, replace later with TC data
unitINFO.Info.Q10=[];
unitINFO.Info.Comment=[];
unitINFO=calib_parse(unitINFO);
unitINFO.Info.TCindToUse=NaN;
unitINFO.Info.QUALITY=[];

%%%%%%%%%%%%%
%%% Tuning-Curve Info
%%%%%%%%%%%%%
unitINFO=tc_parse(unitINFO);

if ~importINFO
    unitINFO=storeBADlines(unitINFO);  %% Need to do before SR in case of IgnorePICS
    unitINFO=SR_parse(unitINFO)
end

disp(['   *** Saving new file: "' unitINFO_filename '" ...'])
eval(['save ' fullfile(UNITinfo_dir,unitINFO_filename) ' unitINFO'])

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unitINFO=calib_parse(unitINFO) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

picList=findPics('calib');
%%%%%%%%%%%%% Pick which CALIB pic is to be used with this unit
if length(picList)>1
   unitINFO.Info.CalibPIC_ToUse=[];
   while isempty(unitINFO.Info.CalibPIC_ToUse)
      picToUse=input(sprintf('Pick which CALIB picture to use %s: ',mat2str(picList)));
      INDToUse=find(picList==picToUse);
      if ~isempty(INDToUse)
         unitINFO.Info.CalibPIC_ToUse=picList(INDToUse);
      end
   end
elseif length(picList)==1
   unitINFO.Info.CalibPIC_ToUse=picList;
else
   error('NO CALIB FILE FOUND')
end

return;  % CALIB


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unitINFO=tc_parse(unitINFO) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TrackNum=str2num(unitINFO.Info.Unit(1));
UnitNum=str2num(unitINFO.Info.Unit(3:4));

picList=findPics('tc',[TrackNum UnitNum]);

% ONLY get basics automatically here 
[Thresh_dBSPL_ret,BF_kHz_ret,Q10_ret] = plotTCs(picList,unitINFO.Info.CalibPIC_ToUse);
beep
TEMP=input('Does this TC look OK (ow/mark for later analysis)?:  0:NO; [1]: YES;');
if isempty(TEMP)
   TEMP=1;
end
if TEMP==0
   unitINFO.Info.TC_reanalyze_later=1;
end

unitINFO.TC=cell(1,length(picList));
for ind=1:length(picList)
   unitINFO.TC{ind}.picNUM=picList(ind);
   
   unitINFO.TC{ind}.BF_kHz=BF_kHz_ret(ind);
   unitINFO.TC{ind}.Thr_dBSPL=Thresh_dBSPL_ret(ind);
   unitINFO.TC{ind}.Q10=Q10_ret(ind);
end         

% If multiple TCs found, ask which to use
if length(picList)>1
   TCpicToUse=-999;
   while isempty(find(picList==TCpicToUse))
      TCpicToUse=input(sprintf('Pick which TC picture to use %s: ',mat2str(picList)));
   end
   unitINFO.Info.TCindToUse=find(picList==TCpicToUse);
elseif length(picList)==1
   unitINFO.Info.TCindToUse=1;
end

unitINFO.Info.BF_kHz= ...
   unitINFO.TC{unitINFO.Info.TCindToUse}.BF_kHz;  % Use TCindToUse data
unitINFO.Info.Threshold_dBSPL= ...
   unitINFO.TC{unitINFO.Info.TCindToUse}.Thr_dBSPL;  % Use TCindToUse data
unitINFO.Info.Q10= ...
   unitINFO.TC{unitINFO.Info.TCindToUse}.Q10;  % Use TCindToUse data

return;  % TC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unitINFO=storeBADlines(unitINFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM: storeBADlines.m
% Date: 26Apr2006 (M. Heinz)
% For: R03 Experiments
%
% Goes through all PICS for a unit and stores lists of BADlines for each
% PICture.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CHECK FOR: DataList.picBADlines  SKIP IF ALREADY DONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataListDir=dir('DataList*');
if ~isempty(DataListDir)
   eval(['load ' DataListDir.name])
   if isfield(DataList,'picBADlines')
      beep
      disp(sprintf('\n***** BADlines already stored in DataList - WILL USE THOSE *****\n'))
   end
else
   DataList=[];
end

unitINFO.PICS.IgnorePicNums=[];

TrackNum=str2num(unitINFO.Info.Unit(1));
UnitNum=str2num(unitINFO.Info.Unit(3:4));

PIClist = findPics('*',[TrackNum,UnitNum]);
unitINFO.PICS.PIClist = PIClist;

for PICind=PIClist
   if isfield(DataList,'picBADlines')   % Use Stored values - if available
      unitINFO.PICS.picBADlines{find(PIClist==PICind)} = DataList.picBADlines{PICind};
   else
      if ~sum(strcmp(getTAG(getFileName(PICind)),{'calib','tc'}))  % only process RASTER data files

         x=loadPic(PICind);
         disp(sprintf('Picture #: %d, filename: %s',PICind,getfileName(PICind)))

         if isfield(x.General,'trigger')
            disp(sprintf('   Trigger: %s',upper(x.General.trigger)))
            if sum(strcmp(deblank(x.General.trigger),{'Poor','Fair'}))
               beep
            end
         end

         if isfield(x.General,'comment')
            if ~isempty(x.General.comment)
               beep
            end
            disp(sprintf('   Comment: %s',upper(x.General.comment)))
         end

         if isfield(x.General,'run_errors')
            for i=1:length(x.General.run_errors)
               if ~sum(strcmp(x.General.run_errors{i}, ...
                     {'In function ''DALinloop_NI_wavfiles'': Input waveform ', ...
                     'has been truncated to fit requested duration. ', ...
                     'has been repeated to fill requested duration. '}))
                  beep
                  disp(sprintf('   Run_errors: %s',x.General.run_errors{i}))
               end
            end
         end

         %%% Need to look at RASTER HERE
         excludeLines=[];
         PICview(PICind,excludeLines,100)  % looks at all lines in FIG 100
         reviewON = 1;
         %          MAXlines = max(x.spikes{1}(:,1));    %%%%% CHANGED
         %          6/7/07
         MAXlines = x.Stimuli.fully_presented_stimuli;     
         
         while reviewON
            excludeLines
            temp = input(sprintf('Enter vector of lines to exclude (MAXlines = %d)\n(Press Enter to accept current exclusions; -999 to IGNORE PIC): ',MAXlines));
            if isempty(temp)
               reviewON = 0;
            else
               if temp ==-999
                  unitINFO.PICS.IgnorePicNums=[unitINFO.PICS.IgnorePicNums PICind];
                  excludeLines=1:MAXlines;
                  reviewON = 0;
               elseif sum(temp>MAXlines) | sum(temp<1) | ~isnumeric(temp)
                  beep
                  disp('NEED TO RE-ENTER vector!!')
               else
                  excludeLines=temp;
                  PICview(PICind,excludeLines,101)  % Look at PICview with proposed excluded line
               end
            end
         end
         unitINFO.PICS.picBADlines{find(PIClist==PICind)} = excludeLines;

      else
         disp(sprintf('\n**  Skipping PIC = %d (%s) because it is not a raster data file\n',PICind,getFileName(PICind)))
      end

   end
end

return;   % storeBADlines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unitINFO=SR_parse(unitINFO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FROM: UnitCalc_SR.m
%
% Estimates SR for the unit.
% ***LATER TO DO***: Need a more accurate way to estimate SR than taking all PICs, because most PST pics are at HLs and therefore
% have reduced SRs.  For now, store all data, and take AVG, but later select by hand for each unit which pictures to average.
%

TrackNum=str2num(unitINFO.Info.Unit(1));
UnitNum=str2num(unitINFO.Info.Unit(3:4));

SRpics=[];
SRests_sps=[];
SRdata={};

UnitPicList=findPics('*',[TrackNum,UnitNum]);
TCpicList=findPics('tc',[TrackNum,UnitNum]);
for PICnum=setdiff(setdiff(UnitPicList,TCpicList),unitINFO.PICS.IgnorePicNums)  % All pics except TC and any in the IGNORE list
   [SR_sps,lineSRs_sps]=calcSRpic(PICnum);  % No excludelines for Now

   if ~isnan(SR_sps)
      SRpics=[SRpics PICnum];
      SRests_sps=[SRests_sps SR_sps];
      SRdata{length(SRdata)+1}=lineSRs_sps;
   end
end

unitINFO.SRdata.SR_sps=mean(SRests_sps);
unitINFO.SRdata.SRpics=SRpics;
unitINFO.SRdata.SRests=SRests_sps;
unitINFO.SRdata.SRdata=SRdata;

unitINFO.Info.SR_sps=unitINFO.SRdata.SR_sps;

return;  %SR_parse




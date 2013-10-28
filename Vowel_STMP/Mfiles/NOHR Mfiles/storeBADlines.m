function storeBADlines(ExpDate,CHOOSEone)
% File: storeBADlines.m
% Date: 26Apr2006 (M. Heinz)
% For: R03 Experiments
%
% Goes through all PICS in an experiment (after makeNOHRdataList is run)
% and stores lists of BADlines for each PICture.
%
% ExpDate: e.g., '070804' (converted later)
%
% Modified From: makeNOHRdataList.m 
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
%
%

global NOHR_dir NOHR_ExpList NOHR_IMPExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

% Globals for functions
global TUlist UNITind DataList OLD_DataList

if ~exist('ExpDate','var')
   ExpDate=0;
   ExpDate='041805';
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

disp(['Processing Experiment: ''' ExpName ''' --- Storing BADlines ... '])

%%%%%% Verify that DataList File already exist 
if length(dir(FileName))
   disp(['   *** Loading DataList file: "' FileName ''])
   eval(['load ' FileName])
else
	error(['   *** DataList file: "' FileName '" DOES NOT EXIST YET - run makeNOHRdatalist'])
end

% Setup cell array to store BADlines for each picture
% Save old copy if exists to be careful
if isfield(DataList,'picBADlines')
	DataList.OLDpicBADlines = DataList.picBADlines;
else
	DataList.picBADlines = cell(1,DataList.General.TOTALpics);
end



%%%%%%%%%%%%%%%%%%
% WORKING HERE 
% - get loop to process each pic, ignore ones we don't need, show raster and
% error otherwise
% - Store DataList.BADpiclist
% - Resave DataList
%%%%%%%%%%%%%%%%%%%

if CHOOSEone
	Track=input('Enter Track (e.g., 1): ');
	Unit=input('Enter Unit (e.g., 5): ');
	PIClist = findPics('*',[Track,Unit]);
else
	PIClist = 1:DataList.General.TOTALpics;
end


for PICind=PIClist
	if isempty(find(DataList.IgnorePicNums==PICind))  % don't bother with IgnorePICS
		if ~sum(strcmp(getTAG(getFileName(PICind)),...
                {'calib','tc','500','1000','2000','4000','8000','dpoae'}))  % only process RASTER data files

			x=loadPic(PICind);
			disp(sprintf('Picture #: %d, filename: %s',PICind,getFileName(PICind)))

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
         MAXlines = x.Stimuli.fully_presented_stimuli
			
			while reviewON
				excludeLines
				temp = input(sprintf('Enter vector of lines to exclude (MAXlines = %d)\n(Press Enter to accept current exclusions): ',MAXlines));
				if isempty(temp)
					reviewON = 0;
				else
					if sum(temp>MAXlines) | sum(temp<1) | ~isnumeric(temp)
						beep
						disp('NEED TO RE-ENTER vector!!')
					else
						excludeLines=temp;
						PICview(PICind,excludeLines,101)  % Look at PICview with proposed excluded line
					end
				end
			end
			DataList.picBADlines{PICind} = excludeLines;
			
			%%%% NEED TO SAVE DATALIST after every 10th picture
			if round(PICind/10)*10 == PICind
				disp(['   *** Saving new file: "' FileName '" ...'])
				eval(['save ' FileName ' DataList'])
			end
			
		else
			disp(sprintf('\n**  Skipping PIC = %d (%s) because it is not a raster data file\n',PICind,getFileName(PICind)))
		end
	else
		disp(sprintf('\n*****   Skipping PIC = %d because it is in IgnorePIClist\n',PICind))
	end
end

%%%% NEED TO SAVE DATALIST at end
disp(['   *** Saving new file: "' FileName '" ...'])
eval(['save ' FileName ' DataList'])

return;

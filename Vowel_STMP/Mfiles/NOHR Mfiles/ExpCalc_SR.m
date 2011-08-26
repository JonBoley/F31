function ExpCalc_SR(ExpDate)
% File: ExpCalc_SR.m
% Date: 04Jan2005 (M. Heinz)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Modified From: makeNOHRdataList.m 
%
% LATER: Need to get more accurate than AVGing across all PICs, bc HL PSTs have reduced SR re: LL RLVs



global NOHR_dir NOHR_ExpList NOHR_IMPExpList

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

disp(['Processing Experiment (SR calcs): ''' ExpName ''' ... '])

%%%%%% Load DataList file to get list of units, and save at end
disp(['   *** Loading old file: "' FileName '";'])
eval(['load ' FileName])


%%%%%%%%%%%%%%%% Get Track/Unit List for Experiment
TUlist=getTrackUnitList;

for UNITind=1:size(TUlist,1)
   % for UNITind=36
   %    beep
   %    disp(sprintf('HARD CODED: UNITind: %d',UNITind))
   UnitName=sprintf('%d.%02d',TUlist(UNITind,:));
   disp(sprintf('   Unit: %s',UnitName))

   SR_sps=UnitCalc_SR(ExpDate,UnitName);   

   % update DataList
   DataList.Units{TUlist(UNITind,1),TUlist(UNITind,2)}.Info.SR_sps=SR_sps;

end

%%%%%%%%%%%% Save DataList at end
disp(['   *** RE-Saving updated file: "' FileName '" ...'])
eval(['save ' FileName ' DataList'])

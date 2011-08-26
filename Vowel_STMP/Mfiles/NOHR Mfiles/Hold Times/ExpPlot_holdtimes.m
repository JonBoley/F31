function ExpPlot_holdtimes(ExpDate)
% File: ExpPlot_holdtimes.m
% Date: 08Nov2004 (M. Heinz) (Modified from ExpPlot_EHFF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Plots distribution of holding times for all units in an Experiment
%

global NOHR_dir NOHR_ExpList

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='090204'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end

%%%% Find the full Experiment Name 
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
disp(sprintf('Plotting  (holdtimes) Experiment: ''%s''',ExpName))

eval(['cd ''' fullfile(data_dir,ExpName) ''''])

TrackUnitList=getTrackUnitList;
holdtimes_min=zeros(1,size(TrackUnitList,1))+NaN;

for UNITind=1:length(holdtimes_min)
   
   holdtimes_min(UNITind)=getUnitHoldTime(TrackUnitList(UNITind,1),TrackUnitList(UNITind,2));

   if holdtimes_min(UNITind) > 45
      disp(sprintf('Unit: %d.%02d has BIG holdtime = %.2f min',TrackUnitList(UNITind,1),TrackUnitList(UNITind,2),holdtimes_min(UNITind)))
   elseif isnan(holdtimes_min(UNITind))
      disp(sprintf('Unit: %d.%02d has holdtime set to NaN',TrackUnitList(UNITind,1),TrackUnitList(UNITind,2)))
   end
   
end


disp('*** Analysis ignores only-TC units ***')
MEANholdtime_min=mean(holdtimes_min(find((~(isnan(holdtimes_min)|(holdtimes_min==0))))));
MEDIANholdtime_min=median(holdtimes_min(find((~(isnan(holdtimes_min)|(holdtimes_min==0))))));

figure(1); clf

hist(holdtimes_min(find((~(isnan(holdtimes_min)|(holdtimes_min==0))))))
title(sprintf('Distribution of Hold Times for Experiment: ''%s''',ExpName),'Interpreter','none')
xlabel('Hold Time (min)')
ylabel('Number of Units')
text(.7,.9,sprintf('MEAN = %.1f min',MEANholdtime_min),'units','norm')
text(.7,.85,sprintf('MEDIAN = %.1f min',MEDIANholdtime_min),'units','norm')

return;
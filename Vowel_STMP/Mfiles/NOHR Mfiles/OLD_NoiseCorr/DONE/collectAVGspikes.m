function [Spikes,levels]=collectAVGspikes(unit,STIMind,paramsIN);
% File: collectAVGspikes.m
%
% 7/26/03: collects all spikes for a given RLF (stimulus and unit) from AVGdata in a cell array for further analysis
% No picking (except Start/End times) or pooling of spikes, all spikes and REPs returned that were used in AVGdata{STIMind}
%
% From: PERhist.m
%
% Spikes: cell array{NUMreps,NUMlevels} with ALL spikes (within Start/End-Time) at a given level
% levels: vector with levels in dB SPL (from AVGdata)

global Recruit_dir

ExpName=strcat('Exp',unit.ExpDate);
OLDdir=pwd;
eval(['cd ''' fullfile(Recruit_dir,'NewData',ExpName,'rawdata') ''''])

%%%%%%%%% Take out parameters if passed
if exist('paramsIN','var')
   if isfield(paramsIN,'StartTime')
      StartTime=paramsIN.StartTime;
   else
      StartTime=-Inf;
   end
   if isfield(paramsIN,'EndTime')
      EndTime=paramsIN.EndTime;
   else
      EndTime=Inf;
   end
else
   StartTime=-Inf;
   EndTime=Inf;
end

%%%%%%%%%%%% Find pictNUMS used in AVG
pictNUMS=unit.REPdata.pictNUMS(find(unit.AVGdata.AVGinds{STIMind}),STIMind);
NUMreps=length(pictNUMS);

%%%%%%%%% Find levels in AVG rate-level (dB SPL)
levels=unit.AVGdata.levels{STIMind};
NUMlevels=length(levels);



%%%%%%%%%%%% Get raw spikes (includes ALL abvals recorded, so need to check later which levels used in REPdata!)
for REPind=1:NUMreps
   if(unit.system=='PDP')
      RLfile=strcat('a',num2str(pictNUMS(REPind)));
      [spikesREP{REPind},abvalsREP{REPind}]=getspikes_PDP(RLfile);
   else
      d=dir(sprintf('*p%04d*',pictNUMS(REPind)));
      for Dind=1:length(d)
         if strcmp(d(Dind).name(1),'p')
            RLfile=d(Dind).name(1:findstr('.m',d(Dind).name)-1);
         end
      end
      [spikesREP{REPind},abvalsREP{REPind}]=getspikes_TDT(RLfile);
   end
end

%%%%%%%%%% Find abvals for each level
cal_level=unit.Calib.cal_level{STIMind};
abvals=levels-cal_level;

%%%%%%%%%% Put all spikes together for each level
Spikes=cell(NUMreps,NUMlevels); % ALL spikes at a given level (ALL times within Start/End-Time)
for LEVind=1:NUMlevels  %Only look at levels in AVG-RLF
   for REPind=1:NUMreps
      %%% VERIFY THIS LEVEL WAS INCLUDED IN REPdata (i.e., after cleanup)!!!
      if(find(unit.REPdata.levels{find(unit.REPdata.pictNUMS(:,STIMind)==pictNUMS(REPind)),STIMind}==levels(LEVind)))
         seqIND=find(abvalsREP{REPind}==abvals(LEVind));
         Spikes{REPind,LEVind}=spikesREP{REPind}(find((spikesREP{REPind}(:,1)==seqIND)& ...
            (spikesREP{REPind}(:,2)>=StartTime)&(spikesREP{REPind}(:,2)<=EndTime)),2)';
      end
   end
end

eval(['cd ''' OLDdir ''''])
return;


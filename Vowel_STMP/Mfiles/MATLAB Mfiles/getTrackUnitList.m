function TrackUnitList=getTrackUnitList
% function TrackUnitList=getTrackUnitList
% Created: M. Heinz 19Mar2004
% For: CNexps (GE/MH)
%
% Returns List of all Units for an Experiment: TrackUnitList=[TrackList,UnitList]

d=dir('p*');
TrackUnitList=NaN+zeros(length(d),2);

for IND=1:length(d)
   if ~d(IND).isdir
      TrackUnitList(IND,:)=getTrackUnit(d(IND).name);
   end   
end

TrackUnitList=unique(TrackUnitList(~isnan(TrackUnitList(:,1)),:),'rows');
return;

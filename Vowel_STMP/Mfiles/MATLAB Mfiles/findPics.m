function picList=findPics(STRING,varargin)
% function picList=findPics(STRING,varargin)
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
% Returns picture-number list from current directory of all files with STRING in their name
% varargin: [can be in any order]
%         : 'notTAG': indicates to perform a general search for STRING, i.e., not required to be a TAG
%         : TrackUnitName=[Track, Unit]

ADDstring1='_'; % Adds to the beginning of STRING to make it a TAG
ADDstring2='.'; % Adds to the end of STRING to make it a TAG

for i=1:length(varargin)
   if ischar(varargin{i})
      if strcmp(varargin{i},'notTAG')
         ADDstring1='';  % Remove TAG requirement
         ADDstring2='';
      else
         error(sprintf('varargin(%d)="%s" is an unidentified string',i,varargin{i}))
      end
   elseif isnumeric(varargin{i})
      if length(varargin{i})==2
         TrackUnitName=varargin{i};
      else
         error(sprintf('Numeric varargin(%d) is length=%d, MUST BE length=2 [Track,Unit]',i,length(varargin{i})))
      end
   else
      error(sprintf('varargin(%d) is neither char no numeric'))
   end
end      
   
picList=[];

d=dir(strcat('*',ADDstring1,STRING,ADDstring2,'*'));
picList=NaN+zeros(1,length(d));

for IND=1:length(d)
   if ~d(IND).isdir
      if exist('TrackUnitName','var')
         if sum(TrackUnitName==getTrackUnit(d(IND).name))==2
            picList(IND)=getPicNum(d(IND).name);
         end
      else
         picList(IND)=getPicNum(d(IND).name);
      end   
   end
end
picList=picList(~isnan(picList));

return;




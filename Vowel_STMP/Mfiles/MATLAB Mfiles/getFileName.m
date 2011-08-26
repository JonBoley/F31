function filename=getFileName(picIND)
% function x=getFileName(picIND)
% Created: M. Heinz 19Mar2004
% For: CNexps (GE/MH)
%
% Returns filename for given picture number

d=dir(strcat('*',sprintf('p%04d',picIND),'*'));
if length(d)>1
   warning('More than 1 file with this picture number');
   filename='';
elseif isempty(d)
   warning('Picture does not exists');
   filename='';
else
   filename=d.name;
end

return;
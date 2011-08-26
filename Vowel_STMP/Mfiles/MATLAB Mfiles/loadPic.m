function [x,errorMSG]=loadPic(picNum)
% function x=loadPic(picIND)
% Created: M. Heinz 18Mar2004
% Modified from GE version, 29Jul2004
%
% Loads picture based on picture number

picSearchString = sprintf('*p%04d*.m', picNum);
picMFile = dir(picSearchString);
errorMSG = '';
if isempty(picMFile)
   errorMSG = sprintf('Picture file p%04d*.m not found.', picNum);
   x = [];
   return;
elseif length(picMFile)>1
   errorMSG = sprintf('More than 1 file with this picture number (p%04d)', picNum);
   x = [];
   return;
else
   eval(['x = ' picMFile.name(1:end-2) ';']);  % Load data for this picture
end

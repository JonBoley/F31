function TAG=getTAG(filename)
% function TAG=getTAG(filename)
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
% Returns TAG from picture filename

if ~isempty(filename)
   if strcmp(filename(1),'p')
      ULINEinds=findstr(filename,'_');
      PERinds=findstr(filename,'.');
      
      TAG=filename(max(ULINEinds)+1:min(PERinds)-1);
      return;
   end
end

TAG='';

return;

function [a,fid] = fill_marked_val(fname,mark,a)
%

% AF 9/4/01

if (mark(1) == '%')
   fmt = ['%' mark ' ' repmat('%f ',1, size(a,2)) '\n'];
else
   fmt = [mark ' ' repmat('%f ',1, size(a,2)) '\n'];
end  
if (isnumeric(fname))
   fid = fname;
else
   fid = fopen(fname);
end
ll = fgets(fid);
while (~isempty(ll) & (ll > 0))
   if (strncmp(ll,mark,length(mark)))
      a = fscanf(fid,fmt,size(a'))';
      return;
   end
   ll = fgets(fid);
end
fclose(fid);
function browse_struct(s, title,f_units, max_width, max_struct_elem)
%

% AF 11/8/01

struct_name = inputname(1);
if ((exist('title','var') ~= 1) | isempty(title))
   title = 'browse structure';
   if (~isempty(struct_name))
      title = [title ' ''' struct_name ''''];
   end
end
if (exist('f_units','var') ~= 1)
   f_units = struct([]);
end
if (exist('max_width','var') ~= 1)
   max_width  = Inf;
end
if (exist('max_struct_elem','var') ~= 1)
   max_struct_elem = 3;
end

str = struct2str(s,f_units,max_width,max_struct_elem);
if (~isempty(struct_name))
   str = cat(1,{struct_name},str);
end

strdlg(str, title);
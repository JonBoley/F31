function convert_oldnel_name(oldname)

startp = findstr('_p',oldname);
if (~isempty(startp))
   i = startp+2;
   while ((i<=length(oldname)) & ~isnan(str2double(oldname(i))))
      i = i+1;
   end
   if (oldname(i) == '_' | ~isnan(str2double(oldname(i))))
      endp = i;
      extra_ = '';
   else
      endp = i-1;
      extra_ = '_';
  end
  newname = ['p0' oldname(startp+2:endp) extra_ oldname(1:startp) oldname(endp+1:end)];  % Adds extra 0
  dos(['rename ' oldname ' ' newname]);
end

disp(sprintf('%s\n%s\n',oldname,newname))


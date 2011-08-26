function dir_traverse(dirname,cmd,parameter)
% DIR_TRAVERSE executes a command on each file in a directory (recursively).
%              dir_traverse(dirname,cmd). cmd gets one parameter which is a file name 
%              Example: dir_traverse('d:\temp\old-phd','lower_case_fname')

% Alon Fishbach 13/8/99

do_recurse = 1;
verbose = 1;

if (exist('parameter') ~= 1)
   parameter = [];
end
parrent = pwd;
cd(dirname);
fprintf('--------------\n\t\tGoing into %s Directory\n',dirname);
%eval(cmd);
files = dir;
for i = 1:length(files)
   if (files(i).isdir)
      if (do_recurse)
         if (files(i).name(1) ~= '.')
            eval([cmd '( files(i).name )']);
            dir_traverse(files(i).name,cmd,parameter); %% Rrecursive call;
         end
      end
   else
      if (verbose)
         fprintf('%s\n', files(i).name);
      end
      if (isempty(parameter))
         eval([cmd '( files(i).name );']);
      else
         eval([cmd '( files(i).name,parameter );']);
      end
   end
end
cd(parrent)

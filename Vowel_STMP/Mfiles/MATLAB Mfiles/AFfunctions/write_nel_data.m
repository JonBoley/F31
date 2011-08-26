function rc = write_nel_data(fname,x,save_spikes)
% write_nel_data  - writes a standard nel output file in a form of an m-file.
%
%       Usage: rc = write_nel_data(fname,x,save_spikes)
%                   'x' should contain the data to be saved in 'fname'.
%                   'save_spikes' is an optional boolean flag 
%                   (by default save_spikes = 1).
%                   'write_nel_data' returns 1 on success and -1 on failure.

% AF 9/5/01

global spikes

if (exist('save_spikes') ~= 1)
   save_spikes = 1;
end
[dummy1 dummy2 ext] = fileparts(fname);
if (strcmp(ext,'.m') ~= 1)
   fname = [fname '.m'];
end
fid = fopen(fname,'wt+');
if (fid < 0)
   rc = -1;
   return;
end
[dirname filename] = fileparts(fname);
fprintf(fid,'function x = %s\n', filename);
mat2text(x,fid);

if (save_spikes)
   for i = 1:length(spikes.times)
      % Code for initializing the spike matrix
      fprintf(fid, '\nx.spikes{%d} = zeros(%d,%d);\n', i, spikes.last(i), size(spikes.times{i},2));
      % Code for self-extracting the data stored in the file comments.
      if (i==1)
         fprintf(fid, '[x.spikes{%d},fid] = fill_marked_val(which(mfilename),''%%%d'',x.spikes{%d});\n', i,i,i);
      else
         fprintf(fid, '[x.spikes{%d},fid] = fill_marked_val(fid,''%%%d'',x.spikes{%d});\n', i,i,i);
      end
   end
   fprintf(fid, 'fclose(fid);');
   % Writing the spike data as comments to the data mfile.
   for i = 1:length(spikes.times)
      fprintf(fid, '\n%%%d\n',i);
      fmt = ['%%' int2str(i) ' %8d %4.10f \n'];
      fprintf(fid,fmt,spikes.times{i}(1:spikes.last(i),:)');
   end
   fprintf(fid,'\n');
   % Add the source code for the subfunction 'fill_marked_val' to the data file.
   subfunc_file = textread(which('fill_marked_val'),'%s','delimiter','\n','whitespace','');
   fprintf(fid,'%s\n',subfunc_file{:});
end
fclose(fid);
rc = 1;
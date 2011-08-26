function mat2text(x,fid,level)
% MAT2TEXT - writes a matlab variable to a text file in a matlab mfile format.
%            The variable can be retrieved by evaluating (executing) the file.
%            The variable can be of any intrinsic matlab format (struct, cell,
%            ND array, string and sparse).
%            USAGE: mat2text(x,fid).  
%                   fid is an output file descriptor. 
%                   Use fid=1 to direct the output to the screen.

% AF 8/6/01

if (exist('fid') ~= 1)
   fid = 1;
end
if (exist('level') ~= 1)
   level = 0;
end
level = level+1;

if (level == 1)
   fprintf(fid,'x = ');
end

switch (class(x))
case 'struct'
   vals = struct2cell(x);
   names = fieldnames(x);
   fprintf(fid,'struct(');
   for i = 1:length(names)
      if (i>1)
         fprintf(fid,',');
      end
      fprintf(fid,'''%s'', ',names{i});
      mat2text(vals(i,:),fid,level);
   end
   fprintf(fid,')');
   
case 'char'
   % fprintf(fid,'''');
   fprintf(fid,'[');
   for ii = 1:size(x,1) % we replace ' with ''
      tmpstr = strrep(x(ii,:),'''',''''''); %% the loop is due to strrep limitation
      line_limit_fprintf(fid,'''%s ''', deblank(tmpstr),100);
      % fprintf(fid,'%s', deblank(tmpstr));
      if (ii < size(x,1))
         fprintf(fid,' ');
      end
   end
   fprintf(fid,']');
   % fprintf(fid,'''');
   
case 'cell'
   if (ndims(x) > 2)
      fprintf(fid, 'cat(%d, ', ndims(x));
      for i = 1:size(x,ndims(x))
         if (i>1)
            fprintf(fid,',');
         end
         eval(['tmpmat = squeeze(x(' repmat(':,',[1 ndims(x)-1]) int2str(i) '));']);
         mat2text(tmpmat,fid,level);
      end
      fprintf(fid,') ...\n ');
   else
      if (size(x,2) == 1)
         fprintf(fid,'{');
         for i = 1:size(x,1)
            mat2text(x{i},fid,level);
            if (i < size(x,1))
               fprintf(fid,';');
            end
         end
         fprintf(fid,'} ...\n');
      else
         fprintf(fid, '[');
         for i = 1:size(x,2)
            mat2text(x(:,i),fid,level);
         end
         fprintf(fid, '] ...\n');
      end
   end
      
case 'sparse'
   fprintf(fid, 'sparse(');
   mat2text(full(x),fid,level);
   fprintf(fid, ')\n');
   
case 'double'
   if (isempty(x))
      fprintf(fid,'[]');
      return;
   end
   if (ndims(x) > 2)
      fprintf(fid, 'cat(%d, ', ndims(x));
      for i = 1:size(x,ndims(x))
         if (i>1)
            fprintf(fid,',');
         end
         eval(['tmpmat = squeeze(x(' repmat(':,',[1 ndims(x)-1]) int2str(i) '));']);
         mat2text(tmpmat,fid,level);
      end
      fprintf(fid,') ...\n ');
   else
      if (size(x,2) == 1)
         if (size(x,1) == 1)
            fprintf(fid,'%1.10g ',x);
         else
            fprintf(fid,'[');
            line_limit_fprintf(fid,'%1.10g;',x);
            % fprintf(fid,'%1.10g;',x);
            fprintf(fid,'] ...\n');
         end
      else
         if (size(x,1) == 1)
            fprintf(fid,'[');
            line_limit_fprintf(fid,'%1.10g ',x);
            % fprintf(fid,'%1.10g ',x);
            fprintf(fid,'] ...\n');
         else
            fprintf(fid, '[');
            %frmt = ['[' repmat('%1.10g ', 1, size(x,2)) ']; \n'];
            %fprintf(fid, frmt, x');
            for i = 1:size(x,2)
               mat2text(x(:,i),fid,level);
            end
            fprintf(fid, '] ...\n');
         end
      end
   end
end
   
   
if (level == 1)
   fprintf(fid,';\n');
end

function line_limit_fprintf(fid,frmt,x,N)
if (exist('N','var') ~= 1)
   N = 10;
end
for i = 1:ceil(length(x)/N)
   from = (i-1)*N+1;
   to   = min(i*N,length(x));
   fprintf(fid,frmt,x(from:to));
   fprintf(fid,' ...\n');
end

function [fname, nfils, binset, binwset, autofilt, filtw] = parscomm(str);

% [fname,nfils,binset,binwset,autofilt,filtw] = PARSCOMM(str)
% Parse command line string. This can take the form:
%   'a4-7,12,15-17b2f7'
% where
%    a4-7,12,15-17 means process P. Nos. 4-7, 12, 15-17
%    b2 means use a binwidth of 2 ms (optional)
%    f7 means filter with a 7-bin wide triangular filter (optional)
% fname(1:nfils) return the file names as 'a4', 'a5', . . .
% binset~=0 if binwidth specified in binwset
% autofilt~=0 if filtering desired, filter bw in filtw

% fnames holds the filenames to be processed (up to 20)
	fname={'', '', '', '', '', '', '', '', '', '', ...
         	'', '', '', '', '', '', '', '', '', ''};
			

% binset = 1 if binwidth set in call string, otherwise determined from file
% autofilt = 1 if filtering specified in call string.
	binset = 0;
	autofilt = 0;
	binwset = 0.1;
	filtw = 1;
	
% Parse out the elements of the command string
	n1 = 1;
	n2 = size(str,2);
	jfn = 1;
	for j=1:22
	   if n1>=n2
			break
	   else
	      [A, nwds, ermsg, nch] = sscanf(str(n1:n2), '%c%g', 2);
	      n1 = n1 + nch -1;
	      if nwds<2
		      break
	      else
		     if j==1
			   fname(jfn) = strcat(fname(jfn), sprintf('%c%g',A(1),A(2)));
			   jfn = jfn + 1;
			   B = A(1);
			   C = A(2)+1;
		     else
			   if char(A(1))=='-'
				   for j=C:A(2)
					   fname(jfn) = strcat(fname(jfn), sprintf('%c%g',B,j));
				       jfn = jfn + 1;
				   end
			   elseif char(A(1))==','
				   fname(jfn) = strcat(fname(jfn), sprintf('%c%g',B,A(2)));
				   jfn = jfn + 1;
				   C = A(2)+1;
			   elseif char(A(1))=='b' | char(A(1))=='B'
				   binset = 1;
				   binwset = A(2);
			   elseif char(A(1))=='f' | char(A(1))=='F'
				   autofilt = 1;
				   filtw = A(2);
			   end
		     end
	      end
	   end
   end
   nfils = jfn - 1;
return

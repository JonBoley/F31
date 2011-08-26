
function pstm(str)
% PSTM('a4') - computes and displays PST histograms of data in the file 'a4.DAT'
% which should be in the lab standard format.
% PSTM('a4b0.5f7') computes a PST histogram with binwidth 0.5 ms and filters
% it with a triangular window of width 7 bins.
% Filename must be first, but binwidth and filterwidth can be in either order.
% SLOW matlab code.

global PSTDATXYZ NLINFXYZ

	legend = {'', '', '', '', '', '', '', '', '', '', ...
	          '', '', '', '', '', '', '', '', '', ''};


% Parse the input command string:
	[fname, nfils, binset, binwset, autofilt, filtw] = parscmpst(str);
	
	if nfils>0		% IF #1
	   nfilsin = 0;
	   for jf=1:nfils
% Read parameter block, determine file type, make sure files are more or less
% the same (note this doesn't check very much); initialize storage.
	      pb = labfile(char(fname(jf)));
	      if strcmp(pb.Error, 'No error')~=1     % IF #2
		     fprintf(1,'\n*** Error = ''%s'' opening %s; check the directory.***\n', ...
			      pb.Error,char(fname))
		     fclose('all');
	      else

	         [ftype, binw, nbins, npsts] = ppblpst(pb, binset, binwset);
			 skipit = 0;
			 if nfilsin==0
				 PSTDATXYZ = zeros(nbins, npsts);
				 NLINFXYZ = zeros(1,npsts);
				 compare = [ftype, binw, nbins, npsts];
		     else
				 if compare ~= [ftype, binw, nbins, npsts];
					 fprintf(1,'\n***File %s is different, skipping it.***\n', ...
						char(fname(jf)))
					 skipit = 1;
				 end
			 end

% Make PST histogram(s) and legend text
		     if skipit==0
	            comppst(pb.resol, ftype, binw, nbins, npsts);
				nfilsin = nfilsin + 1;
				legend(nfilsin) = strcat(legend(nfilsin), ...
                    sprintf('File %s, %s.', char(fname(jf)), ...
				 	pb.ltext(1:lens(pb.ltext))));
			 end
	      end			% ENDIF #2
	   end				% ENDFOR jf=1:nfils
	  

% Compute display title
	   dtit = sprintf('Bin width = %g ms.',binw);
	   
% Convert PST counts to spikes/s
%	   PSTDATXYZ = PSTDATXYZ/diag(NLINFXYZ*binw/1000.);
	   for j=1:npsts
		   if NLINFXYZ(j)>0
		  PSTDATXYZ(:,j) =PSTDATXYZ(:,j)/(NLINFXYZ(j)*binw/1000.);
		   end
	   end

% Autofilter
	   if autofilt==1
	      for j=1,npsts
	         PSTDATXYZ(:,j) = trifilt(PSTDATXYZ(:,j)', filtw)';
	      end
	   end

% Compute abscissa
%	   absc=zeros(nbins,1);
%	   absc(1) = binw/2;
%	   for j=2:nbins
%		  absc(j) = absc(j-1) + binw;
%	   end
	   absc = binw*(tril(ones(nbins,nbins))*ones(nbins,1)-0.5*ones(nbins,1));

% If there are multiple PSTs, offset them for display
	   if npsts>1
		   mxrate = max(max(PSTDATXYZ));
		   addto = mxrate/2*ones(nbins,1)* ...
		        (ones(1,npsts)*triu(ones(npsts,npsts),1));
		   PSTDATXYZ = PSTDATXYZ + addto;
	   end
   
% Display PSTs. NOTE THAT THE FIRST BIN IS NOT PLOTTED. This is because of
% the artifactual triggers at time 0.  
subplot(311);
	   plot(absc(2:nbins), PSTDATXYZ(2:nbins,:),'b-')
	   xlabel('Milliseconds')
       ylabel('Rate, spikes/s')
	   title(dtit)

	   axlim = axis;
	   xleg = axlim(1)+0.5*(axlim(2)-axlim(1));
       yleg = axlim(3) + 0.95*(axlim(4)-axlim(3));
	   ydlt = 0.05*(axlim(4)-axlim(3));
	   for j=1:nfilsin
		  text(xleg, yleg, legend(j));
		  yleg = yleg - ydlt;
	   end
 
		

	end			% ENDIF #1 
return			% from function






function [ftype, binw, nbins, npsts] = ppblpst(pb, binset, binwset);

% Parse parameter block in pb for PST histogram computation. Returns:
% ftype = file type = 0 - normal PST, constant stimulus over all lines
%                   = 1 - sweep stimulus, with large number of different
%						  stimuli (i.e. freq. sweep, rate-level, etc.)
%					= 2 - file rotation stimulus, from STR with smaller
%						  number of different stimuli.
% binw - recommended PST binwidth.
% nbins - recommended number of bins.
%     if binset = 1, binw = binwset (at least 0.1 ms) and nbins is
%                    min(2000, 2*dur/binw)
%               = 0, binw = dur/250 and nbins = 500.
% npsts - number of psts in the file. There is one PST for each different
%    stimulus. NOTE that if ftype=0, then npsts successive lines should be
%    combined into each PST, whereas if ftype=2, then every npsts'th line
%    should be combined.

% Get stimulus duration and period for number of bins computation
	perd = pb.mode_line(4);
	dur = 0;
	if pb.stim_par(1) ~= 11
	   dur = pb.stim_par(9);
	end
	if pb.stim_par(3) ~= 11
		dur = max([dur, pb.stim_par(13)]);
	end

% Binwidth is externally set or is dur/250.
% number of bins is sufficient to make PST cover 2*dur, limited to 4000 bins
	if binset~=0
		binw = max([binwset, 0.01]);
	else
		binw = round(10*dur/250)/10;
	end
	nbins = round(min([4000, 2*dur/binw, (perd-2)/binw]));

% Compute file type and number of psts in the file
	if pb.mode_line(7)==0
		ftype = 0;
		npsts = 1;
	elseif pb.mode_line(7)~=5
		ftype = 1;
		npsts= fix(100/pb.mode_line(1));
	else
		ftype = 2;
		npsts = fix(100/pb.mode_line(1));
	end
return







function comppst(resol, ftype, binw, nbins, npsts)

% Compute pst histograms for file open on FILEPNTRXYZ. Histogram has binwidth
% binw ms, data have time resolution resol microsec. Each histogram has nbins
% bins with latencies 0 through nbins*binw. npsts histograms are made.
% If ftype=1,
%     lines are combined in chunks of 100/npsts (1-5, 6-10, etc.)
% if ftype=2,
%     every npsts'th line is combined (1,6,11, ...; 2, 7, 12, ...; etc.)
% Each PST histogram is a column vector of pstdat.
% nlinf(j) is row vector of the number of lines of data included in the
% jth pst, i.e. in pstdat(:,j). Usually, nlinf = 100/npsts.


global FILEPNTRXYZ PSTDATXYZ NLINFXYZ


% Binwidth and max latency in a PST, in clock cycles
	binwcl = 1000*binw/resol;
	pstleng = binwcl*nbins;

	nlinpst = fix(100/npsts);
	jpst = 0;
	in = 0;
	while in ~= -16384
		in = fread(FILEPNTRXYZ, 1, 'int16');
		if in == -32768
			cum = cum + 32768;
		elseif in < 0
			in1 = fread(FILEPNTRXYZ, 1, 'int16');
			if jpst>0
			   lineleng = cum + (in1 + 32768);
			   if lineleng >= pstleng
				   PSTDATXYZ(:,jpst) = PSTDATXYZ(:,jpst) + thpst;
				   NLINFXYZ(jpst) = NLINFXYZ(jpst) + 1;
			   end
			end
			thpst = zeros(nbins,1);
		   	if in~=-16384
			   line = in + 32768;
			   if ftype~=2
				   jpst = 1 + fix((line-1)/nlinpst);
			   else
			       jpst = 1 + rem((line-1), npsts);
			   end
			end
			cum = 0.;
		else
			if jpst <= 0 
				fprintf(1,'***ERROR spikes before a line-start.***')
				in = -16384;
			else
				lat = cum + in;
				index = 1 + fix(lat/binwcl);
				if index>0 & index<=nbins
				   thpst(index) = thpst(index) + 1;
				end
			end
		end
	end
	fclose(FILEPNTRXYZ);
return



function [fname, nfils, binset, binwset, autofilt, filtw] = parscmpst(str);

% [fname,nfils,binset,binwset,autofilt,filtw] = PARSCOMM(str)
% Parse command line string. This can take the form:
%   'a4-7,12,15-17b2f7'
% where
%    a4-7,12,15-17 means process P. Nos. 4-7, 12, 15-17
%    b2 means use a binwidth of 2 ms (optional)
%    f7 means filter with a 7-bin wide triangular filter (optional)
% fname(1:nfils) return the file names as 'a4', 'a5', . . .
% binset~=0 if binwidth specified in binwset
% autofilt~=0 if filtering desired, filter bw in filtw. NOTE filtw is returned
% odd, so if 4 is entered, 5 is returned and used.

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
   if autofilt~=0
	   filtw = max([1, 1+2*fix(filtw/2)])
   end
return




function nch = lens(str)

% Returns the length of the string, up to the last printable character

	for j=size(str,2):-1:1
		nch = j;
		if abs(str(j))>abs(' ')
			break
		end
	end
return

	function pst(str);
% PST('a4-6,9') - computes and displays PST histograms of data in the files
% a4.DAT, a5.DAT, a6.DAT, a9.DAT which should be in the lab standard format.
% PST('a4b0.5f7') computes a PST histogram with binwidth 0.5 ms and filters
% it with a triangular window of width 7 bins.
% Filename must be first, but binwidth and filterwidth can be in either order.
% Provides routines to handle Hafter TDT click data.

% NOTE THIS PROGRAM ASSUMES THAT TDT CLICK TRAIN DATA ARE TAKEN WITH 4 REPS
% OF 25 STIMULI. IF NOT TAKEN THIS WAY, THEN ANALYSIS BY SIGMA GROUP FAILS.

global PSTDATXYZ NLINFXYZ

	legend = {'', '', '', '', '', '', '', '', '', '', ...
	          '', '', '', '', '', '', '', '', '', ''};

% DATA STRUCTURES containing layout of click trains for Hafter expt.
% sigmav are row vectors of sigma values for various stimulus sets:
% QV row 1; RW row 2; . . . UZ row5.
	sigmav = [0.1,    0.215, 0.464, 1, 2.15
	          0.1,   0.0464, 0.464, 1, 2.15
			  0.0215,0.0464, 0.215, 1, 2.15
			  0.01,  0.0464, 0.215, 1, 2.15
			  0.01,  0.0464, 0.215, 1, 2.15];
			  
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

	         [ftype, binw, nbins, nbnstim, npsts] = ppblpst(pb, binset, binwset);
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
	   dtit = strcat(pb.date, ', ', pb.unit, ...
	            sprintf(', bin width = %g ms.',binw));
	   
% Convert PST counts to spikes/s for normal PSTs, but to spikes/stimulus
% for TDT click trains
	   if ftype~=2
	      for j=1:npsts
		      if NLINFXYZ(j)>0
		  PSTDATXYZ(:,j) =PSTDATXYZ(:,j)/(NLINFXYZ(j)*binw/1000.);
	          end
	      end
	   else
		  for j=1:npsts
		      if NLINFXYZ(j)>0
		  PSTDATXYZ(:,j) =PSTDATXYZ(:,j)/NLINFXYZ(j);
	          end
		  end	  
	   end
  
% Set bin 1 of all PSTs to zero to eliminate the 0-triggering artifact
	   PSTDATXYZ(1,:) = zeros(1,npsts);

% Autofilter
	   if autofilt==1
	      for j=1:npsts
	         PSTDATXYZ(:,j) = trifilt(PSTDATXYZ(:,j)', filtw)';
	      end
	   end

% Compute abscissa
	   absc = binw*([1:nbins]'-0.5*ones(nbins,1));
%	   absc=zeros(nbins,1);
%	   absc(1) = binw/2;
%	   for j=2:nbins
%		  absc(j) = absc(j-1) + binw;
%	   end
%	   absc = binw*(tril(ones(nbins,nbins))*ones(nbins,1)-0.5*ones(nbins,1));

	   showndex = 7;
	   while showndex~=6
		   
% If there are multiple PSTs, offset them for display. Set the first bin in
% the top multiple-PST to the max rate, to control scaling.
	       if npsts>1
			   mxrate = max([0.5, max(max(PSTDATXYZ))]);
			   addto = (3*mxrate/4)*ones(nbins,1)*[0:npsts-1];
%			   addto = (3*mxrate/4)*ones(nbins,1)* ...
%			        (ones(1,npsts)*triu(ones(npsts,npsts),1));
			   displdat = PSTDATXYZ + addto;
			   displdat(1,npsts) = mxrate + addto(1,npsts);
		   else
			   displdat = PSTDATXYZ;
		   end
   
% Display PSTs. NOTE THAT THE FIRST BIN WAS SET TO ZERO ABOVE.
		   plot(absc, displdat,'b-')
		   xlabel('Milliseconds')
		   if ftype==2
			   ylabel('Number of spikes/stimulus')
		   else
	           ylabel('Rate, spikes/s')
		   end
		   title(dtit)
	
		   axlim = axis;
		   xleg = axlim(1)+0.5*(axlim(2)-axlim(1));
	       yleg = axlim(3) + 0.95*(axlim(4)-axlim(3));
		   ydlt = 0.05*(axlim(4)-axlim(3));
		   for j=1:nfilsin
			  text(xleg, yleg, legend(j));
			  yleg = yleg - ydlt;
		   end
	   
% For TDT filerot stimuli, more detailed displays. Determine setimulus set
% (Q, R, . . U or V, W, . . . Z) and sampling rate, get sigma value to display.
% tmult is the multiplier by which times (e.g. sigma and T) change due to
% changes in the sampling rate from its nominal value (80 or 160 kHz).
		   if ftype==2			% IF #3
		      [srate, stimset] = getsrate(pb.ltext);
			  setndex = stimset - 'P';
			  nomsamp = 80.;
			  if setndex>5
			     setndex=setndex-5;
				 nomsamp = 160.;
			  end
			  tmult = nomsamp/srate;
			  
			  showndex = 7;
			  while showndex~=0 & showndex~=6
				 if showndex<1 | showndex>5
			        fprintf(1,'\nFor stimulus set %c:\n',stimset)
			        fprintf(1,'  sigma = %g(1), %g(2), %g(3), %g(4), %g(5).\n', ...
			           tmult*sigmav(setndex,:))
			        while showndex<0 | showndex>6
			           showndex = input('Show PSTs for which sigma (1-5, 0 whole set, 6 quit)? ');
				    end
				 end
			
		         if showndex>0 & showndex<6		% IF #4
	                ret = plotpart(ftype, setndex, showndex, nbins, nbnstim, nfilsin, ...
				        legend, tmult, sigmav(setndex,showndex), absc, binw);
					showndex = ret;
					if ret=='i' | ret=='I' | ret=='l' | ret=='L'
						showndex = 7;
					elseif ret>='0' & ret<='6'
						showndex = ret - '0';
				    else
					    showndex = 6;
					end
				 end	 % ENDIF #4
			  end		 % END WHILE showndex~=0 & ~=6
		   else
			  showndex = 6;
		   end	   		 % ENDIF #3
	   end				 % END WHILE showndex~=6
	end				 % ENDIF #1 
return				 % from function





function [srate, stimset] = getsrate(str);

% srate = GETSRATE(str) returns sampling rate for TDT filrot stimuli from
% lab parameter specification line. e.g.
% gtsrate('S.XR4W(152.), . . .')       returns srate=152, stimset='W'

	n1 = 1;
	n2 = size(str, 2);
	A = 0;
	while char(A)~='('
		stimset = A;
		[A,nwd,errm,nch]=sscanf(str(n1:n2),'%c',1);
		n1 = n1 + nch - 1;
	end
	srate = sscanf(str(n1:n2), '%g', 1);
return







function [ftype, binw, nbins, nbnstim, npsts] = ppblpst(pb, binset, binwset);

% Parse parameter block in pb for PST histogram computation. Returns:
% ftype = file type = 0 - normal PST, constant stimulus over all lines
%                   = 1 - sweep stimulus, with large number of different
%						  stimuli (i.e. freq. sweep, rate-level, etc.)
%					= 2 - TDT file rotation stimulus, from STR with smaller
%						  number of different stimuli.
% binw - recommended PST binwidth.
% nbins - recommended number of bins.
%     if binset = 1, binw = binwset (at least 0.1 ms) and nbins is
%                    min(2000, 2*dur/binw)
%               = 0, binw = dur/250 and nbins = 500.
% nbnstim - number of bins during the stimulus (dur).
% npsts - number of psts in the file. There is one PST for each different
%    stimulus. NOTE that if ftype=1, then 100/npsts successive lines should be
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
% Also limited to stimulus period - 10 ms to avoid an error in computing the
% PST (see lineleng<pstleng in cmppst).
	if binset~=0
		binw = max([binwset, 0.01]);
	else
		binw = round(10*dur/250)/10;
	end
	nbins = round(min([4000, 2*dur/binw, (perd-10)/binw]));
	nbnstim = min([nbins, round(dur/binw)]);

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







function comppst(resol, ftype, binw, nbins, npsts);

% Compute pst histograms for file open on FILEPNTRXYZ. Histogram has binwidth
% binw ms, data have time resolution resol microsec. Each histogram has nbins
% bins with latencies 0 through nbins*binw. npsts histograms are made.
% If ftype=1,
%     lines are combined in chunks of 100/npsts (1-5, 6-10, etc.)
% if ftype=2,
%     every npsts'th line is combined (1,6,11, ...; 2, 7, 12, ...; etc.)
% Each PST histogram is a column vector of PSTDATXYZ.
% NLINFXZY(j) is row vector of the number of lines of data included in the
% jth pst, i.e. in PSTDATXYZ(:,j). Usually, nlinf = 100/npsts. NOTE that a
% line is not included in the PST unles its length exceeds the abscissa length
% of the PST.


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
	   filtw = max([1, 1+2*fix(filtw/2)]);
   end
return




function nch = lens(str);

% Returns the length of the string, up to the last printable non-blank character

	for j=size(str,2):-1:1
		nch = j;
		if abs(str(j))>abs(' ')
			break
		end
	end
return





function ret = plotpart(ftype, setndex, showndex, nbins, nbnstim, nfilsin, legend, ...
    tmult, thsigma, absc, binw);

% Displays a portion of the PSTs in PSTDATXYZ, selecting those with one sigma
% value, as set by showndex.
% ret returns the character entered in response to the what next prompt.

global PSTDATXYZ

% DATA STRUCTURES containing layout of click trains for Hafter expt.
% XXcols are row vectors that select the PST columns in PSTDATXYZ that
% contain data for each sigma value. XXnums are the meaningful dimensions
% of the row vectors in XXcols, XXtvals, and XXfcli. XXtvals are the T values
% for each set of PSTs. The XXfcli are the latencies of the center of the
% first click in each stimulus.
    QVcols = [  1,  3,  6, 10, 15, 20, 25,
	            2,  5,  9, 14, 19, 24, 25,
		        4,  8, 13, 18, 23, 25, 25,
		  	    7, 12, 17, 22, 25, 25, 25
	           11, 16, 21, 25, 25, 25, 25];
	QVnums = [7 6 5 4 3];
	QVtvals = [0.5,    1,    2,  4.5,   10, 21.5, 46.5,
	             1,    2,  4.5,   10, 21.5, 46.5,    0,
	             2,  4.5,   10, 21.5, 46.5,    0,    0,
	           4.5,   10, 21.5, 46.5,    0,    0,    0,
                10, 21.5, 46.5,    0,    0,    0,    0];
	QVfcli = [10.5,  11.,  12., 14.5,  20., 20.75, 33.25,
			   11.,  12., 14.5,  20., 20.75,33.25,    0,
			   12., 14.5,  20., 20.75,33.25,    0,    0,
			  14.5,  20., 20.75,33.25,    0,    0,    0,
			   20., 20.75,33.25,    0,    0,    0,    0];

	RWcols = [  1,  3,  6, 10, 15, 20, 25,
                2,  4,  7, 11, 16, 21, 25,
                5,  9, 14, 19, 24, 25, 25,
                8, 13, 18, 23, 25, 25, 25,
	           12, 17, 22, 25, 25, 25, 25];
	RWnums = [7 6 5 4 3];
	RWtvals = [ 0.5,    1, 2.25, 4.75,   10, 21.5, 46.5,
	            0.5,    1, 2.25, 4.75,   10, 21.5,    0,
	           2.25, 4.75,   10, 21.5, 46.5,    0,    0,
	           4.75,   10, 21.5, 46.5,    0,    0,    0,
	             10, 21.5, 46.5,    0,    0,    0,    0];
	RWfcli = [ 10.5,  11., 12.25,14.75, 20., 20.75,33.25,
	           10.5,  11., 12.25,14.75, 20., 20.75,   0.,
			  12.25,14.75,   20.,20.75,33.25,   0.,   0.,
			  14.75,  20., 20.75,33.25,   0.,   0.,   0.,
			    20.,20.75, 33.25,   0.,   0.,   0.,   0.];

	SXcols = [  2,  5,  8, 12, 17, 22,
                1,  4,  7, 11, 16, 21,
                3,  6, 10, 15, 20, 25,
                9, 14, 19, 24, 25, 25,
	           13, 18, 23, 25, 25, 25];
	SXnums = [6 6 6 4 3];
	SXtvals = [0.5,     1, 2.125, 4.625,   10,   21.5,
	           0.5,     1, 2.125, 4.625,   10,   21.5,
	             1, 2.125, 4.625,    10, 21.5, 46.375,
	         4.625,    10,  21.5,46.375,    0,      0,
	            10,  21.5,46.375,     0,    0,      0];
	SXfcli = [ 10.5,  11., 12.125,14.625,  20., 20.75,
	           10.5,  11., 12.125,14.625,  20., 20.75,
			    11.,12.125,14.625,   20.,20.75,33.1875,
			 14.625,  20.,  20.75,33.1875,   0.,   0.,
			    20.,20.75, 33.1875,   0.,    0.,   0.];
				
	TYcols = [  2,  5,  8, 12, 17, 22,
                1,  4,  7, 11, 16, 21,
                3,  6, 10, 15, 20, 25,
                9, 14, 19, 24, 25, 25,
	           13, 18, 23, 25, 25, 25];
	TYnums = [6 6 6 4 3];
	TYtvals = [0.4375,     1, 2.125, 4.625,   10,   21.5,
	           0.4375,     1, 2.125, 4.625,   10,   21.5,
	                1, 2.125, 4.625,    10, 21.5, 46.375,
	            4.625,    10,  21.5,46.375,    0,      0,
	               10,  21.5,46.375,     0,    0,      0];
	TYfcli = [10.4375,   11.,12.125,14.625,  20., 20.75,
	          10.4375,   11.,12.125,14.625,  20., 20.75,
			      11.,12.125,14.625,   20.,20.75,33.1875,
			  14.625,    20., 20.75,33.1875,  0.,   0.,
			     20.,  20.75,33.1875,   0.,   0.,   0.];

	UZcols = [  2,  5,  8, 12, 17, 22,
                1,  4,  7, 11, 16, 21,
                3,  6, 10, 15, 20, 25,
                9, 14, 19, 24, 25, 25,
	           13, 18, 23, 25, 25, 25];
	UZnums = [6 6 6 4 3];
	UZtvals = [0.4375,     1, 2.125, 4.625,   10, 21.5,
	           0.4375,     1, 2.125, 4.625,   10, 21.5,
	                1, 2.125, 4.625,    10, 21.5, 46.4,
	            4.625,    10,  21.5,  46.4,    0,    0,
	               10,  21.5,46.375,     0,    0,    0];
	UZfcli = [10.4375,   11.,12.125,14.625,  20., 20.75,
	          10.4375,   11.,12.125,14.625,  20., 20.75,
			      11.,12.125,14.625,   20.,20.75,33.1875,
			  14.625,    20., 20.75,33.1875,  0.,   0.,
			     20.,  20.75,33.1875,   0.,   0.,   0.];


% numpsts is the number of PSTs for this sigma set. savdat contains only
% those PSTs. tvals contains the inter-click intervals (T). fclat contains
% the first-click latencies
     switch setndex
         case 1
			numpsts = QVnums(showndex);
	        savdat = PSTDATXYZ(:,QVcols(showndex,1:numpsts));
			tvals = tmult*QVtvals(showndex, 1:numpsts);
			fclat = QVfcli(showndex, 1:numpsts);
		 case 2
			numpsts = RWnums(showndex);
	        savdat = PSTDATXYZ(:,RWcols(showndex,1:numpsts));
			tvals = tmult*RWtvals(showndex, 1:numpsts);
			fclat = RWfcli(showndex, 1:numpsts);
		 case 3
			numpsts = SXnums(showndex);
	        savdat = PSTDATXYZ(:,SXcols(showndex,1:numpsts));
			tvals = tmult*SXtvals(showndex, 1:numpsts);
			fclat = SXfcli(showndex, 1:numpsts);
		 case 4
			numpsts = TYnums(showndex);
	        savdat = PSTDATXYZ(:,TYcols(showndex,1:numpsts));
			tvals = tmult*TYtvals(showndex, 1:numpsts);
			fclat = TYfcli(showndex, 1:numpsts);
		 case 5
			numpsts = UZnums(showndex);
	        savdat = PSTDATXYZ(:,UZcols(showndex,1:numpsts));
			tvals = tmult*UZtvals(showndex, 1:numpsts);
			fclat = UZfcli(showndex, 1:numpsts);
	 end

% Offset chosen PSTs for plotting. Set one display value in top histogram
% to the max value, to control scaling (not necessary any more, see axis
% call after plot).
	 mxrate = max([0.5, max(max(savdat))]);
	 addto = mxrate*ones(nbins,1)*[0:numpsts-1];
	 displdat = savdat + addto;

% Make a matrix of abscissae with the first click time lined up on time 0.
% zroindx points to bin with 0 time, used later.
	 [absc1, zroindx] = absczero(fclat, binw, nbins, numpsts);

% Display title
	  dtit = sprintf('For sigma = %g (nominally %g) ms, binw = %g ms.', ...
	      tmult*thsigma, thsigma, binw);

% Display the PSTs
	 nbindisp = fix(1.5*nbnstim);
	 ret = 'x';
	 while ret=='x' | ret=='X'
	    dispPST(numpsts, absc1, displdat, nbindisp, nbnstim, binw, mxrate, ...
           ftype, tmult, thsigma, tvals, dtit)

% Make integrated plot, if desired
	   ret = input('\nInteg plot, Log-log integ plot, Xpand, or next sigma(0-6): ','s');

	   if ret=='i' | ret=='I' | ret=='l' | ret=='L'
	      dispinteg(ret, savdat, numpsts, absc1, nbins, nbnstim, binw, ftype, ...
             dtit, tvals, zroindx);
	   elseif ret=='x' | ret=='X'
	      newmax = input('\nNew max abscissa (ms) = ');
		  nbindisp = fix(newmax/binw);
	   end
	end
return





function dispPST(numpsts, absc1, displdat, nbindisp, nbnstim, binw, mxrate, ...
    ftype, tmult, thsigma, tvals, dtit)
	 
% Display PSTs, labelled with T values. nbindisp bins are displayed.
	  for j=1:numpsts
         plot(absc1(1:nbindisp,j), displdat(1:nbindisp,j),'b-')
		 if j==1; hold on; end
	  end
	  hold off
	  axis([-10. min([nbindisp*binw 1.5*nbnstim*binw]) 0 numpsts*mxrate]);
      xlabel('Milliseconds re first-click center')
	  if ftype==2
	     ylabel('Number of spikes/stimulus')
	  else
	     ylabel('Rate, spikes/s')
	  end
      title(dtit)

      axlim = axis;
	  xleg = axlim(1) + 0.80*(axlim(2)-axlim(1));
	  ydlt = mxrate;
	  yleg = axlim(3) + 0.8*ydlt;
	  for j=1:numpsts
		  time = sprintf('T=%g ms.', tvals(j));
		  text(xleg, yleg, time)
		  yleg = min([yleg + ydlt, axlim(3)+0.95*(axlim(4)-axlim(3))]);
	  end
	  disptick(axlim, mxrate, numpsts, tvals, nbnstim, binw);
return




function dispinteg(ret, savdat, numpsts, absc1, nbins, nbnstim, binw, ftype, ...
     dtit, tvals, zroindx)

% Display integrated plot (ret='i' or 'I') or log-log integ plot
% (ret='l' or 'L')

% Plot colors
	plcol = {'b-', 'g-', 'r-', 'c-', 'm-', 'k-', 'y-', 'b:', 'g:'};
	lncol = {'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g'};

% Compute integrated spike count function.
	cumdat = cumsum(savdat);

	if ret=='i' | ret=='I'
       for j=1:numpsts
          plot(absc1(:,j), cumdat(:,j), char(plcol(j)))
	      if j==1; hold on; end
       end
       hold off
	   axlim = axis;
	   axis([-10. 1.5*nbnstim*binw 0 axlim(4)]);
       xlabel('Milliseconds re first click center')
       if ftype==2
          ylabel('Number of spikes/stimulus')
       else
          ylabel('Rate, spikes/s')
	   end
	   title(dtit)
	   axlim = axis;
	   xleg = axlim(1) + 0.80*(axlim(2)-axlim(1));
	   xdlt = xleg - 0.05*(axlim(2) - axlim(1));
       ydlt = 0.05*(axlim(4) - axlim(3));
       yleg = axlim(3) + ydlt;
       for j=1:numpsts
	      time = sprintf('T=%g ms.', tvals(j));
		  line([xdlt, xleg], [yleg, yleg], 'Color', char(lncol(j)));
	      text(xleg, yleg, time)
	      yleg = yleg + ydlt;
       end
    else
% Loglog plot plots integrated number of spikes vs log(1+(t-t1st_click)/T)
% which is like #spikes vs. click number. Note that absc1 is already t-t1st_click.
       absc2 = ones(nbins,numpsts) + absc1./(ones(nbins,1)*tvals);
	   for j=1:numpsts
		   loglog(absc2(zroindx(j):nbins,j), ...
		      cumdat(zroindx(j):nbins,j), char(plcol(j)))
	       if j==1; hold on; end
       end
       hold off
       xlabel('log[1 + (ms re first click)/T]')
       if ftype==2
          ylabel('log(Number of spikes/stim.)')
       else
          ylabel('log(Rate, spikes/s)')
	   end
	   title(dtit)
	   axlim = axis;
	   xleg = axlim(1)*exp(0.80*log(axlim(2)/axlim(1)));
	   xdlt = axlim(1)*exp(0.75*log(axlim(2)/axlim(1)));
       ydlt = exp(0.05*log(axlim(4)/axlim(3)));
       yleg = axlim(3)*ydlt;
       for j=1:numpsts
	      time = sprintf('T=%g ms.', tvals(j));
		  line([xdlt, xleg], [yleg, yleg], 'Color', char(lncol(j)));
	      text(xleg, yleg, time)
	      yleg = yleg*ydlt;
	  end
   end
return




function [absc1, zroindx] = absczero(fclat, binw, nbins, numpst);

% fclat is row vector of first click latencies for a set of PSTs with same sigma.
% binw is the histogram binwidth.
% The PSTs have nbins bins and there are numpst of them (i.e. the
% PST data matrix is nbins x numpst). It is assumed that the true abscissa takes
% the form
%      {bw/2, 3*bw/2, . . . (2n+1)*bw/2, . . .
% This routine returns a matrix of new abscissae, absc1, for each PST such that
% time zero on each abscissa is the center of the first click. Also returns a
% row vector zroindx of indices of the bin containing these 0 ms.

	absc1 = zeros(nbins, numpst);
	zroindx = fix(fclat/binw) + ones(1,size(fclat,2));
	absc1(1,:) = -binw*fix(fclat/binw);
	for j=2:nbins
	   absc1(j,:) = absc1(j-1,:) + binw;
    end
return
		
		
		
		
		
function disptick(axlim, mxrate, numpsts, tvals, nbnstim, binw)

% Display small ticks below the abscissae of PST plot at the latencies of
% the centers of the clicks.

	dur = nbnstim*binw;
	for j=1:numpsts
		ncl = fix(dur/tvals(j));
		xmat = [1 1]'*tvals(j)*[0:ncl-1];
		ymat = [axlim(3)+(j-1)*mxrate, ...
		   axlim(3)+(j-1)*mxrate-0.02*(axlim(4)-axlim(3))]'*ones(1,ncl);
		line(xmat, ymat, 'Color', 'm')
	end
return

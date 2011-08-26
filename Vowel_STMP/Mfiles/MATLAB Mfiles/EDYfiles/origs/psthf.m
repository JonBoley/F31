	function [absc] = psthf(str);
% PSTHF('a4-6,9') - computes and displays PST histograms of data in the files
% a4.DAT, a5.DAT, a6.DAT, a9.DAT which should be in the lab standard format.
% PST('a4b0.5f7') computes a PST histogram with binwidth 0.5 ms and filters
% it with a triangular window of width 7 bins.
% PST('a4h') assumes data are Hafter TDT click data.
% Filename must be first, but binwidth and filterwidth can be in either order.

% NOTE: in order to use the TDT click part of the program, it must be called
% as follows:  global INDEXXYZ;psthf(...).
% The easiest way is to call the script psth.

% Modified 2/24/98 to get click train parameters from a script genhafter
% which is in the same folder as the data.

% Modified 2/11/99 to assume is not TDT click data, unless "h" is in the call.

% NOTE THIS PROGRAM ASSUMES THAT TDT CLICK TRAIN DATA ARE TAKEN WITH 4 REPS
% OF 25 STIMULI. IF NOT TAKEN THIS WAY, THEN ANALYSIS BY SIGMA GROUP FAILS.

global PSTDATXYZ NLINFXYZ INDEXXYZ

	legend = {'', '', '', '', '', '', '', '', '', '', ...
	          '', '', '', '', '', '', '', '', '', ''};
			  
% Parse the input command string:
	[fname, nfils, binset, binwset, autofilt, filtw, dohaf] = parscmpst(str);
	
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
% for (Hafter) TDT click trains
	   if ftype~=2 | dohaf==0
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
		   figure(1);
		   plot(absc, displdat,'b-')
		   xlabel('Milliseconds')
		   if ftype==2 & dohaf==1
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

		   if ftype==2 & dohaf==1		% IF #3
			   
% First, get the sigmav and hafterfudge values for these data.
		      ishafter = 1;
			  genhafter;
			  if ishafter~=0 | strcmp(hafterdate, pb.date)~=1
				  fprintf(1,'***WARNING: click train params are probably wrong.\n');
				  fprintf(1,'            most likely, the genhafter file is not in the data folder.\n');
			  end
		  
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
			           tmult*hafterfudge*sigmav(setndex,:))
			        while showndex<0 | showndex>6
			           showndex = input('Show PSTs for which sigma (1-5, 0 whole set, 6 quit)? ');
				    end
				 end
			
		         if showndex>0 & showndex<6		% IF #4
	                ret = plotpart(ftype, setndex, showndex, nbins, nbnstim, nfilsin, ...
				        legend, tmult, sigmav(setndex,showndex), absc, binw, ...
						hafterfudge, dtit);
					showndex = ret;
					if ret=='i' | ret=='I' | ret=='l' | ret=='L' | ret=='n' | ret=='N'
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
%     if binset = 1, binw = binwset (at least 0.01 ms) and nbins is
%                    min(20000, 2*dur/binw)
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
% number of bins is sufficient to make PST cover 2*dur, limited to 20000 bins
% Also limited to stimulus period - 10 ms to avoid an error in computing the
% PST (see lineleng<pstleng in cmppst).
	if binset~=0
		binw = max([binwset, 0.01]);
	else
		binw = round(10*dur/250)/10;
	end
	nbins = round(min([20000, 2*dur/binw, (perd-10)/binw]));
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
% line is not included in the PST unless its length exceeds the abscissa length
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



function [fname, nfils, binset, binwset, autofilt, filtw, dohaf] = parscmpst(str);

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
% dohaf = 0 if not doing TDT clicks; 1 if doing TDT clicks

% fnames holds the filenames to be processed (up to 20)
	fname={'', '', '', '', '', '', '', '', '', '', ...
         	'', '', '', '', '', '', '', '', '', ''};
			

% binset = 1 if binwidth set in call string, otherwise determined from file
% autofilt = 1 if filtering specified in call string.
	binset = 0;
	autofilt = 0;
	binwset = 0.1;
	filtw = 1;
	dohaf = 0;
	
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
			   elseif char(A(1))=='h' | char(A(1))=='H'
				   dohaf = 1;
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
    tmult, thsigma, absc, binw, hafterfudge, dtite);

% Displays a portion of the PSTs in PSTDATXYZ, selecting those with one sigma
% value, as set by showndex.
% ret returns the character entered in response to the what-next prompt.
% ftype = file type (2 for TDT stimuli)
% setndex = {1:5} for stimulus sets {Q:U} or {V:Z}
% showndex = indexes sigma value for which data are to be shown
% nbins = number of bins in PST
% nbnstim = duration of stimulus in PST bins
% nfilsin = number of files averaged into the PST
% legend = legend text for displays
% tmult = multiplies click times to correct for sampling rate changes
% thsigma = the sigma value, ms (uncorrected for sampling rate issues)
% absc = abscissa of the PST
% binw = PST binwidth, ms
% hafterfudge = multiplies click times to correct for hw/sw error
% dtite is the display title: '4/10.97, 1.04, P.14, bin width = 0.2 ms'


global PSTDATXYZ INDEXXYZ

% Get the data structures containing layout of click trains for Hafter expt.
% XXcols are row vectors that select the PST columns in PSTDATXYZ that
% contain data for each sigma value. XXnums are the meaningful dimensions
% of the row vectors in XXcols, XXtvals, and XXfcli. XXtvals are the T values
% for each set of PSTs. The XXfcli are the latencies of the center of the
% first click in each stimulus.
	 ishafter = 2;
	 genhafter;
	 if ishafter~=0
		 fprintf(1,'***WARNING: click train params are probably wrong.\n');
		 fprintf(1,'            most likely, the genhafter file is not in the data folder.\n');
	 end



% numpsts is the number of PSTs for this sigma set. savdat contains only
% those PSTs. tvals contains the inter-click intervals (T). fclat contains
% the first-click latencies
     switch setndex
         case 1
			numpsts = QVnums(showndex);
	        savdat = PSTDATXYZ(:,QVcols(showndex,1:numpsts));
			tvals = hafterfudge*tmult*QVtvals(showndex, 1:numpsts);
			fclat = hafterfudge*tmult*QVfcli(showndex, 1:numpsts);
		 case 2
			numpsts = RWnums(showndex);
	        savdat = PSTDATXYZ(:,RWcols(showndex,1:numpsts));
			tvals = hafterfudge*tmult*RWtvals(showndex, 1:numpsts);
			fclat = hafterfudge*tmult*RWfcli(showndex, 1:numpsts);
		 case 3
			numpsts = SXnums(showndex);
	        savdat = PSTDATXYZ(:,SXcols(showndex,1:numpsts));
			tvals = hafterfudge*tmult*SXtvals(showndex, 1:numpsts);
			fclat = hafterfudge*tmult*SXfcli(showndex, 1:numpsts);
		 case 4
			numpsts = TYnums(showndex);
	        savdat = PSTDATXYZ(:,TYcols(showndex,1:numpsts));
			tvals = hafterfudge*tmult*TYtvals(showndex, 1:numpsts);
			fclat = hafterfudge*tmult*TYfcli(showndex, 1:numpsts);
		 case 5
			numpsts = UZnums(showndex);
	        savdat = PSTDATXYZ(:,UZcols(showndex,1:numpsts));
			tvals = hafterfudge*tmult*UZtvals(showndex, 1:numpsts);
			fclat = hafterfudge*tmult*UZfcli(showndex, 1:numpsts);
	 end

% Offset chosen PSTs for plotting. Set one display value in top histogram
% to the max value, to control scaling (not necessary any more, see axis
% call after plot).
	 mxrate = max([0.5, max(max(savdat))]);
	 displdat = savdat + mxrate*ones(nbins,1)*[0:numpsts-1];

% Make a matrix of abscissae with the first click time lined up on time 0.
% zroindx points to bin with 0 time, used later.
% minsigndx points to bin that is two sigmas back from first click, used later.
% minonendx points to bin that is 1 ms back from first click, used later.
	 [absc1, zroindx, minsigndx, minonendx] = absczero(fclat, ...
	       tmult*hafterfudge*thsigma, binw, nbins, numpsts);
		   
% Display title
	  dtit = strcat(dtite, sprintf('For sigma = %g (nom. %g) ms.', ...
	      tmult*hafterfudge*thsigma, thsigma));

% Display the PSTs
	 nbindisp = fix(1.5*nbnstim);
	 ret = 'x';
	 while ret=='x' | ret=='X'
		figure(1);
	    dispPST(numpsts, absc1, displdat, nbindisp, nbnstim, zroindx, binw, ...
		   mxrate, ftype, thsigma, tvals, dtit)

% Make integrated plot, if desired
	   ret = input('\nLog-log integ plot, log(N) plot, Xpand, or next sigma(0-6): ','s');

	   if ret=='l' | ret=='L'
		  figure(1);		  
		  dispinteg(savdat, numpsts, nbins, binw, ftype, dtit, tvals, fclat, mxrate)
	   elseif ret=='x' | ret=='X'
	      newmax = input('\nNew max abscissa (ms) = ');
		  nbindisp = fix(newmax/binw);
	   elseif ret=='n' | ret=='N'
		  ret = logncalc(savdat, numpsts, fclat, binw, tvals, nbnstim, dtit, ...
		        ftype, mxrate);
	   end
	end
return





function dispPST(numpsts, absc1, displdat, nbindisp, nbnstim, zroindx, binw, ...
     mxrate, ftype, thsigma, tvals, dtit)
	 
% Display PSTs, labelled with T values. nbindisp bins are displayed.
	  for j=1:numpsts
         plot(absc1(1:zroindx(j)+nbindisp,j), ...
		            displdat(1:zroindx(j)+nbindisp,j),'b-')
		 if j==1; hold on; end
	  end
	  hold off
	  axis([-10. min([nbindisp*binw 1.5*nbnstim*binw]) ...
	       -0.1*mxrate numpsts*mxrate]);
      xlabel('Milliseconds re first-click center')
	  if ftype==2
	     ylabel('Number of spikes/stimulus')
	  else
	     ylabel('Rate, spikes/s')
	  end
      title(dtit)

      axlim = axis;
	  xleg = axlim(1) + 0.8*(axlim(2)-axlim(1));
	  ydlt = mxrate;
	  yleg = axlim(3) + 0.8*ydlt;
	  for j=1:numpsts
		  time = sprintf('T=%g ms.', tvals(j));
		  text(xleg, yleg, time)
		  yleg = min([yleg + ydlt, axlim(3)+0.95*(axlim(4)-axlim(3))]);
	  end
	  disptick(axlim, mxrate, numpsts, tvals, nbnstim, binw);
return




function dispinteg(savdat, numpsts, nbins, binw, ftype, dtit, tvals, fclat, ...
     mxrate)

% Display integrated rate plot on log-log coordinates. Latency compensation is
% done by user interaction.

% Plot colors
	plcol = {'b-', 'g-', 'r-', 'c-', 'm-', 'k-', 'y-', 'b:', 'g:'};
	lncol = {'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g'};

% Use screen cursor to select latency of the click response. latindex is a vector
% of indices of bins containing the cursor. The actual response latencies should be
% the right-hand edges of these bins.
	fprintf(1,'\n>>>>> Set cursor at start latency and click. <<<<<')
	[sttime, yval] = cursor(gcf);
%	axlim = axis;
%	line([sttime sttime], [axlim(3) axlim(4)])
	axlim = axis;
	disptick2(axlim, mxrate, numpsts, sttime, tvals, binw, 0, 'g')
	drawnow

% Compute integrated spike count function and integrated count starting with 0 count
% at first-click latency, entered above (note this latency can be <0).
	latindex = 1+fix((fclat+sttime)/binw);
	cumdat = cumsum(savdat);
	cumdatzro = zeros(nbins,numpsts);
	for j=1:numpsts
		cumdatzro(:,j) = cumdat(:,j) - cumdat(latindex(j),j)*ones(nbins,1);
	end

% Loglog plot plots integrated number of spikes vs log(1+(t-t_latency)/T), computed
% in absc2. This is like cumulative #spikes vs. click number.
	absc2 = zeros(nbins, numpsts);
	absc2(1,:) = binw*(0.5*ones(1,numpsts)-latindex);
	for j=2:nbins
		absc2(j,:) = absc2(j-1,:) + binw;
	end
	absc2 = ones(nbins,numpsts) + absc2./(ones(nbins,1)*tvals);

	fprintf(1,'\n>>>>> Hit CR to continue. <<<<<\n')
	pause

	for j=1:numpsts
	   loglog(absc2(latindex(j):nbins,j), ...
	      cumdatzro(latindex(j):nbins,j), char(plcol(j)))
	   if j==1; hold on; end
	end
	hold off
	axlim = axis;
	axis([0.9 axlim(2) min([0.9 max([0.1 axlim(3)])]) axlim(4)])
	xlabel('log[1 + (ms re 1st resp. latency)/T]')
	if ftype==2
	  ylabel('log(Cumulative number of spikes/stim.)')
	else
	  ylabel('log(Cumulative rate, spikes/s)')
	end
	title(dtit)
	axlim = axis;
	xleg = axlim(1)*exp(0.80*log(axlim(2)/axlim(1)));
	xdlt = axlim(1)*exp(0.75*log(axlim(2)/axlim(1)));
	ydlt = exp(0.05*log(axlim(4)/axlim(3)));
	yleg = axlim(3)*ydlt;
	for j=1:numpsts
	  line([xdlt, xleg], [yleg, yleg], 'Color', char(lncol(j)));
	  text(xleg, yleg, sprintf('T=%g ms.', tvals(j)))
	  yleg = yleg*ydlt;
    end
    line([1 min([axlim(2) axlim(4)])], [1 min([axlim(2) axlim(4)])], ...
	    'Color', 'k', 'LineStyle', '--');
return




function [absc1, zroindx, minsigndx, minonendx] = absczero(fclat, thsigma, binw, ...
       nbins, numpst);

% fclat is row vector of first click latencies for a set of PSTs with same sigma.
% thsigma is that sigma value.
% binw is the histogram binwidth.
% The PSTs have nbins bins and there are numpst of them (i.e. the
% PST data matrix is nbins x numpst). It is assumed that the true abscissa takes
% the form
%      {bw/2, 3*bw/2, . . . (2n+1)*bw/2, . . .
% This routine returns a matrix of new abscissae, absc1, for each PST such that
% time zero on each abscissa is the center of the bin containing the first click.
% Also returns a row vector zroindx of indices of the bin containing the first
% clicks and a row vector minsigndx containing the point two sigmas back from
% the first click (i.e. fclat-thsigma) and a row vector minonendx which points
% to -1 ms re the first click.

	absc1 = zeros(nbins, numpst);
	zroindx = fix(fclat/binw) + ones(1,size(fclat,2));
	minsigndx = fix((fclat-2.*thsigma*ones(1,size(fclat,2)))/binw) + ...
	        ones(1,size(fclat,2));
	minonendx = fix((fclat-ones(1,size(fclat,2)))/binw) + ...
	        ones(1,size(fclat,2));
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
		ymat = [axlim(3)+(j-1)*mxrate+0.1*mxrate, ...
		   axlim(3)+(j-1)*mxrate+0.1*mxrate-0.02*(axlim(4)-axlim(3))]'*ones(1,ncl);
		line(xmat, ymat, 'Color', 'm')
	end
return



function ret = logncalc(savdat, numpsts, fclat, binw, tvals, nbnstim, dtit, ...
      ftype, mxrate)

% Set markers in on-screen PST display of responses at a sigma value. Compute
% number of spikes/stimulus versus click number between the markers.

global INDEXXYZ

% maxnpts is the maximum number of clicks to which the regression line is fit
	maxnpts = 100;
	
% Plot colors and legend space
	plcol = {'bo', 'gx', 'r*', 'cs', 'md', 'k^', 'yp', 'bh', 'g+'};
	plcol2 = {'b-', 'g-', 'r-', 'c-', 'm-', 'k-', 'y-', 'b:', 'g:'};
	thleg = {'','','','','','','','',''};

% Use screen cursor to select latencies for the click response.
	fprintf(1,'\n>>>>> Set cursor at start latency and click. <<<<<')
	[sttime, yval] = cursor(gcf);
%	axlim = axis;
%	line([sttime sttime], [axlim(3) axlim(4)])
	
	fprintf(1,'\n>>>>> Set cursor at end latency and click. <<<<<\n')
	[endtime, yval] = cursor(gcf);
%	line([endtime endtime], [axlim(3) axlim(4)])

% Compute the spike counts in response to each click in the train
	counts = zeros(100, numpsts);
	ncliks = zeros(1, numpsts);
	asttime = zeros(1, numpsts);
	aendtime = zeros(1, numpsts);
	axlim = axis;
	for j=1:numpsts
	   [counts(:,j), ncliks(j), asttime(j), aendtime(j)] = ncalc(savdat, j, sttime, ...
		     endtime, fclat(j), binw, tvals(j), nbnstim);	 
	   disptick2(axlim, mxrate, numpsts, asttime(j), tvals, binw, j, 'g')
	   disptick2(axlim, mxrate, numpsts, aendtime(j), tvals, binw, j, 'r')
	end
	
	fprintf(1,'\n>>>>> Hit CR to continue. <<<<<\n')
	pause

% Plot log(cum. number of spikes) versus log(number of clicks), along with line
% fit to first maxnpts data points.
	figure(1);
	ret = 'r';
	while ret=='r' | ret=='R'
		fitcount = zeros(100,numpsts);
		nclkfit = zeros(1, numpsts);
		slope = zeros(1, numpsts);
		offset = zeros(1, numpsts);
		regress = zeros(1, numpsts);
		for j=1:numpsts
		    cumcount = cumsum(counts(1:ncliks(j),j));
		    [fitcount(:,j), nclkfit(j), slope(j), offset(j), regress(j)] = ...
			                       logfit(cumcount, ncliks(j), maxnpts);
			loglog([1:ncliks(j)]', cumcount, char(plcol(j)))
			if j==1; hold on; end
			thleg(j) = strcat(thleg(j), sprintf('T=%g, lats=[%g,%g]ms, sl=%g', ...
			      tvals(j),asttime(j),aendtime(j),slope(j)));
		end
%		legend(char(thleg(1:numpsts)),2)
%		hold on
		for j=1:numpsts
		    loglog([1:nclkfit(j)]', fitcount(1:nclkfit(j),j), char(plcol2(j)))
		end  
		hold off
		doleg(thleg,numpsts);
		title(dtit)
		xlabel('Number of clicks')
		if ftype==2
	       ylabel('Cumulative number of spikes/click')
		else
	       ylabel('Cumulative rate, spikes/s')
	    end

	    fprintf(1,'\n           latency')
	    fprintf(1,'\n   T      start   end   slope  offset    r   fitting %g points.',maxnpts)
	    for j=1:numpsts
	      fprintf(1,'\n%7.3f%8.3f%c%6.3f%7.3f%8.3f%7.3f',tvals(j),asttime(j),'-', ...
		      aendtime(j),slope(j),offset(j),regress(j))
	    end
	    fprintf(1,'\n')
		
	    ret = input('Reset number of clicks fit, X-start again, or next sigma(0-6): ','s');
	    if ret=='r' | ret=='R'
		   maxnpts = input('Number of clicks fit by line: ');
		   thleg = cell(9);
	    end
	end
return
	
	
function [counts, ncliks, asttime, aendtime] = ncalc(savdat, nthPST, sttime, ...
     endtime, fclik, binw, thtval, nbnstim)

% Calculate cumulative spike counts in response to clicks for the PST savdat(:,nthPST).
% Spikes are cumulated over intervals
%    [sttime+j*thtval,  min([endtime sttime+thtval])+j*thtval   j=0, nmax
% (except that these intervals are adjusted to the left and right edges of bins)

% Returns:
%   counts - cumulative cound in savdat for each click (column)
%   ncliks - dimension of counts
%   [asttime aendtime] actual latency of first count interval. Successive ones
%        are spaced by ‰thtval, quantized by binw. NOTE that these are not the same
%	     as sttime and endtime.

% ncliks is the number of clicks in this PST
% aendtime is the endtime for counting, which will be reduced from endtime to prevent
%     overlap of counting windows.
% nbinst and nbinend are row vectors of pointers to PST bins over which counts are
%     accumulated
% asttime and aendtime are actual start and end times for first count, i.e. the
%     left and right edges of the bins

	ncliks = min([100, fix((nbnstim*binw-2*fclik)/thtval)]);
	aendtime = min([endtime sttime+thtval]);
	
	nbinst = 1 + fix(((fclik+sttime)*ones(1,ncliks) + thtval*[0:ncliks-1])/binw);
	asttime = (nbinst(1)-1)*binw - fclik;
	nbinend = 1 + fix(((fclik + aendtime)*ones(1,ncliks) + ...
	     thtval*[0:ncliks-1])/binw);
	if nbinend(1)>=nbinst(2)
		nbinend = nbinend - 1;
	end
	aendtime = nbinend(1)*binw - fclik;

	cumdat = cumsum(savdat(:,nthPST));
	counts = zeros(100,1);
	counts(1:ncliks) = cumdat(nbinend) - cumdat(nbinst-1);
return


function [fitcount, nclka, slope, offset, regress]=logfit(counts, ncliks, maxnpts)

% Fits a regression line to log(counts) versus log([1:ncliks]) for the first maxnpts
% points in counts. Returns:
%     fitcount - row vector of best fitting line evaluated at [1:ncliks] (fitcount
%           is linear as for counts, i.e. it is not log(counts))
%     nclkact - number of elements in fitcount. This is the number of data points
%           used in the fit.
%     slope, offset - line is
%         log(counts) = offset + slope*log(click number) + error
%	  regress - r value for the line in log-log space

	nclka = min([ncliks maxnpts]);
	xvals = log([1:nclka]);
	xmean = mean(xvals);
	sumx2 = (xvals-xmean)*(xvals-xmean)';
	ymean = mean(log(counts(1:nclka)));
	goob = (log(counts(1:nclka))-ymean);
	sumy2 = goob'*goob;
	slope = (xvals-xmean)*log(counts(1:nclka))/sumx2;
	offset = ymean - slope*xmean;
	regress = slope*sqrt(sumx2/sumy2);
	fitcount=zeros(100,1);
	fitcount(1:nclka) = exp(offset + slope*xvals)';
return


function doleg(thleg, numpsts)

% Plot the legend text on the screen. Assumes the screen contains a log-log plot.

	plcol = {'b', 'g', 'r', 'c', 'm', 'k', 'y', 'b', 'g'};

	axlim = axis;
	xleg1 = axlim(1)*(axlim(2)/axlim(1))^0.02;
	xleg2 = xleg1*(axlim(2)/axlim(1))^0.04;
	xleg3 = xleg2*(axlim(2)/axlim(1))^0.02;
	yleg = axlim(3)*(axlim(4)/axlim(3))^0.97;
	ydlt = (axlim(4)/axlim(3))^0.04;
	
	for j=1:numpsts
		line([xleg1, xleg2], [yleg, yleg], 'Color', char(plcol(j)));
		text(xleg3, yleg, char(thleg(j)))
		yleg = yleg/ydlt;
	end
return




function disptick2(axlim, mxrate, numpsts, sttime, tvals, binw, iwh, thcolor)

% Display ticks in the PST plot to mark the latencies [sttime, sttime+T, ...]
% For iwh=0, displays ticks for all PSTs
%        =1-numpsts, displays for only one PST

	dur = axlim(2) - sttime;
	for j=1:numpsts
		if iwh==0 | iwh==j
			ncl = 1 + fix(dur/tvals(j));
			xmat = [1 1]'*(sttime*ones(1,ncl)+tvals(j)*[0:ncl-1]);
			ymat = [axlim(3)+(j-1)*mxrate+0.6*mxrate, ...
			   axlim(3)+(j-1)*mxrate+0.6*mxrate-0.02*(axlim(4)-axlim(3))]'*ones(1,ncl);
			line(xmat, ymat, 'Color', thcolor)
		end
	end
return

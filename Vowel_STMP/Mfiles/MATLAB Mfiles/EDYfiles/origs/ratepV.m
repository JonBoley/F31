	function [drat,srat,abval,npts,nralvs]=ratepV(str)
% RATEPV('a4-7,15[s][f]') - computes and displays  rate versus level functions
% for driven and spontaneous rates from the files specified, which should
% be in standard lab format. Also reads PRE_GNUSYS files.
% Options:   s - superimposes files, otherwise are averaged together
%            f - automatically filters rate functions
% Uses matlab code (slow), but processes all file types correctly.
% Returns driven, spontaneous rates, abscissa values as 100xnralvs matrix,
% with number of points in jth column equal to npts(j). nralvs is the number
% of rate functions computed (columns in drat,srat,abval).
% NOTE: for convenience, call this as ratv.

global INDEXXYZ

% ispropen~=0 when log file is open for printing.
ispropen = 0;

plcolor=['k',  'r',  'g',  'm',  'b', 'c', 'k',  'r',  'g',  'm',  'b', 'c'];

% Allocate space for data

drda = zeros(100,20);
spda = zeros(100,20);
abda = zeros(100,20);
nlavg = zeros(20,1);
logs = zeros(20,1);
delts = zeros(20,1);
legend = {'', '', '', '', '', '', '', '', '', '', ...
     '', '', '', '', '', '', '', '', '', ''};

% thatsall is a loop-breakout flag
% nplots counts plots superimposed, controls plot appearance
% logabs~=0 if abscissa (of all plots) is logarithmic

thatsall = 0;
nplots = 0;
logabs = 0;
hold off

while thatsall==0           % #1

% Parse the input command string:

    [fname, nfils, autofilt, filtw, issuper] = parscmrat(str);

% Flag assumes non pre-GNUSYS data
	ispregnu = 0;
	
	if nfils<=0                               % #2
		thatsall = 1;
	else
		nplthis = nplots+1;      % Keep track of which files are in current group
		
% Open file and compute type of abscissa from the parameter block:

	for jf=1:nfils                            % #3
	   pb = labfileV(char(fname(jf)));        % <--- PREGNUSYS VERSION
	   if strcmp(pb.Error, 'No error')~=1     % #4
		  fprintf('\n*** Error = ''%s'' opening %s; check the directory.***\n', ...
		            pb.Error,char(fname(jf)))
		  fclose('all');
	   else
		
	      [absc, logabsth, delth, navrow, navskip, nplt, abname] = parsepb(pb);
		  if pb.prog_name=='TU' | pb.prog_name=='XC'
			  delth=pb.stim_par(17);          % <-- PREGNUSYS VERSION
		  end
	
	      if logabsth~=logabs & nplots>0    % Next plot is #nplots+1
		     fprintf(1,'*** WARNING: files have both linear and log abscissae.***')
	      else
		     logabs = logabsth;
	      end

% Compute rate windows using default values

		  dur = max([pb.stim_par(9), pb.stim_par(13)]);
		  timew = [10, dur+10, (dur + pb.mode_line(4))/2., pb.mode_line(4)];
   
% Compute rates (file closed here)
		  [rrate, rsrate, rlinf, lstlin] = comprate(pb.resol, timew);
	
% Average rates from same stimuli together

		  nplots = nplots + 1;		% COUNT THIS PLOT HERE
	
		  [rate, spont, linf, nlavg(nplots)] = ...
	           avgrate(rrate, rsrate, rlinf, lstlin, navrow, navskip, nplt);

% Construct title during first file processed
		  if nplots==1
			  titline=strcat(pb.date,', ', pb.unit(1:6));
	      end

%
% If pre-GNUSYS, separate into two units. Note that this involves incrementing
% nplots.
		  if pb.prog_name=='TU' | pb.prog_name=='XC'
			  
			  issuper = 1;        % <-- NOTE the use of this flag!
			  ispregnu = 1;       % Flag pre-GNUSYS data
			  
			  nlavg(nplots) = nlavg(nplots)-50;
			  drda(1:nlavg(nplots),nplots) = rate(1:nlavg(nplots));
			  spda(1:nlavg(nplots),nplots) = spont(1:nlavg(nplots));
			  logs(nplots) = logabsth; delts(nplots) = delth;
			  nplots = nplots +1;
			  nlavg(nplots) = nlavg(nplots-1);
			  drda(1:nlavg(nplots),nplots) = rate(51:50+nlavg(nplots));
			  spda(1:nlavg(nplots),nplots) = spont(51:50+nlavg(nplots));
			  logs(nplots) = logabsth; delts(nplots) = delth;
			  if pb.prog_name=='TU'
				  abda(1:nlavg(nplots-1),nplots-1) = ...
				     [absc(1):delth:absc(1)+(nlavg(nplots-1)-1)*delth]';
			  else
				  abda(1:nlavg(nplots-1),nplots-1) = [1:1:nlavg(nplots-1)]'; 
			  end
			  abda(1:nlavg(nplots),nplots) = abda(1:nlavg(nplots-1),nplots-1);
			  if autofilt==1             % <--- have to filter PREGNUSYS separately
				 drda(1:nlavg(nplots-1),nplots-1) = ...
				    trifilt(drda(1:nlavg(nplots-1),nplots-1)', filtw)';
				 spda(1:nlavg(nplots-1),nplots-1) = ...
				    trifilt(spda(1:nlavg(nplots-1),nplots-1)', filtw)';
				 drda(1:nlavg(nplots),nplots) = ...
				    trifilt(drda(1:nlavg(nplots),nplots)', filtw)';
				 spda(1:nlavg(nplots),nplots) = ...
				    trifilt(spda(1:nlavg(nplots),nplots)', filtw)';
			  end
			  legend(nplots-1) = strcat(legend(nplots-1), ...
	             sprintf('File %s, LW, %s,   R[%g,%g]&[%g,%g] ms', char(fname(jf)), ...
				 pb.ltext(1:lens(pb.ltext)), timew));
			  legend(nplots) =  strcat(legend(nplots), ...
	             sprintf('File %s, UP, %s,   R[%g,%g]&[%g,%g] ms', char(fname(jf)), ...
				 pb.ltext(1:lens(pb.ltext)), timew));
		      if autofilt==1
				 legend(nplots-1) = strcat(legend(nplots-1),', F3');
				 legend(nplots) = strcat(legend(nplots),', F3');				 
			  end
		  else			  
	%
	% If autofiltering,
			  if autofilt==1
				 rate = trifilt(rate', filtw)';
				 spont = trifilt(spont', filtw)';
			  end
	%
	% Collect results into data-pile
			  drda(1:nlavg(nplots),nplots) = rate(1:nlavg(nplots));
			  spda(1:nlavg(nplots),nplots) = spont(1:nlavg(nplots));
			  abda(1:nlavg(nplots),nplots) = absc(1:nlavg(nplots));
			  logs(nplots) = logabsth;
			  delts(nplots) = delth;
	% Including minimal identification of figure
			  legend(nplots) = strcat(legend(nplots), ...
	             sprintf('File %s, %s,   R[%g,%g]&[%g,%g] ms', char(fname(jf)), ...
				 pb.ltext(1:lens(pb.ltext)), timew));
		      if autofilt==1
				 legend(nplots) = strcat(legend(nplots),', F3');
			 end
		  end
	   end          % #4 end of processing for this file
	end             % #3 end of loop processing the files requested

% Start of display portion. First check if there are no files. Then check
% if files in the current group [nplthis,nplots] are to be averaged.

	if nplots>0     % #5 (no data, e.g. errors in opening all files)
		if nplots-nplthis>=1 & issuper==0
			[arate, asrate, aabsc, nlevs, aleg] = avgplots(drda, spda, ...
				       abda, nlavg, delts, logs, legend, nplthis, nplots);
			drda(1:nlevs, nplthis) = arate(1:nlevs);
			spda(1:nlevs, nplthis) = asrate(1:nlevs);
			abda(1:nlevs, nplthis) = aabsc(1:nlevs);
			nlavg(nplthis) = nlevs;
			legend(nplthis) = aleg;
			for jl=nplthis+1:nplots; legend(jl)={''}; end
			nplots = nplthis;
		end
	
% Start of display loop. Display until displayloop~=0. iscursrate~=0
% when there is a cursor display line.

		displayloop = 0;
		iscursrate = 0;
		
		while displayloop==0               % #6 displayloop
		   for j=1:nplots
		      xax = [100000 -100000];
	          if logabs == 0
                 plot(abda(1:nlavg(j),j), drda(1:nlavg(j),j), strcat(plcolor(j),'-'), ...
			       abda(1:nlavg(j),j), spda(1:nlavg(j),j), strcat(plcolor(j),'--'))
	          else
		         semilogx(abda(1:nlavg(j),j), drda(1:nlavg(j),j), strcat(plcolor(j),'-'), ...
			       abda(1:nlavg(j),j), spda(1:nlavg(j),j), strcat(plcolor(j),'--'))
			     xax(1) = min([xax(1), abda(1,j)]);
			     xax(2) = max([xax(2), abda(nlavg(j),j)]);
	          end
		      if j==1; hold on; end
		   end
		   hold off
		   xlabel(abname)
		   ylabel('Rate, spikes/s')
		   title(titline)
		   if logabs~=0
			  axv = axis;
			  axis([xax(1), xax(2), axv(3), axv(4)]);
		   end
		   axlim = axis;
		   xleg = axlim(1)+0.05*(axlim(2)-axlim(1));
		   xdlt = 0.04*(axlim(2) - axlim(1));
		   if logabs~=0
			  xleg = axlim(1)*10.^(0.05*log10(axlim(2)/axlim(1)));
			  xdlt = xleg - axlim(1)*10.^(0.01*log10(axlim(2)/axlim(1)));
		   end
		   yleg = axlim(3) + 0.95*(axlim(4)-axlim(3));
		   ydlt = 0.05*(axlim(4)-axlim(3));
		   for j=1:nplots
			  text(xleg, yleg, legend(j));
			  line([xleg-xdlt, xleg], [yleg, yleg], 'Color', char(plcolor(j)));
			  yleg = yleg - ydlt;
		   end
		   if iscursrate~=0
			   text(xleg,yleg-ydlt,cursrate)
			   yleg = yleg - 2*ydlt;
		   end
			   
  
	  
		   what = input( ...
    '\nFn-filter, Hn-shift n dB right, Cursor, type 2 anal(@), Sa5-superimpose, Na5-next: ','s');
		   if what(1)=='F' | what(1)=='f'
			  if size(what,2) <= 1; what='f3'; end
			  nfw = sscanf(what(2:size(what,2)), '%g');
			  nfw = 2*floor(nfw/2) + 1;
			  legend(nplots) = strcat(legend(nplots), ', F', sprintf('%g',nfw));
			  drda(1:nlavg(nplots),nplots) = ...
					trifilt(drda(1:nlavg(nplots),nplots)',nfw)';
			  spda(1:nlavg(nplots),nplots) = ...
					trifilt(spda(1:nlavg(nplots),nplots)',nfw)';
			  if ispregnu~=0
				  legend(nplots-1) = strcat(legend(nplots-1), ', F', ...
				      sprintf('%g',nfw));
				  drda(1:nlavg(nplots-1),nplots-1) = ...
						trifilt(drda(1:nlavg(nplots-1),nplots-1)',nfw)';
				  spda(1:nlavg(nplots-1),nplots-1) = ...
						trifilt(spda(1:nlavg(nplots-1),nplots-1)',nfw)';
			  end
			  yleg = yleg + ydlt;
		   elseif what(1)=='H' | what(1)=='h'
			  nsh = round(sscanf(what(2:size(what,2)), '%g'));
			  uplo=input('Upper or Lower unit? ','s');
			  npl = nplots;
			  if uplo=='L' | uplo=='l';npl=nplots-1;end
			  abda(1:nlavg(npl),npl) = abda(1:nlavg(npl),npl)+nsh;
			  legend(npl) = cellstr(shleg(char(legend(npl)), nsh));
		   elseif what(1)=='N' | what(1)=='n'
			  str = what(2:size(what,2));
			  nplots = 0;
			  for j=1:10
				 legend(j) = cellstr('');
			  end
			  displayloop = 1;
		   elseif what(1)=='S' | what(1)=='s'
			  str = what(2:size(what,2));
			  displayloop = 1;
		   elseif what(1)=='C' | what(1)=='c'
			  fprintf('\n>>>> Position cursor, click mouse. <<<<\n')
			  [xval,rval]=cursor(gcf);
			  cursrate = sprintf('Rate = %g at %g on abscissa.',rval,xval);
			  iscursrate = 1;
		   elseif what(1)=='2' | what(1)=='@'
			  istypeII = 1;
			  uplo=input('Upper or Lower unit? ','s');
			  npl = nplots;
			  if uplo=='L' | uplo=='l';npl=nplots-1;end			  
			  if what(1)=='@'
				 [thrin,abmaxin,abminin]=showme;
			  else
			     thrin = 1000;abmaxin=1000;abminin=1000;
			  end
			  while istypeII==1
				  [nthr,spont,nd1,nd2,dynr,rmax,hislope,hioff,nhi1,nhi2] = ...
				      getII(drda(1:nlavg(npl),npl),abda(1:nlavg(npl),npl),...
					  thrin,abmaxin,abminin);
				  if nthr>0
					  hold on
					  axlim=axis;
					  line([abda(nthr,npl) abda(nthr,npl)], ...
					     [axlim(3) axlim(3)+0.2*(axlim(4)-axlim(3))],'Color','m');
					  line([abda(1,npl) abda(nthr,npl)],...
					     [spont spont],'Color','m');
					  plot([abda(nd1,npl) abda(nd2,npl)], ...
					       [drda(nd1,npl) drda(nd2,npl)], 'm*')
					  text(abda(nd1,npl)+0.02*(axlim(2)-axlim(1)), ...
					      drda(nd1,npl), ...
					      sprintf('dyn. range = %g dB',dynr),'Color','m')
					  if nhi1>0
						  hiy= hioff + hislope*abda(nhi1:nhi2,npl);
						  plot(abda(nhi1:nhi2,npl), hiy, 'm-')
						  text(abda(round((nhi1+nhi2)/2),npl)+0.02*(axlim(2)-axlim(1)), ...
						     drda(round((nhi1+nhi2)/2),npl), ...
							 sprintf('slope = %g /s*dB\n      = %g /dB', ...
							 hislope,hislope/rmax),'Color','m')
					  end
					  hold off
				  end
				  fprintf('Dyn range = %g; slope = %g /s*dB = %g /dB.\n', ...
				          dynr, hislope, hislope/rmax);
				  what=input('Print, Redo or Other processing? ','s');
				  if what=='p' | what=='P'
					 if ispropen==0
						fid = fopen(getlfn, 'at');
						fseek(fid, 0, 'eof');
						if ftell(fid)<=0;
						   fprintf(fid, ...
       'Date,trun     P.No.  dynrng   rmax    slope  norm_slope slope range\n');
					    end
					    ispropen = 1;
					 end
				     fprintf(fid, ...
					     '%s %s %5.1f %7.2f %8.4f %g %6.1f %6.1f\n', ...
					     nobl(titline),simplfy(char(legend(npl))),dynr, ...
						 rmax,hislope,hislope/rmax, ...
					     abda(nhi1,npl),abda(nhi2,npl));
				  elseif what=='r' | what=='R'
					  [thrin,abmaxin,abminin]=showme;
					  erasem(abda,drda,nthr,npl,spont,nd1,nd2,dynr, ...
							nhi1,nhi2,hioff,hislope);
				  else
					  istypeII=0;
				  end
			  end		  
		   else
			  displayloop = 1;
			  thatsall = 1;
		   end
		end       % #6  Display loop (while displayloop==0)
	else          % #5  No data in pile
	   thatsall = 1;
	end           % #5  Averaging and display
	end           % #2  Terminal end for nfils<=0
end               % #1  Loops to beginning from here until thatsall~=0

% On the way out, report the results

drat = [];srat=[];abval=[];nralvs=0;
if nplots>0
	drat = drda(:,1:nplots);
	srat = spda(:,1:nplots);
	abval = abda(:,1:nplots);
	npts = nlavg(1:nplots);
	nralvs = nplots;
end

return






function [fname, nfils, autofilt, filtw, issuper] = parscmrat(str)

% Parse command line string. This can take the form:
%   'a4-7,12,15-17[f3][s]'
% where
%    a4-7,12,15-17 means process P. Nos. 4-7, 12, 15-17
%    f3 means automatically filter with binwidth 3
%    s means superimpose, otherwise are averaged
% fname(1:nfils) return the file names as 'a4', 'a5', . . .
% autofilt~=0 if filtering desired, filter bw in filtw. NOTE filtw is returned
%      odd, so if 4 is entered, 5 is returned and used.
% issuper~=0 if rate files are to be superimposed; otherwise averaged

% fnames holds the filenames to be processed (up to 20)
	fname={'', '', '', '', '', '', '', '', '', '', ...
         	'', '', '', '', '', '', '', '', '', ''};
			

% autofilt = 1 if filtering specified in call string.
% issuper = 1 if superimposing files ('s' in call string).
	issuper = 0;
	autofilt = 0;
	filtw = 1;
	
% Parse out the elements of the command string
	n1 = 1;
	n2 = size(str,2);
	jfn = 1;
	for j=1:22
	   if n1>n2
			break
	   else
	      [A, nwds, ermsg, nch] = sscanf(str(n1:n2), '%c%g', 2);
	      n1 = n1 + nch -1;
	      if nwds<1 | (nwds<2 & (j==1 | char(A(1))=='-' | char(A(1))==','))
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
			   elseif char(A(1))=='f' | char(A(1))=='F'
				   autofilt = 1;
				   filtw = 3;
				   if nwds>=2; nfiltw = A(2); end
			   elseif char(A(1))=='s' | char(A(1))=='S'
				   issuper =1;
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







function [absc, logabs, delt, navrow, navskip, nplt, abname] = parsepb(pb)

% Computes abscissa and rate averaging parameters from the pb.
%    absc(1:nplt) is the abscissa data, after averaging
%    logabs=0 if linear abscissa
%          =1 if logarithmic
%    delt is the change from one abscissa point to the next. Equals the
%		   amount added (logabs=0) or the multiplier (logabs=1).
%    navrow is the number of lines of data to average into each plot point,
%         equal to the number of reps of each stimulus
%	 navskip is the number of rates to skip when averaging filrot stimuli
%    nplt is the number of plot points
%    abname is abscissa label name for rate plot

	logabs = 0;
	whrang = pb.mode_line(7);
	if whrang == 1
%  Frequency range in left (index=1) or right (index=6) channel
		navrow = max([1, pb.mode_line(1)]);
		navskip = 0;
		nplt = fix((99 + navrow)/navrow);
		abname = 'Frequency, kHz.';
		if pb.stim_par(1)==0 & pb.stim_par(2)==-1
			index = 1;
		elseif pb.stim_par(3)==0 & pb.stim_par(2)==-1
			index = 6;
	    elseif pb.stim_par(1)==7 & pb.ranges(1)~=0
			index = 1;
		else
			index = 6;
		end
		islog = pb.ranges(index+4);
		delt = pb.ranges(index+3);
		absc = zeros(nplt, 1);
		absc(1) = pb.ranges(index);
		if islog==0
			delt = delt/1000.;
			for j=2:nplt
				absc(j) = absc(j-1) + delt;
			end
		else
			logabs = 1;
			for j=2:nplt
				absc(j) = absc(j-1)*delt;
			end
		end
	elseif whrang == 2
%  Attenuation range in left (index=1) or right (index=6) channel
		navrow = max([1, pb.mode_line(1)]);
		navskip = 0;
		nplt = fix((99 + navrow)/navrow);
		abname = 'Attenuation, dB';
		if pb.stim_par(6) ~= -1
			fatt = pb.stim_par(5);
			latt = pb.stim_par(6);
		else
			fatt = pb.stim_par(7);
			latt = pb.stim_par(8);
		end
%     The very complicated attenuation delta from SELECT. delt is the change in
%     atten from line to line. NOTE that attenuations are made negative here.
		inpr = max([1, fix(100/pb.mode_line(1))]);
		delt = -sign(latt-fatt)*max([1,fix(abs((latt - fatt)/max([1,inpr-1]))+0.5)]);
	    absc = zeros(nplt, 1);
	    absc(1) = -fatt;
	    for j=2:nplt
		   absc(j) = absc(j-1) + delt;		
	    end
    elseif whrang == 5
%  File rotation
	    navrow = 0;
		abname = 'File number';
		if pb.stim_par(1) == 7
		   navskip = fix(0.5 + 100/pb.stim_par(2));
	    else
		   navskip = fix(0.5 + 100/pb.stim_par(4));
		end
		nplt = navskip;
		absc = zeros(nplt, 1);
		absc(1) = 1;
		delt = 1;
		for j=2:nplt
			absc(j) = absc(j-1) + delt;
		end
	else
%  No range or a range not understood
		navrow = 1;
		navskip = 0;
		nplt = 100;
		abname = 'Line number';
		absc = ones(100, 1);
		delt = 1;
		for j=2:100
			absc(j) = absc(j-1) + delt;
		end
	end
return
		
		
		

		
		
function [rate, srate, linf, lstlin] = comprate(resol, timew);

% Comute rate vs. line number for the file opened with labfile(). Contents
% of the file read with fread(), file closed on exit.

global FILEPNTRXYZ

% Convert time windows to clock cycles. Driven spikes are in [win1, win2],
% spont spikes are in [win3, win4].
	win1 = 1000.*timew(1)/resol;
	win2 = 1000.*timew(2)/resol;
	win3 = 1000.*timew(3)/resol;
	win4 = 1000.*timew(4)/resol;
	conv = resol/1.E6;              % Clock tick in seconds
	
	lstlin = 0;
	rate = zeros(100, 1);
	srate = zeros(100, 1);
	linf = zeros(100, 1);
	
	line = 0;
	in = 0;

	while in ~= -16384
		in = fread(FILEPNTRXYZ, 1, 'int16');
		if in == -32768
			cum = cum + 32768;
		elseif in < 0
			in1 = fread(FILEPNTRXYZ, 1, 'int16');
			if line > 0
			   lineleng = cum + (in1 + 32768);
			   durd = conv*(min([win2, lineleng]) - win1);
			   durs = conv*(min([win4, lineleng]) - win3);
			   if durd > 0
				  rate(line) = count/durd;
				  linf(line) = -1;
			   end
			   if durs > 0
				  srate(line) = scount/durs;
				  linf(line) = 1;
			   end
		    end
		   	if in~=-16384
			   line = in + 32768;
			   lstlin = max([line, lstlin]);
			end
			cum = 0.;
			count = 0;
			scount = 0;
		else
			if line <= 0 
				fprintf(1,'***ERROR spikes before a line-start.***')
				in = -16384;
			else
				lat = cum + in;
				if lat >= win1 & lat <= win2
					count = count + 1;
				end
				if lat >= win3 & lat <= win4
					scount = scount + 1;
				end
			end
		end
	end
	fclose(FILEPNTRXYZ);
return
		
		
		
		
function [rate, srate, llinf, nlavg] = avgrate(rrate, rsrate, linf, lstlin, navrow, navskip, nplt);

% Takes raw rates from comprate() and averages together lines with the same
% stimulus value, as determined in parsepb(). rrate, rsrate, linf are rates and
% line-found vectors from comprate. lstlin is the largest line number found in
% comprate. nplt is the theoretical number of different stimuli in the picture,
% equal to the theoretical number of lines at the output of this routine. nlavg
% is the actual number of lines at the output of this routine (smaller than
% nplt because lstlin<100 and some stimuli are missing).
% NOTE:RATE/RRATE, SRATE/RSRATE, AND LINF/LLINF MUST BE SEPARATE VECTORS.

	rate = zeros(nplt, 1);
	srate = zeros(nplt, 1);
	llinf = zeros(nplt, 1);
	nlavg = 0;
	
	if navrow>= 2
		j1 = 1;
		for j=1:nplt
			nd = 0;
			ns = 0;
			for j2=1:navrow
				if linf(j1)~=0
					rate(j) = rate(j) + rrate(j1);
					nd = nd + 1;
				end
				if linf(j1)>0
					srate(j) = srate(j) + rsrate(j1);
					ns = ns + 1;
				end
				j1 = j1 + 1;
				if j1>lstlin 
					break
				end
			end
			if nd>0 
				rate(j) = rate(j)/nd;
				nlavg = max([nlavg, j]);
				llinf(j) = -1;
			end
			if ns>0
				srate(j) = srate(j)/ns;
				llinf(j) = 1;
			end
			if j1>lstlin
				break
			end
		end
	elseif navskip > 0
		for j=1:nplt
			nd = 0;
			ns = 0;
			for j1=j:navskip:lstlin
				if linf(j1)~=0
					rate(j) = rate(j) + rrate(j1);
					nd = nd + 1;
				end
				if linf(j1)>0
					srate(j) = srate(j) + rsrate(j1);
					ns = ns + 1;
				end
			end
			if nd>0 
				rate(j) = rate(j)/nd;
				nlavg = max([nlavg, j]);
				llinf(j) = -1;
			end
			if ns>0
				srate(j) = srate(j)/ns;
				llinf(j) = 1;
			end
		 end
	else
		nlavg = min([lstlin, nplt]);
		rate(1:nlavg) = rrate(1:nlavg);
		srate(1:nlavg) = rsrate(1:nlavg);
		llinf(1:nlavg) = linf(1:nlavg);
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
		



function legout = shleg(legin, nsh)

% Correct legend report of shift. (add ',sh-4' or ',sh5' to end of legend).
% First check if a shift report is already present

	nshc = nsh;
	indsh = findstr(',sh', legin);
	if ~isempty(indsh)
		n1 = indsh;
		n2 = size(legin,2)+1;
		indsh = findstr(',', legin(n1+3:n2-1));
		if ~isempty(indsh); n2=n1+2+indsh; end
		nshc = round( nsh + sscanf(legin(n1+3:n2-1),'%i',1));
		if n2<=size(legin,2)
			holdleg = strcat(legin(1:n1-1),legin(n2:size(legin,2)));
		else
			holdleg = legin(1:n1-1);
	    end
	else
		holdleg = legin;
	end
	if nshc~=0
	   legout = strcat(holdleg, sprintf(',sh%i', nshc));
    else
	   legout = holdleg;
	end
return





function [arate, asrate, aabsc, nlevs, aleg] = avgplots(drda, spda, ...
       abda, nlavg, delts, logs, legend, nplthis, nplots)

% Average together plots in data pile locations [nplthis, nplots] to form
% a single plot. Averaging can be done iff
%   All are log or all are linear
%   All have same delt
% Points with same abscissa are averaged and the resulting abscissa is the
% portion that overlaps in all input plots only.
% If averaging fails, nothing is returned and nlevs = 0.

	canavg = 1;
	for j=nplthis+1:nplots
		if logs(j)~=logs(nplthis); canavg=0; end
		if delts(j)~=delts(nplthis); canavg=0; end
	end

% Summate input files

	if canavg==1
		arate = drda(1:nlavg(nplthis),nplthis);
		asrate = spda(1:nlavg(nplthis),nplthis);
		aabsc = abda(1:nlavg(nplthis),nplthis);
		nlevs = nlavg(nplthis);
		[bleg, endleg] = avleg1(char(legend(nplthis)));
		for j=nplthis+1:nplots
		   stavg = 1;
		   stnew = 1;
		   endavg = min([nlevs, nlavg(j)]);
		   endnew = endavg;
		   if abda(1,j)~=aabsc(1)
			  shnew = round((abda(1,j)-aabsc(1))/delts(nplthis));
			  if shnew>0
				  stavg = 1 + shnew;
				  endavg = min([nlevs, nlavg(j)+shnew]);
				  endnew = min([nlavg(j), nlevs-shnew]);
			  else
				  stnew = 1 - shnew;
				  endnew = min([nlavg(j), nlevs-shnew]);
				  endavg = min([nlevs, nlavg(j)+shnew]);
			  end
		   end
		   arate = arate(stavg:endavg) + drda(stnew:endnew, j);
		   asrate = asrate(stavg:endavg) + spda(stnew:endnew, j);
		   aabsc = aabsc(stavg:endavg);
		   nlevs = endavg - stavg + 1;
		   bleg = avleg2(char(legend(j)), bleg);
	    end
	    navgd = nplots - nplthis + 1;
	    arate = arate/navgd;
		asrate = asrate/navgd;
		aleg = cellstr(strcat(bleg,endleg));
	else
	   nlevs = 0;
    end
return
 



function [aleg, endleg] = avleg1(legin)

% Set up legend for files averaged. Legend contains comma-delimited list
% of file names (in aleg) terminated by qualifiers like ', F3' (in endleg).
% avleg2 adds the file names.

	indend = findstr(', F', legin);
	if isempty(indend)==0
		endleg = legin(indend:size(legin,2));
	else
		endleg = ' ';
	end
	
	indfn = findstr(',', legin);
	if isempty(indfn)~=0
		indfn = size(legin,2);
	end
	aleg = strcat('Files', legin(5:indfn-1));
return


function aleg = avleg2(legin, bleg)

	indfn = findstr(',', legin);
	if isempty(indfn)~=0
		indfn = size(legin,2);
	end
	aleg = strcat(bleg, ',', legin(5:indfn-1));
return





% String functions for printout only

function strout=nobl(strin)
%Output string is same as input with no blanks
	strout = strin;
	j1=0;
	for j=1:max(size(strin))
		if strin(j)~=' ';j1=j1+1;strout(j1)=strin(j);end
	end
	strout=strout(1:j1);
return

function strout=simplfy(strin)
% Simplify legend strings to contain only picture identifications

	indx=findstr(', S',strin);
	if isempty(indx)~=0; indx=findstr(', L',strin);end
	if isempty(indx)==0
		strout = nobl(strin(5:indx(1)-1));
		indx = findstr(', F',strin);
		if isempty(indx)==0
			strout = nobl(strcat(strout,strin(indx(1):size(strin,2))));
		end
	else
	   strout = nobl(strin(6:size(strin,2)));
    end
return



function [thrin,abmaxin,abminin]=showme
	fprintf('>>> Show me the threshold, click. <<<\n')
	[thrin,junk]=cursor(gcf);
	fprintf('>>> Show me the inflection, click. <<<\n')
	[abmaxin,junk]=cursor(gcf);
	fprintf('>>> Show me the right side of hislope, click. <<<\n')
	[abminin,junk]=cursor(gcf);
return


function erasem(abda,drda,nthr,npl,spont,nd1,nd2,dynr, ...
    nhi1,nhi2,hioff,hislope)
% Erase the lines (sort of)
	hold on
	axlim=axis;
	line([abda(nthr,npl) abda(nthr,npl)], ...
	   [axlim(3) axlim(3)+0.2*(axlim(4)-axlim(3))],'Color','w');
	line([abda(1,npl) abda(nthr,npl)],...
	   [spont spont],'Color','w');
	plot([abda(nd1,npl) abda(nd2,npl)], ...
	   [drda(nd1,npl) drda(nd2,npl)], 'w*')
	text(abda(nd1,npl)+0.02*(axlim(2)-axlim(1)), drda(nd1,npl), ...
	   sprintf('dyn. range = %g',dynr),'Color','w')
	if nhi1>0
	   hiy= hioff + hislope*abda(nhi1:nhi2,npl);
	   plot(abda(nhi1:nhi2,npl), hiy, 'w-')
	   text(abda(round((nhi1+nhi2)/2),npl)+0.02*(axlim(2)-axlim(1)), ...
	      drda(round((nhi1+nhi2)/2),npl), ...
		  sprintf('slope = %g',hislope),'Color','w')
	end
	hold off
return



function lfnam = getlfn
% Gets logfile name of form 'Hostname:Desktop Folder:type_II_log.txt'.
% Uses matlabroot to discover the hostname.

	str = matlabroot;
	ind = find(':'==str);
	if isempty(str)~=0
		lfnam = str;
	else
		lfnam = str(1:ind(1)-1);
	end
	lfnam = strcat(lfnam,':Desktop Folder:type_II_log.txt');
return

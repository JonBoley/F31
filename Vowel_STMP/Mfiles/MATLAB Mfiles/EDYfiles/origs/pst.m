function [psth_dat,time,stim_string]=pst(fn)
% PSTH maker:      PK 12/20/98
%   'a4-7,12,15-17b2f7t400z50-51z100-102'
% where
%    a4-7,12,15-17 means process P. Nos. 4-7, 12, 15-17
%    b2 means use a binwidth of 2 ms (optional)
%    f7 means filter with a 7-bin wide triangular filter (optional)
%	 t400 means force PST abscissa to be 400 ms.
%	 z50-51 means zero from 50 ms to 51 ms
%    q1  quiet mode, no output plot 


[fname, nfils, binset, binwset, autofilt, ...
     filtw, lengset, pstleng,zerocnt,zeroval,qmode] = parscmpst(fn);
for i=1:nfils
	[pst_dat,time,stimu]=psth_pk (fname(i,:),binwset,0,pstleng,zerocnt,zeroval);
	pst_dat(1)=pst_dat(2);
	if filtw>1
		stimu=sprintf('%s  F%i',stimu,filtw);
	end
	if i==1
		psth_dat=pst_dat;
		stim_string=stimu;
	else
		psth_dat=psth_dat+pst_dat;
		stim_string=sprintf('%s\n%s',stim_string,stimu);
	end
		
end
psth_dat=psth_dat/nfils;
psth_dat=trifilt(psth_dat, filtw);
if qmode
	return
end
pts=length(psth_dat);
plot(time(2:pts),psth_dat(2:pts))
axis([0,time(pts),0,max(psth_dat)*1.2])
text(10,max(psth_dat)*1.1,sprintf('%s',stim_string),'FontSize',9);
xlabel('ms')
ylabel('spikes/s')
stimu=findstr(stim_string,':');
title(sprintf('Unit %s',stim_string(1:stimu(1)-1)))


return


function [fname, nfils, binset, binwset, autofilt, ...
     filtw, lengset, pstleng,zeroing,zeroval,qmode] = parscmpst(str);
% taken from EDY and changed a bit.....
%
% [fname,nfils,binset,binwset,autofilt,filtw] = PARSCOMM(str)
% Parse command line string. This can take the form:
%   'a4-7,12,15-17b2f7t400z50z100'
% where
%    a4-7,12,15-17 means process P. Nos. 4-7, 12, 15-17
%    b2 means use a binwidth of 2 ms (optional)
%    f7 means filter with a 7-bin wide triangular filter (optional)
%	 t400 means force PST abscissa to be 400 ms.
%	 z50z51 means zero from 50 ms to 51 ms
%	d1 for dot raster output
% fname(1:nfils) return the file names as 'a4', 'a5', . . .
% binset~=0 if binwidth specified in binwset
% autofilt~=0 if filtering desired, filter bw in filtw. NOTE filtw is returned
% odd, so if 4 is entered, 5 is returned and used.


% binset = 1 if binwidth set in call string, otherwise determined from file
% autofilt = 1 if filtering specified in call string.
	qmode=0;
	dotflag=0;
	binset = 0;
	autofilt = 0;
	binwset = 0.1;
	filtw = 1;
	lengset = 0;
	pstleng = 600;
	zeroing=0;
	zerocnt =1;
	zeroval = zeros(1,10);
	extension='.DAT';
	
	
% Parse out the elements of the command string
	n1 = 1;
	n2 = size(str,2);  % nr of columns
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
				 
				filn=sprintf('%c%g%s',A(1),A(2),extension);
			   	fname(jfn,1:length(filn)) = filn(1:length(filn));
			   	jfn = jfn + 1;
			   	B = A(1);
			   	C = A(2)+1;
		     else
			  	if char(A(1))=='-'
				   for j=C:A(2)
					   filn=sprintf('%c%g%s',B,j,extension);
					   fname(jfn,1:length(filn)) = filn(1:length(filn));
				       jfn = jfn + 1;
				   end
			   elseif char(A(1))==','
				   filn=sprintf('%c%g%s',B,A(2),extension);
				   fname(jfn,1:length(filn)) = filn(1:length(filn));
				   jfn = jfn + 1;
				   C = A(2)+1;
			   elseif char(A(1))=='b' | char(A(1))=='B'
				   %binwidth
				   binset = 1;
				   binwset = A(2);
			   elseif char(A(1))=='f' | char(A(1))=='F'
				   %filter
				   autofilt = 1;
				   filtw = A(2);
			   elseif char(A(1))=='t' | char(A(1))=='T'
				   %max time
				   lengset = 1;
				   pstleng = A(2);
			   elseif char(A(1))=='q' | char(A(1))=='Q'
				   %max time
				   qmode = 1;
			   elseif char(A(1))=='z' | char(A(1))=='Z'
		   			%zero artifacts
					zeroval(zerocnt) = A(2);
				   	zerocnt = zerocnt + 1;
					if n1>=n2
					  	zeroval(zerocnt) = zeroval(zerocnt-1);
				   		zerocnt = zerocnt + 1;	
					else
						[A, nwds, ermsg, nch] = sscanf(str(n1:n2), '%c%g', 2);
	      				n1 = n1 + nch -1;
			  			if char(A(1))=='-'
				  		 	zeroval(zerocnt) = A(2);
				   			zerocnt = zerocnt + 1;		
						else
						  	zeroval(zerocnt) = zeroval(zerocnt-1);
				   			zerocnt = zerocnt + 1;	
						end
					end
			   end
		     end
	      end
	   end
   end
   nfils = jfn - 1;
   if autofilt~=0
	   filtw = max([1, 1+2*fix(filtw/2)])
   end
	zeroing=zerocnt-1;
return

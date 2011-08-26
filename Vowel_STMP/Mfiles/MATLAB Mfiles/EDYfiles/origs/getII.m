function [nthr,spont,nlocr,nhicr,dynr,rmax,hislope,hioff,nhi1,nhi2] = ...
         getII(drate,absc,thrin,abmxin,abmnin);

% Function computes dynamic range and high frequency slope for type II and similar
% units. Proceeds by:
% 1. Finds threshold as first point 3 sd above average rate
% 2. Finds max value within 50 dB of threshold. Computes dynamic range as 10-90%
%    from threshold to max.
% 3. Finds minimum at least 20% below max at levels above max and computes slope
%    from max to minimum.
% Inputs:
%	thrin - forces threshold to a particular point, unless 1000
%   abmxin - forces inflection point to a particular point, unless 1000
%   abmnin - forces right end of hi level slope region to part. point, unless 1000
% Returns:
%   nthr - absc index at which threshold found (absc(nthr),drate(nthr))
%         returns 0 if no threshold found.
%	spont - average spont when threshold found
%	dynr - 10-90% dynamic range in abscissa units
%	Line fitting rates above rate maximum: rate = hioff + hislope*absc for
%	      absc(nhi1:nhi2).

% Determine abscissa increment

	sab=size(absc,1);
	ainc = round(mean(absc(2:sab)-absc(1:sab-1)));
	n5dB = round(5/ainc);

% Are there forced inputs?
	nthr=0;indmx=0;indmn=0;
	if thrin~=1000; nthr = findab(thrin,absc,ainc); end
	if abmxin~=1000; indmx = findab(abmxin,absc,ainc); end
	if abmnin~=1000; indmn = findab(abmnin,absc,ainc); end

% Find first rate that exceeds 10 dB from start and 2 SDs from spont

	if nthr==0
		for j=n5dB:sab-1
			if drate(j+1) > mean(drate(1:j))+3.*std(drate(1:j))
				nthr=j+1; thr = absc(nthr); spont=mean(drate(1:j)); break;
			end
		end
	else
		thr = absc(nthr);
		spont = mean(drate(1:nthr-1));
	end

	if nthr<=0
		fprintf('*** getII: no threshold in rate data. ***\n')
		spont=mean(drate);dynr = 0; hislope=0;
	else
		
% Find next rate maximum, within 50 dB, then find dynamic range points (10,90%)
% searching outward from the center of the dynamic range.

	   if indmx==0
		   n50dB = min([nthr + round(50/ainc), sab]);
		   [rmax,indmx]=max(drate(nthr:n50dB));
		   indmx = indmx + nthr -1;
	   else
		   rmax = drate(indmx);
	   end
	   rlocr = spont + 0.1*(rmax-spont);
	   rhicr = spont + 0.9*(rmax-spont);
	   nlocr = 0;
	   nhicr = 0;
	   rmid = (rmax+spont)/2;
	   for j=nthr:indmx
		   if drate(j)>rmid;jmid=j;break;end
	   end
	   for j=1:max([indmx-jmid, jmid])
		   if nlocr<=0 & jmid-j>0 & drate(jmid-j)<=rlocr;nlocr=jmid-j;end
		   if nhicr<=0 & jmid+j<=indmx & drate(jmid+j)>=rhicr;nhicr=jmid+j;end
	   end
	   dynr = (nhicr-nlocr)*ainc;

% Find rate minimum above rmax and fits a line to rates from rmax to rmin.
% This checks several things:
%   1. Must be at least 10 data points, otherwise takes all points.
%   2. Minimum is accepted only if it's 30% below rmax

	   nhi1=0;nhi2=0;hislope=0;hioff=0;
	   if sab-indmx>=10 | indmn~=0
		   if indmn==0
			   [rmin,indmn]=min(drate(indmx:sab));
			   indmn = indmn + indmx -1;
			   if rmin>0.8*rmax | indmn-indmx<9;indmn=sab;end
		   else
			   rmin = drate(indmn);
		   end
		   nhi1 = indmx; nhi2=indmn;
		   mabs = mean(absc(nhi1:nhi2));
		   hislope = ((absc(nhi1:nhi2)-mabs)'*drate(nhi1:nhi2))/ ...
		        ((nhi2-nhi1)*std(absc(nhi1:nhi2))^2);
		   hioff = mean(drate(nhi1:nhi2))-hislope*mabs;
	   end
	end
return


function ind = findab(val,data,res)

% Get index of val in data with resolution res.

	inar = data-0.5*res<=val & data+0.5*res>val;
	ind = find(inar);
	
return

function [Q10, frhi, frlo, Q10lev] = findQ10(f,inten,BF)
% File: findQ10.m
% Modified by M. Heinz 8/1/02 from btq2.m
%
% Now BF is passed, because Q10 is specified as 10 dB above BFthreshold,
% and (impaired) BFs may be chosen to be different than the minimum on the TC (e.g., based on Liberman and Dodds, 1984)
%     also, interpolation changed from linear to log on y-axis
%
%function used to find Q10 of TC, requires colomn vectors 'frequency' and 'intensity'
%which uses a slightly different method and takes into account cases in which either or both
%side of TC end before reaching 10dB above threshold.
%original program by Elizateth Allegretto (1998), modified by Teng Ji (1999)
%function [BF, thresh, Q10, frhi, frlo] = btq(f, inten)


%find the threshold based on passed BF (must be a value in f)
temp5 = size(f);
j10=find(f==BF);  %index of BF
if isempty(j10)
   error('BF passed to FindQ10 was not in frequency array')
end
thresh = inten(j10);  % This is based on passed inten (i.e, smoothed TC, so Q10lev~=unit.thresh+10)

%find Q10
Q10lev = thresh + 10;
frhia = f(1);
frhib = f(1);
frloa = f(temp5(1));
frlob = f(temp5(1));
inhia = inten(1);
inhib = inten(1);
inloa = inten(temp5(1));
inlob = inten(temp5(1));

frlo = 10000;
frhi = 20000;

%for high frequency portion of TC
if inten(1) < Q10lev %in case not enough data pts were collected 
   %    slp=(inten(1)-inten(2))/(f(1)-f(2));
   %    frhi=(Q10lev-inten(1))/slp+f(1);
   slp=(inten(1)-inten(2))/(log10(f(1))-log10(f(2)));
   frhi=10^((Q10lev-inten(1))/slp+log10(f(1)));
else %finding pts immediately above(frhib) and below(frhia) Q10lev   
   for j = 2:j10
      if frhia == f(1) 
         if inten(j) <= Q10lev
            frhia = f(j);
            inhia = inten(j);
         end
      end   
      if (inten(j) >= Q10lev) & (frhia == f(1))
         frhib = f(j);
         inhib = inten(j);
      end
   end
   %  slope = (inhib - inhia)/(frhib - frhia);
   % 	frhi=(Q10lev-inhia)/slope+frhia;
   slope = (inhib - inhia)/(log10(frhib) - log10(frhia));
   frhi=10^((Q10lev-inhia)/slope+log10(frhia));
end

%for low frequency portion of TC
if inten(temp5(1))< Q10lev %in case not enough data pts were collected
   %    slp=(inten(temp5(1))-inten(temp5(1)-1))/(f(temp5(1))-f(temp5(1)-1));
   %    frlo=(Q10lev-inten(temp5(1)))/slp+f(temp5(1));
   slp=(inten(temp5(1))-inten(temp5(1)-1))/(log10(f(temp5(1)))-log10(f(temp5(1)-1)));
   frlo=10^((Q10lev-inten(temp5(1)))/slp+log10(f(temp5(1))));
else %finding pts immediately above(frlob) and below(frloa) Q10lev
   for j = temp5(1):-1:j10
      if frloa == f(temp5(1))
         if inten(j) <= Q10lev
            frloa = f(j);
            inloa = inten(j);
         end
      end   
      if (inten(j) >= Q10lev) & (frloa == f(temp5(1)))
         frlob = f(j);
         inlob = inten(j);
      end
   end
   %  slope = (inlob - inloa)/(frlob - frloa);
   % 	frlo=(Q10lev-inloa)/slope+frloa;
   slope = (inlob - inloa)/(log10(frlob) - log10(frloa));
   frlo=10^((Q10lev-inloa)/slope+log10(frloa));
end

temp10 = [0.001 (frhi-frlo)];
Q10 = BF/max(temp10);   
errors=0;

if Q10 == BF/0.001;
   'Error in Q10, frhi-frlo < 0.001'
   if ~exist('j','var'),  j=99;   end
   figure(j)
   plot (f(:,1),inten(:,1));
end   

if frhi == 20000
   'Error in Q10, frhi not defined'
end

if frlo == 10000
   'Error in Q10, frlo not defined'
end   
         
return

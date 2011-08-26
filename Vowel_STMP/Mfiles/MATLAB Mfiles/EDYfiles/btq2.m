%function used to find Q10 of TC, requires colomn vectors 'frequency' and 'intensity'
%which uses a slightly different method and takes into account cases in which either or both
%side of TC end before reaching 10dB above threshold.
%original program by Elizateth Allegretto (1998), modified by Teng Ji (1999)

function [BF, thresh, Q10, frhi, frlo] = btq(f, inten)

%find the BF  and theshold
temp5 = size(f);
thresh = min(inten);
for j = 1:temp5(1)
   if inten(j) == thresh
      BF = f(j);
      j10=j;
   end      
end   

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
   slp=(inten(1)-inten(2))/(f(1)-f(2));
   frhi=(Q10lev-inten(1))/slp+f(1);
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
   slope = (inhib - inhia)/(frhib - frhia);
	frhi=(Q10lev-inhia)/slope+frhia;
end

%for low frequency portion of TC
if inten(temp5(1))< Q10lev %in case not enough data pts were collected
   slp=(inten(temp5(1))-inten(temp5(1)-1))/(f(temp5(1))-f(temp5(1)-1));
   frlo=(Q10lev-inten(temp5(1)))/slp+f(temp5(1));
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
   slope = (inlob - inloa)/(frlob - frloa);
	frlo=(Q10lev-inloa)/slope+frloa;
end

temp10 = [0.001 (frhi-frlo)];
Q10 = BF/max(temp10);   
errors=0;

if Q10 == BF/0.001;
   'Error in Q10, frhi-frlo < 0.001'
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

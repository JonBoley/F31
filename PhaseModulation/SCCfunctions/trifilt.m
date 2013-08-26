function otvect = trifilt(invect, nfw, ignoreNaNs)
% Modified by M. Heinz 04Jun2004 
% Added an option (ignoreNaNs=1) to ignore NaNs within the window if not at center point
% e.g., helps to keep from losing dynamic range on rate-level curves with NaNs at HLs
%
% TRIFILT filters ROW-vector invect with triangular filter of width nfw.
%          otvect = trifilt(invect,nfw)
% nfw must be odd and will behave as if it is nfw+1 if it is even.
% invect is convolved with a zero-centered smoothing function as:
%   (e.g. for nfw=5)  0.11  0.22  0.33  0.22  0.11
%   for n=             -2    -1     0     1     2
% invect is padded with its first and last value in both directions
% but only the output points corresponding to the original input
% points are returned

if isempty(invect)  % Catch empty invect passed
   otvect=[];
   return;
end

if ~exist('ignoreNaNs','var')
   ignoreNaNs=0;
end

%%%%%%%%%%%%%%% 6/4/04 TODO MGH
% uses conv, so hard to fix NaNs as it goes.
% So, need to (if ignoreNaNs), go back at end and re-compute by hand those values in otvect that came up as NaN
%
% May take extra time, but this is only when ignoreNaN=1
%


nfwi = 2*floor(nfw/2) + 1;
filt = zeros(1,nfwi);
summ = 0;
nfw2 = floor(nfwi/2);
for j=1:nfw2
   filt(j) = j;
   filt(nfwi+1-j) = j;
   summ = summ + 2*j;
end
nfw3 = nfw2 + 1;
filt(nfw3) = nfw3;
summ = summ + nfw3;
for j=1:nfwi
   filt(j) = filt(j)/summ;
end
svect = size(invect,2) + 2*nfw2;
vect1 = zeros(1,svect);
vect1(1:nfw2) = invect(1)*ones(1,nfw2);
vect1(nfw3:svect-nfw2) = invect;
vect1(svect-nfw2+1:svect) = invect(size(invect,2))*ones(1,nfw2);
vect2 = conv(vect1, filt);
otvect = vect2(2*nfw2+1:svect);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% if option is set, Recompute filtered values, ignoring incoming NaNs when not at the center point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ignoreNaNs
   Ninds=find(isnan(otvect)&~isnan(invect));
   filtORIG=filt*summ;
   for i=Ninds
      invectTEMP=vect1(i:i+2*nfw2);
      GOODinds=find(~isnan(invectTEMP));
      filtTEMP=filtORIG(GOODinds)/sum(filtORIG(GOODinds));
      otvect(i)=invectTEMP(GOODinds)*filtTEMP';
   end
end

return;

	function otvect = trifilt(invect, nfw)
% TRIFILT filters ROW-vector invect with triangular filter of width nfw.
%          otvect = trifilt(invect,nfw)
% nfw must be odd and will behave as if it is nfw+1 if it is even.
% invect is convolved with a zero-centered smoothing function as:
%   (e.g. for nfw=5)  0.11  0.22  0.33  0.22  0.11
%   for n=             -2    -1     0     1     2
% invect is padded with its first and last value in both directions
% but only the output points corresponding to the original input
% points are returned
	nfwi = 2*floor(nfw/2) + 1;
	filt = zeros(1,nfwi);
	sum = 0;
	nfw2 = floor(nfwi/2);
	for j=1:nfw2
	   filt(j) = j;
	   filt(nfwi+1-j) = j;
	    sum = sum + 2*j;
	end
	nfw3 = nfw2 + 1;
	filt(nfw3) = nfw3;
	sum = sum + nfw3;
	for j=1:nfwi
	   filt(j) = filt(j)/sum;
	end
	svect = size(invect,2) + 2*nfw2;
	vect1 = zeros(1,svect);
	vect1(1:nfw2) = invect(1)*ones(1,nfw2);
	vect1(nfw3:svect-nfw2) = invect;
	vect1(svect-nfw2+1:svect) = invect(size(invect,2))*ones(1,nfw2);
	vect1 = conv(vect1, filt);
	otvect = vect1(2*nfw2+1:svect);
return

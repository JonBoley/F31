function flnum = pffc(ipfn)

% Takes pseudofloating number ipfn, from lab-system parameter block,
% and returns a floating number flnum.
% Pseudofloating numbers are 16-bit strings of the form
%      ee mmmmmmmmmmmmmm
%      12 34567890123456
% and are interpreted as
%    number = mmmmmmmmmmmmmm / 10^(2+ee)
% Thus if mmmmmmmmmmmmmm = 8769 and ee = 1, number is 8.769

	if ipfn>0
		binstr = dec2bin(ipfn);
	else
		binstr = dec2bin(65536+ipfn);
	end
	n1 = size(binstr,2);
	n2 = max([1, n1-13]);
	mant = bin2dec(binstr(n2:n1));
	exp = 0;
	if n1 > 14
		n3 = max([1,n1-14]);
		n4 = max([1,n1-15]);
		exp = bin2dec(binstr(n4:n3));
	end
	flnum = mant/10.^(2.+exp);
	return

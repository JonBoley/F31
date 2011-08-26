function srate = getsrate(str)

% srate = GETSRATE(str) returns sampling rate for TDT filrot stimuli from
% lab parameter specification line. e.g.
% gtsrate('S.XR4W(152.), . . .')       returns 152

	n1 = 1;
	n2 = size(str, 2);
	A = 0;
	while char(A)~='('
		[A,nwd,errm,nch]=sscanf(str(n1:n2),'%c',1);
		n1 = n1 + nch - 1;
	end
	srate = sscanf(str(n1:n2), '%g', 1);
return

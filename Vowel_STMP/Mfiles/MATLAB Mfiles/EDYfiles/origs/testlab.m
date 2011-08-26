% Program to test labfile and nxtword.

global BUFFXYZ BUFFPNTRXYZ FILEPNTRXYZ TOPOFBUFFXYZ FILECOUNTXYZ MAXCOUNTXYZ INPREVXYZ

	fname = input('Gimme a name (a47): ', 's');
	pb = labfile(fname);
	if pb.Error(1:6)~='No err'
		fprintf('*** Error = ''%s'' while opening file, check directory.***\n',pb.Error);
	else
	j = 0;
	in = 0;
	while in ~= -16384
		in = fread(FILEPNTRXYZ, 1, 'int16');
%		in = nxtword;
		if in<0 & in~=-32768
			fprintf(1,'\nLine %o',65536+in);
			in1 =  fread(FILEPNTRXYZ, 1, 'int16');
			in1 = nxtword;
			j = j+1;
		end
		j = j + 1;
	end
	fprintf('\nThere are %d numbers in file %s.\n', j, fname)
	fclose('all');
	end

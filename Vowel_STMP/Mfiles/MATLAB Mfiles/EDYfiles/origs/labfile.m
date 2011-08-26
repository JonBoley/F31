function pb = labfile(filename)

% LABFILE: When called as pb = labfile('a47'), opens file a47.DAT, returns parameter
% block in pb. File is left open so fread(FILEPNTRXYZ, nn, 'int16') can be used
% to return data. On file errors FILEPNTRXYZ is set to -1 and the pb consists of
% only pb.Error = 'File error'. If the parameter block is read correctly, then
% pb.Error = 'No error' (see read_parm()).

global FILEPNTRXYZ

	endian='ieee-le'; %data from little endian => integers flipped

% Open data file
	fname = sprintf('%s.DAT',filename);
	FILEPNTRXYZ = fopen(fname,'r',endian);
	
	if FILEPNTRXYZ <= 0
		fprintf(1,'\n***Error opening %s, probably wrong directory.***\n', fname)
		pb = struct('Error', 'File error');
	else
		
% Get parameter block.
	    pb = read_parm1(FILEPNTRXYZ);
	end
return

	

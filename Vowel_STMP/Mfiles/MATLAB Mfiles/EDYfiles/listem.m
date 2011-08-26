% LISTEM - opens files from lab data, writes text file on the desktop
% containing stim spec lines for requested range of files.

global FILEPNTRXYZ

	fc = input('Enter file ID character ','s');
	jf1 = input('First P. No. (required) ');
	jf2 = input('Last P. No. (=0 if all) ');
	if jf2==0; jf2n=99999;else;jf2n=jf2;end
	
	nlines = 0;
	for j=jf1:jf2n
		fprintf('Opening file %c%g . . . ',fc,j);
		pb = labfile( sprintf('%c%g', fc, j));
		if strcmp('File error', pb.Error)==1
			break
		end
		fclose(FILEPNTRXYZ);
		if j==jf1
         	fid = fopen(strcat('Elburz:Desktop Folder:D',pb.date(1:8)), 'wt')

		end
		fprintf('& writing.\n');
		fprintf(fid, '%s, %s: %s\n', pb.date, pb.unit, pb.ltext);
		if pb.short==0
			fprintf(fid, '                          %s\n',pb.mtext);
			fprintf(fid, '                          %s\n',pb.rtext);
		end
		nlines = nlines + 1;
	end
	if nlines>=1;fclose(fid); end
	
	

function [parameterblock] = read_parm(filename)

% reads and decodes NEL datafile param block and puts them into a struct
% 	pb = read_parm('a3')    <- reads a3.DAT into structure pb.
% Complete contents of structure documented in the m file.
% PK 10/2/97, EDY 12/8/97
% NOTE: this m-file opens and closes the data file. If file handling is
% done externally, use read_parm1(fptr)

% Contents of the output structure:
%           date: '12/04/97'          date of expt.
%           time: '18:10:21'          start time of picture
%           unit: '1.03  , P20   '    track,unit no's, pict. no.
%      prog_name: 'ST'                program used to take data
%     file_descr: [ 4x1  double]      no. 512-byte blocks in the
%										fil;, pict. no.; 1; data
%										format: 1-dot disp, 4-TC,
%										19-AVELN, 30-SRCAL, etc.
%          short: 1                   1 if S. only; 0 if L.,M.,R.
%          ltext: 'S.10.-40./L1,70,200,LE,CEN              '
%          mtext: 'M.                                      '
%          rtext: 'R.                                      '
%       stim_par: [17x1  double]      stim params (see below)
%         ranges: [10x1  double]      tone range params (see below)
%            rom: [10x1  double]      ROM, RMF PARAMS: start block;
%										no. blocks; 10(ROM) or
%										1000(RMF)*syncfreq (Hz);
%										rate divisor; 0 or RMF channel.
%										rom(1:5) for left ch., rom(6:10) for right
%    spect_level: [ 2x1  double]      spect. level of noise re tone
%										for left (1) and right (2) ears
%      mode_line: [ 8x1  double]      number of reps of each stimulus;
%										ear: 1-RI, 2-LE, 3-LE+RI, 4-ROT, 5-RIM, 6-LEM, 7-BO;
%										presentation order 0-SEQ, 1-RND, 2-CEN;
%										stimulus repetition period, ms;
%										sync source: 0-no sync, 1-SYNR, 2-SYNL;
%										0  (1 if EP recorded during spike data);
%										type of stimulus range:   
%										   0-no range  1-frequency  2-atten.
%										   3-delay     4-ear (ROT)  5-Tucker/Davis file rotation
%										   (TDT noise band freq sweeps coded as 1);
%										section count:  1 (RND, CEN)  or
%											no. of reps/stimulus (SEQ)
%          resol: 10                  data resolution, microsec./clock tick
%       esw_word: 36                  electronic switch hardware word
%        PC_line: '200,50/                            ' PC communication line
%											used in ST programs
%         tc_par: [10x1  double]      TC parameters: number of data points in the TC;
%										number of trials meeting criterion for threshold;
%										ratio of search frequency step size to TC frequency step size;
%										index of BF in the TC data;
%										threshold (real*4) dB;
%										BF (real*4);
%										spontaneous rate (real*4);
%										total number of stimuli (real*4);
%										10*spont rate;
%										criterion, driven - spont spikes for a threshold;
%		Error: 'No error' if OK, 'File error' if error opening or reading file.
%	                             'No data' if file contains no data


% Stimulus parameters:

% stim_par(n) = 

% 1 - left stimulus: 0-tone, 3-noise, 4-click, 6-RMF, 7-Tucker-Davis (incl.
%       filrot), 8-ROM farm, 9-noisebands, 10-2bands,
%	  11-no stimulus
% 2 - left stimulus: tone freq. (kHz), -1 if freq. range (see ranges(n) below),
%      ROM or RMF sync freq (Hz), =nn for 'XRnn', otherwise 0
% 3,4 - same for right stimulus
% 5,6 - left attenuation: (low atten, high atten) or (atten, -1).
%	     atten set to 140 dB if no stimulus.
% 7,8 - same for right attenuation
% 9-12 - left timing: (duration ms, PF low delay (ms),
%         PF high delay or 0 if no range, gate (-1 NG, 0 UG, 1 gated)
%13-16 - same for right timing
% 17 - 0 or attenuation step size in TC


% ranges(n) =  (specify frequency ranges in left and right channels)

% 1-5 - left ch. frequency range: low freq (kHz), high freq, center
%         freq, step size (Hz if linear, multiplier if log), 1 if log
%		  or 0 if linear



	endian='ieee-le'; %data from little endian => integers flipped

% Open data file
	fname = sprintf('%s.DAT',filename);
	fptr = fopen(fname,'r',endian);
	if fptr < 0
		fprintf(1, '*** Error opening file %s, probably directory error.***', fname)
		parameterblock = struct('Error', 'No data');
	else
% Terminal of the if above is at the end of the program


% Read
	junk1 = fread(fptr, 1, 'int16');
% short = 1 if stim text is S.xxxx, in ltext only
%       = 0 if 3 lines of stim text in ltext, mtext, and rtext
	short = fread(fptr, 1, 'int16');
	junk1 = fread(fptr, 1, 'int16');
% prog_name is 2-char name of program taking data
	prog_name = char(fread(fptr,2,'char')');
	junk1 = fread(fptr, 1, 'int16');
% date as 'mm/dd/yy'
	date = char(fread(fptr,8,'char')');
	junk1 = fread(fptr, 1, 'int16');
% track, unit, and picture number as  '  3.04, P103  ' (not adusted any way)
	unit = labstr(char(fread(fptr,14,'char')'));
	junk1 = fread(fptr, 1, 'int16');
% L. or S. stimulus text line ('XR4W(145.),90,200,LE,SEQ,1000,SYNL')
%   terminated with 0 character.
	ltext = labstr(char(fread(fptr,40,'char')'));
% Stimulus parameters (see parm.m for definitions)
	stim_par = fread(fptr, 17, 'int16');
	filnam = char(fread(fptr, 10, 'char')');
% File description:
%	1 - Number of blocks (512 byte chunks) in the file, counting param bl.
%	2 - Picture number
%	3 - 1
%	4 - data format: 1-dotdisplay, 4-tuning curve, 19-analog averager,
%	                30-acoustic calib
	file_descr = fread(fptr,4,'int16');
	junk4 = fread(fptr, 4, 'int16');
% Left and right tone range parameters and ROM/RMF params  (see parm.m)
	ranges = fread(fptr, 10, 'int16');
	rom = fread(fptr, 10, 'int16');
% Spectrum levels of noise stimulus, re 0-dB tone
	spect_level = fread(fptr, 2, 'int16');
% Mode-line parameters
	mode_line = fread(fptr, 8, 'int16');
	junk10 = fread(fptr, 10, 'int16');
% Spike timing resolution (microseconds/clock tick)
	resol = fread(fptr, 1, 'int16');
% TC program parameters (see parm.m). Note floating point params
% should be returned correctly.
	tc_par1 = fread(fptr,4,'int16');
	tc_par2 = fread(fptr,4,'real*4');
	tc_par3 = fread(fptr,2,'int16');
	tc_par = [tc_par1' tc_par2' tc_par3']';
	junk5 = fread(fptr, 5, 'int16');
% Mode and Right stimulus text
	mtext = labstr(char(fread(fptr,40,'char')'));
	rtext = labstr(char(fread(fptr,40,'char')'));
	junk1 = fread(fptr, 1, 'int16');
% PC text line, for controlling commains and filrot
	PC_line = labstr(char(fread(fptr, 38, 'char')'));
	junk57 = fread(fptr, 56, 'int16');
% Electronic switch control words
	esw_word = fread(fptr, 1, 'int16');
	junk6 = fread(fptr, 6, 'int16');
% Time of day at start of picture ('hh:mm:ss')
	time = char(fread(fptr, 8, 'char')');
	junk1 = fread(fptr, 1, 'int16');

% Convert pseudofloating frequencies. Frequencies are left as
% floating kHz. Delays are left as floating ms.
	if stim_par(1)==0 & stim_par(2)~=-1
		stim_par(2) = pffc1(stim_par(2));
	end
	if stim_par(3)==0 & stim_par(4)~=-1
		stim_par(4) = pffc1(stim_par(4));
	end
	for j=10:14:4
	   stim_par(j) = 10.*pffc1(stim_par(j));
	   if stim_par(j+1)~=-1
		  stim_par(j+1) = 10.*pffc1(stim_par(j+1));
	   end
	end
	for j=1:3
		ranges(j) = pffc1(ranges(j));
		ranges(j+5) = pffc1(ranges(j+5));
	end

% Compute the frequency-step multiplier in case of log frequency
% range
	if ranges(5)==1
		ranges(4) = 1 + ranges(4)/32768.;
	end
	if ranges(10)==1
		ranges(9) = 1 + ranges(9)/32768.;
	end

% Construct the output structure

	parameterblock=struct('date', date, 'time', time, 'unit', unit, ...
	   'prog_name', prog_name, 'file_descr', file_descr, 'short', short, ...
	   'ltext', ltext, 'mtext', mtext, 'rtext', rtext, 'stim_par', stim_par, ...
	   'ranges', ranges, 'rom', rom, 'spect_level', spect_level, ...
	   'mode_line', mode_line, 'resol', resol, 'esw_word', esw_word', ...
	   'PC_line', PC_line, 'tc_par', tc_par, 'Error', 'No error');

% OK, that's it.
	fclose(fptr);
	
% Following end terminates the if statement after the fopen() command.
	end

return

	
	
function flnum = pffc1(ipfn)

% Takes pseudofloating number ipfn, from lab-system parameter block,
% and returns a floating number flnum.
% Pseudofloating numbers are 16-bit strings of the form
%      ee mmmmmmmmmmmmmm
%      12 34567890123456
% and are interpreted as
%    number = mmmmmmmmmmmmmm / 10^(2+ee)
% Thus if mmmmmmmmmmmmmm = 8769 and ee = 1, number is 8.769

% dec2bin() doesn't handle negative numbers, so deal with the sign
% bit separately
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

	
	
function ostr = labstr(str)

% Searches for lab parameter block end of line (mod(ichar('x'),64)=0)
% and blank pads characters from there on.

	nc = size(str,2);
	for j=1:nc
		if mod(fix(str(j)),64) == 0
			for j1=j:nc
				ostr(j1) = ' ';
			end
			break
		else
			ostr(j) = str(j);
		end
	end
return
	

%[tcdata,bf,thresh]=read_tc_cal(filename,calfile)
%	Usage: type [tcdata,bf,thresh]=read_tc_cal('a3','a2') to read tuning curve 
%	data from tuning curve file a3.DAT, with a2.DAT as the 
%	calibration file.  Program outputs tuning curve values, bf, threshold

function [tcdata_cal,bf,thresh]=read_tc_cal(filename,calfile)

	[pb,FILEPNTRXYZ]=labfile2(filename);
   [Adec,count]=fread(FILEPNTRXYZ,256,'uint16');
   Abin=dec2bin(Adec,16);
   
   %convert data to correct format
   j=0;							
   for i=1:2:count
      j=j+1;
      exponent=Abin(i,3:9);
      if exponent=='0000000'
         rawval=sprintf('.0%s%s',Abin(i,10:16),Abin(i+1,1:16));
      else
         rawval=sprintf('.1%s%s',Abin(i,10:16),Abin(i+1,1:16));
      end
      val(j)=bin2decfloat(rawval)*2^(bin2dec(exponent));
   end
   num1=j;
   
	%convert single line of data to 2 columns (BF,atten)
   j=0;				
   for i=1:2:num1
      if val(i)~=0
         j=j+1;
         tcdata(j,1:2)=[val(i) -1*val(i+1)];   
      end
   end
   num2=j;
   
   %adjust attenuations to absolute values with calibration file
   calfile=sprintf('%s.DAT',calfile);
   fid=fopen(calfile,'rt');
   if fid <= 0
		fprintf(1,'\n***Error opening %s\n', calfile)
   end
   
   go=0;
   while go==0
      junk=fgets(fid);
      go=strncmp(junk,'   Freq      SPL    phase(/PI)',30);
   end
   calib = (fscanf(fid,' %f',[3, inf]))';
   for i=1:num2
      cal_level=0;
      for j=1:size(calib,1)-1
         if tcdata(i,1)>calib(j,1) & tcdata(i,1)<=calib(j+1,1)
            cal_level=(tcdata(i,1)-calib(j,1))/(calib(j+1,1)-calib(j,1)) ...
               *(calib(j+1,2)-calib(j,2)) + calib(j,2);	%linear interpolation of frequency
            															% to get the appropriate calib.level
            break;
         end
      end
      tcdata_cal(i,1:2)=[tcdata(i,1) cal_level+tcdata(i,2)];	%calculate absolute threshold
   end
   
   %find bf, threshold
   [thresh,index]=min(tcdata_cal(:,2));
   bf=tcdata_cal(index,1);
   
   %plot tuning curve
   %figure;
   subplot(313);
   semilogx(tcdata_cal(:,1),tcdata_cal(:,2),'o:');  
   ylabel('dB SPL'); xlabel('frequency (kHz)');
   axis([0.03 39 0 110]);
   %figure;
   %subplot(313);
   %semilogx(calib(:,1),calib(:,2),'o:');
   %axis([0.03 39 0 110]);
   %ylabel('dB SPL'); xlabel('frequency (kHz)');
   %fclose('all');

return
   
   
   
   
function [pb,FILEPNTRXYZ] = labfile2(filename)

% LABFILE: When called as pb = labfile('a47'), opens file a47.DAT, returns parameter
% block in pb. File is left open so fread(FILEPNTRXYZ, nn, 'int16') can be used
% to return data. On file errors FILEPNTRXYZ is set to -1 and the pb consists of
% only pb.Error = 'File error'. If the parameter block is read correctly, then
% pb.Error = 'No error' (see read_parm()).

%global FILEPNTRXYZ

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


function [decval] = bin2decfloat(rawval)
%function for converting a non-integer binary number to decimal format
%useful only for numbers of form 0.xxxxxxxxxxxx

	num_length=size(rawval,2);
   decval=0;
   for i=2:num_length
      digit_exp = -1*(i-1);
      digit_val = str2num(rawval(i))*(2^digit_exp);
      decval = decval + digit_val;
   end
   
return



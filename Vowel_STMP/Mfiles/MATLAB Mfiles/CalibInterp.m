function max_dBSPL=CalibInterp(cal_freq_kHz,calib)
% FILE: CalibInterp.m
% Created 6/25/02 M. Heinz
% Usage: max_dBSPL=CalibInterp(cal_freq_kHz,calib)
%
%%%% Find Max dB SPL at cal_freq
%
% Does simple interpolation to find the calibration level
% at any freq from a calibration file 
% calib: 1st col=freq (kHz), 2nd col=dB SPL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% If cal_freq if out of calibratino range, use end points of calibration data and display warning
if cal_freq_kHz<min(calib(:,1))  % 
    cal_freq_kHz=min(calib(:,1));
    warndlg(sprintf('Calib_Interp: cal_freq (%.1f kHz) is out of range - using endpoints of calibDATA',cal_freq_kHz),'CALIB Extrapolation','modal');
elseif cal_freq_kHz>max(calib(:,1))  % 
    cal_freq_kHz=max(calib(:,1));
    warndlg(sprintf('Calib_Interp: cal_freq (%.1f kHz) is out of range - using endpoints of calibDATA',cal_freq_kHz),'CALIB Extrapolation','modal');
end

% Find closest frequency for which calibration data exists
[x,i]=min(abs(cal_freq_kHz-calib(:,1)));  % i=index of closest point
    
if(calib(i,1)<cal_freq_kHz)
   %linear interpolation of frequency
   max_dBSPL=(cal_freq_kHz-calib(i,1))/(calib(i+1,1)-calib(i,1))* ...
      (calib(i+1,2)-calib(i,2)) + calib(i,2);  
elseif(calib(i,1)>cal_freq_kHz)
   max_dBSPL=(cal_freq_kHz-calib(i-1,1))/(calib(i,1)-calib(i-1,1))* ...
      (calib(i,2)-calib(i-1,2)) + calib(i-1,2);	
elseif(calib(i,1)==cal_freq_kHz)
   max_dBSPL=calib(i,2);
end

return;

    

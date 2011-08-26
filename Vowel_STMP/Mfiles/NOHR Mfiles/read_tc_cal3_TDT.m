function [timeX,tcdata_cal,bf,thresh,handTCbf,handTCthresh,calib]=read_tc_cal3_TDT(TCfile,CALfile,plotYES)
% [timeX,tcdata_cal,minTCbf,minTCthresh,handTCbf,handTCthresh,calib]=read_tc_cal3_TDT(TCfile,CALfile,plotYES)
% Program outputs tuning curve values, bf, threshold, calibration data, and time of TC
%
% Modified by MGH to return tcDATE and tcTIME
% Modified by MGH 102501 for TDT Data
% Modified 12/6/01 to call TCcommon.m to do things common to both PDP and TDT


eval(['tc=' TCfile ';'])
eval(['cal=' CALfile ';'])

%%%tcdata == [freq Atten]
%% Only take those frequencies for which a response was found
tcdata=tc.TcData(find((tc.TcData(:,1)~=0)&(tc.TcData(:,2)>tc.Stimuli.file_attlo)),1:2);
tcdata(:,2)=-tcdata(:,2);   % convert to negative attenuations, consistent with PDP
num2=size(tcdata,1);   %%%num2=length tcdata

%adjust attenuations to absolute values with calibration file
calib = cal.CalibData(:,1:3);

% Return time
if(isfield(tc.General,'time'))
    timeX=tc.General.time;
else
    timeX='NOT SAVED';
end

TCcommon

if isfield(tc,'Thresh')
   handTCbf=tc.Thresh.BF;
   handTCthresh=tcdata_cal(find(tcdata_cal(:,1)==handTCbf),2);
   % This hack is needed for AN model data, where TC is fake and so the BF is not a real point on curve
   if isempty(handTCthresh)
      handTCthresh=CalibInterp(handTCbf,calib)-tc.Thresh.thresh;
      disp('Threshold (dB SPL) is artificial for this unit, bc BF does not lie on TCdata (e.g., AN model data with mad up TC)')
   end 
else
	handTCbf=[];
	handTCthresh=[];
end

return
   
   

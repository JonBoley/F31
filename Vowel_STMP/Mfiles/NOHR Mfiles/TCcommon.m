% File: TCcommon.m
%
% Common TC processing for both PDP and TDT systems
% Called from: read_tc_cal3_PDP.m and read_tc_cal3_TDT.m
%
% Created: 12/6/01 to ce called to do things common to both PDP and TDT


%%%%%%%% Start general: assumes tcdata, calib, and num2 stored

for i=1:num2
    cal_level=0; tag=0;
    for j=1:size(calib,1)-1
        if tcdata(i,1)>calib(j,1) & tcdata(i,1)<=calib(j+1,1)
            cal_level=(tcdata(i,1)-calib(j,1))/(calib(j+1,1)-calib(j,1)) ...
                *(calib(j+1,2)-calib(j,2)) + calib(j,2);	%linear interpolation of frequency
            % to get the appropriate calib.level
            tag=1;
        end
    end
    if tag==1
        tcdata_cal(i,1:2)=[tcdata(i,1) cal_level+tcdata(i,2)];	%calculate absolute threshold
    elseif tag==0
        tcdata_cal(i,1:2)=[tcdata(i,1) calib(j,2)+tcdata(i,2)];
    end
end

%find bf, threshold
[thresh,index]=min(tcdata_cal(:,2));
bf=tcdata_cal(index,1);

%plot tuning curve if asked for
if(plotYES)
    set(0,'DefaultTextInterpreter','none');

    h10=figure(10); clf
    set(h10,'Position',[228 4 560 694])
    h2=subplot(212);
    semilogx(tcdata_cal(:,1),tcdata_cal(:,2),'x:');  
    ylabel('dB SPL'); xlabel('frequency (kHz)'); title(strcat('Tuning Curve:   ',TCfile))
    axis([0.03 39 0 120]);
    grid on
    %%%%%% ADDED 2/19/01 from read_tc_cal (MGH)
    h1=subplot(211);
    semilogx(calib(:,1),calib(:,2),'x:');
    axis([0.03 39 60 140]);
    grid on
    ylabel('dB SPL'); title(strcat('Calibration:   ',CALfile)); %xlabel('frequency (kHz)');
    disp(sprintf('Mean Calibration (%s) Level = %.2f dB SPL',CALfile,mean(calib(:,2))))
    text(.1,65,sprintf('Mean Calibration Level = %.2f dB SPL',mean(calib(:,2)))) 
end


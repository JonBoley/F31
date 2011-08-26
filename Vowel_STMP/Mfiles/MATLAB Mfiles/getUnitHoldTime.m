function holdtime_min=getUnitHoldTime(track,unit);
% function holdtime_min=getUnitHoldTime(track,unit);
% Created: M. Heinz 08Nov2004
% For: NOHR
%
% Returns Hold Time (difference between 1st and last data file for a given unit)
%
% 

datafiles=dir(sprintf('*u%d_%02d*',track,unit));

time1=datafiles(1).date(end-7:end);
time2=datafiles(end).date(end-7:end);

sec1=str2num(time1(1:2))*60*60+str2num(time1(4:5))*60+str2num(time1(7:8));
sec2=str2num(time2(1:2))*60*60+str2num(time2(4:5))*60+str2num(time2(7:8));

if sec2<sec1
   sec2=sec2+24*60*60;
end

holdtime_min=(sec2-sec1)/60;

holdtime_minCRIT=5*60;
if holdtime_min>(holdtime_minCRIT)
   beep
   disp(sprintf('   ***** Unit: %d.%02d has a holdtime > %d min (=%.2f min) - must be artifact!!',track,unit,holdtime_minCRIT,holdtime_min))
   holdtime_min=input('Enter holdtime by hand (in min), or Return to set to NaN: ');
   disp(sprintf('holdtime_min set to %.2f min BY HAND',holdtime_min))
   if isempty(holdtime_min)
      holdtime_min=NaN;
   end
end


return;

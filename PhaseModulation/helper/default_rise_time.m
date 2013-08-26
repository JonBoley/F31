function rise_fall_time = default_rise_time(stm_duration)
%

% AF 11/9/01

if (stm_duration < 100 & stm_duration > 10)
   rise_fall_time = 5;
elseif (stm_duration <= 10)
   rise_fall_time = 1;
else
   rise_fall_time = 10;
end   

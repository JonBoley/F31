function freq = CalcFixedFreq
% freq = Calcfreq
%
% M Heinz 07July2004
% For use with NOHR experiments.  Calculates a fixed frequency based on the current_unit_bf

% global NelData

if (current_unit_bf>=.25)&(current_unit_bf<1.0)
   freq=0.5;
elseif (current_unit_bf>=1.0)&(current_unit_bf<4.0)
   freq=2.0;
elseif (current_unit_bf>=4.0)&(current_unit_bf<16.0)
   freq=8.0;
else
   freq=[];
end

return
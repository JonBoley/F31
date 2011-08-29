function default_inloop_status(h,msg,plot_info)
%

% AF 10/10/01
global NumDrivenSpikes nstim % MH 29Jun2004

if isempty(nstim)
   nstim=1;
end

str = '';
if (isnumeric(msg))
   % msg is the line number
   frmt = ['%s: ' plot_info.var_frmt '%s'];
   if (isfield(plot_info,'var_labels'))
      str = sprintf(frmt, plot_info.var_name, plot_info.var_labels{msg}, plot_info.var_unit);
   else
      str = sprintf(frmt, plot_info.var_name, plot_info.var_vals(msg), plot_info.var_unit);
   end
elseif (isstr(msg))
   str = msg;
end

if nstim == 0 | nstim ==1
   str=sprintf('%s\nNumber of Driven Spikes: %d',str,NumDrivenSpikes(1)); % MH: 29Jun2004
else
   if length(NumDrivenSpikes)==1
      NumDrivenSpikes=NumDrivenSpikes*ones(1,3);
   end
   str=sprintf('%s\nDriven Spikes/Condition: min=%d, mean=%.f, max=%d',str,NumDrivenSpikes(1),NumDrivenSpikes(2),NumDrivenSpikes(3)); % MH: 29Jun2004
end
set(h,'String',str);

% set(h,'FontSize',10,'pos',[124.6 1.07692 68.2 3.15385])  %MH: 11Nov2004 to make room  (THIS IS ORIGINAL SETTINGS)
set(h,'FontSize',9,'pos',[124.6 1.1 68.2 2.4])  %MH: 11Nov2004 to make room

drawnow;
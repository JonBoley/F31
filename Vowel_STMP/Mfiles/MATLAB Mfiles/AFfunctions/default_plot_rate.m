function rate_params = default_plot_rate(index,rate_params,plot_info)
%

% AF 10/10/01

global spikes

% To speed the plot process we don't check rate_params
if (index > 0)
   indices = max(1,index-1):index;  % recalculate the last index in case we missed spikes
   % indices = 1:index;  % recalculate the last index in case we missed spikes
   for index = indices
      ind = find(spikes.times{1}(1:spikes.last(1),1) == index);
      driven = length(find(spikes.times{1}(ind,2) <= rate_params.stim_dur));
      spont  = length(ind) - driven;
      xd = get(rate_params.cache.hdriven,'Xdata');
      yd = get(rate_params.cache.hdriven,'Ydata');
      xs = get(rate_params.cache.hspont,'Xdata');
      ys = get(rate_params.cache.hspont,'Ydata');

      yd(rate_params.plot_order(index)) = driven / rate_params.stim_dur;
      ys(rate_params.plot_order(index)) = spont / (rate_params.line_dur-rate_params.stim_dur);
      % fprintf('index=%d, driven=%d\n', index, driven);
   end
   curr_maxy = max(get(rate_params.cache.haxes,'YLim'));
   maxy = max([max(yd) max(ys)]);
   if (maxy > 0.9*curr_maxy)
      set(rate_params.cache.haxes,'YLim',[0 maxy/0.9]);
   end
   % fprintf('%f %f\n', [xd ;yd]);
   set(rate_params.cache.hdriven,'Xdata',xd,'Ydata',yd);
   set(rate_params.cache.hspont ,'Xdata',xs,'Ydata',ys);
   drawnow;
   return;
end

% No spikes to plot. Create axes if necessary and set properties.
if ((exist('plot_info','var') ~= 1))
   plot_info = [];
end
if ((exist('rate_params','var') ~= 1))
   rate_params = [];
end

% Create 
if ((~isfield(rate_params,'cache')) |  (~isfield(rate_params.cache,'haxes')))
   haxes   = axes('position',[ .1 .1 .3 .25]);
   hdriven = plot(-5,-5,'-');
   hold on;
   hspont  = plot(-5,-5,'-');
   hxlabel = xlabel('');
   hylabel = ylabel('Discharge rate (sp/s)');
   rate_params.cache.haxes   = haxes;
   rate_params.cache.hdriven = hdriven;
   rate_params.cache.hspont  = hspont;
   
   %%%% Set stimulus related properties (possible only when creating new axes)
   if (isempty(plot_info))
      plot_info = default_inloop_plot_info;
   end   
   label = plot_info.var_name;
   if (~isempty(plot_info.var_unit))
      label = [label ' (' plot_info.var_unit ')'];
   end
   set(hxlabel,'String',label);
   if (~isempty(plot_info.XYprops))
      set(haxes, strcat('X',fieldnames(plot_info.XYprops))', struct2cell(plot_info.XYprops)');
   end
   rate_params.cache.var_vals = plot_info.var_vals;
   [xdata order] = sort(plot_info.var_vals);
   [dummy rate_params.plot_order] = sort(order);
   set(rate_params.cache.hdriven,'Xdata',xdata,...
      'Ydata',repmat(NaN,size(rate_params.cache.var_vals)));
   set(rate_params.cache.hspont,'Xdata',xdata,...
      'Ydata',repmat(NaN,size(rate_params.cache.var_vals)));
end

%%%% Set general properties 
if (~isfield(rate_params,'props'))
   dflt_params = default_plot_rate_params;
   rate_params.props = dflt_params.props;
end   
props = rate_params.props;
for f = {'axes','driven','spont','xlabel','ylabel'}
   if (isfield(props,f{1}))
      eval(['set(h' f{1} ', props.' f{1} ');']);
   end
end

function raster_params = default_plot_raster(spikes,raster_params,plot_info)
%

% AF 10/3/01

% To speed the plot process we don't check raster_params
if ((exist('spikes','var') == 1) & (~isempty(spikes)))
   if (~isempty(spikes{1}))
      set(raster_params.cache.hline,'xdata',spikes{1}(:,2),'ydata',raster_params.cache.var_vals(spikes{1}(:,1)));
   end
   return;
end

% No spikes to plot. Create axes if necessary and set properties.
if ((exist('plot_info','var') ~= 1))
   plot_info = [];
end
if ((exist('raster_params','var') ~= 1))
   raster_params = [];
end

% Create 
if ((~isfield(raster_params,'cache')) |  (~isfield(raster_params.cache,'haxes')))
   haxes = axes('position',[.2 .40 .6 .55]);
   hline = plot(-5,-5,'.');
   hxlabel = xlabel('Time (sec)');
   hylabel = ylabel('');	
   raster_params.cache.haxes = haxes;
   raster_params.cache.hline = hline;
end

%%%% Set general properties 
if (~isfield(raster_params,'props'))
   dflt_params = default_plot_raster_params;
   raster_params.props = dflt_params.props;
end   
props = raster_params.props;
for f = {'axes','line','xlabel','ylabel'}
   if (isfield(props,f{1}))
      eval(['set(h' f{1} ', props.' f{1} ');']);
   end
end

%%%% Set stimulus related properties 
if (isempty(plot_info))
   plot_info = default_inloop_plot_info;
end   
label = plot_info.var_name;
if (~isempty(plot_info.var_unit))
   label = [label ' (' plot_info.var_unit ')'];
end
set(hylabel,'String',label);
if (~isempty(plot_info.XYprops))
   set(haxes, strcat('Y',fieldnames(plot_info.XYprops))', struct2cell(plot_info.XYprops)');
end
raster_params.cache.var_vals = plot_info.var_vals;

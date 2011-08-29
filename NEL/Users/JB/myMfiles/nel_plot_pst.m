function rate_params = nel_plot_pst(index,rate_params,plot_info,nCh)
%

% adapted by GE (28Jul2003) from "default_plot_rate" written by AF 10/10/01

% Added NumDrivenSpikes (MH: 29Jun2004)

global spikes NumDrivenSpikes
persistent stat_pst

if (index > 0)
   xd = [rate_params.pst.firstBin_sec:rate_params.pst.binWidth_sec:rate_params.pst.lastBin_sec];
                                                                        % time axis data of pst histogram.
   for i = 1:length(spikes.times) % number of input spike channels
      indices = max(1,index-1):index;  % recalculate the last index in case we missed spikes
                                       % line numbers that have passed since last plot update.
      for index = indices
         spike_inds = find(spikes.times{i}(1:spikes.last(1),1) == index); % indices of spikes in lines since last update.
         for j = spike_inds
            spikeTime_sec = spikes.times{i}(j, 2);
            binNumber = floor((spikeTime_sec - rate_params.pst.firstBin_sec) / rate_params.pst.binWidth_sec) + 1;
				if (binNumber<1 | binNumber>stat_pst.high_bin)  % spike time is not in an available bin.
					continue;
				else   % add a count to the appropriate bin.
					stat_pst.histo(binNumber) = stat_pst.histo(binNumber) + 1;
				end
         end
         stat_pst.collected_lines = stat_pst.collected_lines + 1;
      end
      yd = stat_pst.histo / (rate_params.pst.binWidth_sec * stat_pst.collected_lines); % normalize;
% GE 12Nov2003: was used to remove spike artifacts
%       if (rate_params.pst.nullThreshold > 0)
% 			yd(find(yd > rate_params.pst.nullThreshold)) = NaN;
% 		end
      set(rate_params.cache(i).hdriven,'Xdata',xd,'Ydata',yd);
%       set(rate_params.cache(i).hdriven,'Ydata',yd,'Xdata',xd);  % GE: swap axes
%       set(rate_params.cache(i).hdriven,'Xdata',yd);

% MH 29Jun2004
NumDrivenSpikes=round(mean(yd(find(xd<=rate_params.stim_dur)))*rate_params.stim_dur*max(index));

      curr_maxy = max(get(rate_params.cache(i).haxes,'YLim')); % check and re-scale axis if necessary.
      maxy = max(yd);
      if ((maxy > 0.9*curr_maxy) | (maxy < 0.5*curr_maxy))
         set(rate_params.cache(i).haxes,'YLim',[0 max(maxy/0.9, 100)]);
         set(rate_params.cache(i).haxes,'YTick', [0:100:1000]);
         if (maxy/0.9 > 600)
            set(rate_params.cache(i).haxes,'YTick', [0:200:1000]);
         end
      end
   end
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
if ((exist('nCh','var') ~= 1))
   nCh = 1;
end

spikes_fig(nCh);
[dummy rate_params] = spikes_fig('reset_axes',[],rate_params,plot_info);
xd = [rate_params.pst.firstBin_sec:rate_params.pst.binWidth_sec:rate_params.pst.lastBin_sec];
                                                                        % time axis data of pst histogram.
stat_pst.histo = xd;
stat_pst.histo(:) = 0;
stat_pst.collected_lines = 0;
stat_pst.high_bin = size(xd, 2);
NumDrivenSpikes=0; % MH 29Jun2004

for i = 1:nCh
   set(rate_params.cache(i).hdriven,'Xdata',xd, 'Ydata', stat_pst.histo);                                                                        
%    set(rate_params.cache(i).hdriven,'Ydata',stat_pst.histo, 'Xdata', xd);           % GE: swap axes.                                                             
	set(rate_params.cache(i).haxes,'XLim',[0 max(xd)]);
	set(rate_params.cache(i).haxes,'XTick',[0:max(xd)/5:max(xd)]);
	set(rate_params.cache(i).haxes,'XDir','normal');
   set(get(rate_params.cache(i).haxes,'XLabel'),'String','time (sec)');
   set(get(rate_params.cache(i).haxes,'YLabel'),'String','discharge rate (sp/sec)');
end

function rate_params = default_plot_rate(index,rate_params,plot_info,nCh)
%

% AF 10/10/01

global spikes NumDrivenSpikes nstim

if isempty(nstim)
   nstim=1;
end

% To speed the plot process we don't check rate_params
if (index > 0)
   indices = max(1,index-1):index;  % recalculate the last index in case we missed spikes
   for i = 1:length(spikes.times) % # of channels
      for index = indices

         %% GE/MH 06Nov2003.  This "find" line in particular causes conflict problems with the NI6052 board.
         %%    Don't know exactly why, but we suspend activity while the trigger is high (see
         %%    'data_acquisition_loop_NI'.
         ind = find(spikes.times{i}(1:spikes.last(1),1) == index);
         
         driven = length(find(spikes.times{i}(ind,2) <= rate_params.stim_dur));
         spont  = length(ind) - driven;
         xd = get(rate_params.cache(i).hdriven,'Ydata');
         yd = get(rate_params.cache(i).hdriven,'Xdata');
         xs = get(rate_params.cache(i).hspont,'Ydata');
         ys = get(rate_params.cache(i).hspont,'Xdata');
         
         yd(rate_params.plot_order(index)) = driven / rate_params.stim_dur;
         ys(rate_params.plot_order(index)) = spont / (rate_params.line_dur-rate_params.stim_dur);
         % fprintf('index=%d, driven=%d\n', index, driven);
      end
      curr_maxy = max(get(rate_params.cache(i).haxes,'XLim'));
      maxy = max([max(yd) max(ys)]);
      if (maxy > 0.9*curr_maxy)
         set(rate_params.cache(i).haxes,'XLim',[0 maxy/0.9]);
         if (maxy/0.9 > 150)
            set(rate_params.cache(i).haxes,'XTick', [0:100:1000]);
         end
      end
      % fprintf('%f %f\n', [xd ;yd]);
      set(rate_params.cache(i).hdriven,'Ydata',xd,'Xdata',yd);
      set(rate_params.cache(i).hspont ,'Ydata',xs,'Xdata',ys);

      %MH 29Jun2004
      %MH: 11Nov2004 - updated to allow interleaving of conditions (returns min, mean, max)
      if nstim == 0 | nstim ==1
         NumDrivenSpikes=sum(yd(find(~isnan(yd)))*rate_params.stim_dur)*ones(1,3);
      else
         ydINDs=find(~isnan(yd));
         Cond_nspikes=zeros(1,nstim);
         for condIND=1:nstim
            cond_ydINDs=ydINDs(find(mod(ydINDs-1,nstim)==(condIND-1)));
            Cond_nspikes(condIND)=sum(yd(cond_ydINDs)*rate_params.stim_dur);
         end
         NumDrivenSpikes=[min(Cond_nspikes) mean(Cond_nspikes) max(Cond_nspikes)];      
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
NumDrivenSpikes=0; %MH 29Jun2004

spikes_fig(nCh);
[dummy rate_params] = spikes_fig('reset_axes',[],rate_params,plot_info);

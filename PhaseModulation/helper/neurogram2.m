function neurogram2(CHdata,FIGinfo,PARAMInfo)
% FROM: ARO2005_replicate.m
% M. Heinz Jun 09, 2007
%
% General program to plot a neurogram (multi-channl data).  Setup for
% general STMP data, but should more general for all multi-channel data.
%
% CHdata.YData: is a 2-dimensial cell {ChannelIND, Pind} that
% contains the data (vectors) to be plotted on the ORDinate for each
% channel and for each CONDITION PARAMETER (e.g.,  multiple levels) to be
% plotted on top of one another.
%
% CHdata.XData: is a 1-dimensial cell {ChannelIND} that contains the
% data (vectors) to be plotted on the ABScissa for each channel and for each
% CONDITION PARAMETER (e.g.,  multiple levels) to be plotted on top of one
% another.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this data, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(CHdata.Ydata,1);
NumP=size(CHdata.Ydata,2);

%%%% Plot data
PARAMcolors={'b','b','g','g','c','c','k','k','y','y'};
PARAMmarkers={'none','o','none','x','none','s','none','^','none','*'};
PARAMlines={'-','none','-','none','-','none','-','none','-','none'};
FIG.LEGfontSize=12;
FIG.ANNOTfontSize=14;
FIG.LABELfontSize=16;
FIG.TITLEfontSize=10;
FIG.FEATfontSize=14;
FIG.VERTLABfontSize=10;
FIG.HORLABfontSize=14;

% Default is solid line only
if ~isfield(FIGinfo,'LineStyle')   
   FIGinfo.LineStyle='-';
end
if ~isfield(FIGinfo,'Marker')   
   FIGinfo.Marker='none';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET plotting params based on data
%
%% Find appropriate YLIMITS based on channel_values and chGAIN
chGAIN=FIGinfo.chGAIN; % # of channels covered by max PERhist - USED TO SCALE CHdataY
highCH=max(CHdata.CHvals);
lowCH=min(CHdata.CHvals);
if FIGinfo.Ylog==1
   ch_logCHwidth=log10(highCH/lowCH)/(NumCH-1);  % log10 channel width
   ch_YMAX=10^((NumCH-1+chGAIN)*ch_logCHwidth)*lowCH;    % sets an extra (GAIN-1) log channel widths (to allow for chGAIN>1)
else
   ch_CHwidth=(highCH-lowCH)/(NumCH-1);  % linear channel width
   ch_YMAX=(NumCH-1+chGAIN)*ch_CHwidth+lowCH;    % sets an extra (GAIN-1) channel widths (to allow for chGAIN>1)
end
YLIMITS=[lowCH ch_YMAX];  % Used for all plots
YTICKS=sort(round(CHdata.CHvals*100)/100);

%% Find maximum value acros all Ydata o be plotted - FOR SCALING
Y_MAXval=-Inf;
for ChanIND=1:NumCH
   for Pind=1:NumP
      Y_MAXval=max([Y_MAXval max(trifilt(CHdata.Ydata{ChanIND,Pind},CHdata.TriFiltWidth))]);
   end  % end Param
end  % end Chan


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% NEED TO FIND Y_MAXval based on smoothed data!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LEGtext='';

BOLDLINEwidth=4;
REGLINEwidth=2;

for CHind=1:length(CHdata.CHvals)
   if ismember(CHdata.BOLDchs,CHind), LINEwidth=BOLDLINEwidth;
   else               LINEwidth=REGLINEwidth;       end
   for Pind=1:NumP
      if CHind==1
         PARAMval=PARAMInfo.param_values(Pind);
         if mod(Pind,2)
             LEGtext{length(LEGtext)+1}=sprintf('%.f %s',PARAMval,FIGinfo.param_units);
         end
      end
      
      % This normalization plots each signal the same size on a log scale
      if FIGinfo.Ylog==1
         NormFact=(10^(chGAIN*ch_logCHwidth)-1)*CHdata.CHvals(CHind)/Y_MAXval;
         semilogy(CHdata.Xdata{CHind,Pind}, ...
            trifilt(CHdata.Ydata{CHind,Pind},CHdata.TriFiltWidth)*NormFact+CHdata.CHvals(CHind), ...
            'LineWidth',LINEwidth,'LineStyle',PARAMlines{Pind},'Marker',PARAMmarkers{Pind},'Color',PARAMcolors{Pind})
         hold on
      else
         error('SETUP LINEAR PLOTTING!!!!')
      end
   end
   if CHind==1
      hleg=legend(LEGtext,1);
   end
end
set(hleg,'FontSize',FIG.LEGfontSize)
xlabel(FIGinfo.Xlabel_text,'FontSize',FIG.LABELfontSize)
ylabel(FIGinfo.Ylabel_text,'FontSize',FIG.LABELfontSize)
title(FIGinfo.title,'FontSize',FIG.TITLEfontSize)
xlim(FIGinfo.XLIMITS) 
set(gca,'YTick',YTICKS,'YTickLabel',YTICKS)
ylim(YLIMITS)  % Same Ylimits for all plots

text(0,1,FIGinfo.CONDlabel_text,'Units','norm','Vert','bottom','FontSize',FIG.FEATfontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ADD all Vertical Labels and/or lines
if isfield(FIGinfo,'VerticalLabels')
   for ii=1:length(FIGinfo.VerticalLabels.Text)
      if (FIGinfo.VerticalLabels.Xvals{ii}>=FIGinfo.XLIMITS(1))&(FIGinfo.VerticalLabels.Xvals{ii}<=FIGinfo.XLIMITS(2))
         text(FIGinfo.VerticalLabels.Xvals{ii},YLIMITS(1),FIGinfo.VerticalLabels.Text{ii}, ...
            'Units','data','Vert','top','Horiz','center','FontSize',FIG.VERTLABfontSize,'Color',FIGinfo.VerticalLabels.color{ii})
         if ~isempty(FIGinfo.VerticalLabels.linestyle{ii})
            plot(FIGinfo.VerticalLabels.Xvals{ii}*ones(1,2),YLIMITS,'Color',FIGinfo.VerticalLabels.color{ii}, ...
               'LineStyle',FIGinfo.VerticalLabels.linestyle{ii},'LineWidth',FIGinfo.VerticalLabels.linewidth{ii})
         end
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ADD all Horizontal Labels and/or lines
if isfield(FIGinfo,'HorizontalLabels')
   for ii=1:length(FIGinfo.HorizontalLabels.Text)
      if (FIGinfo.HorizontalLabels.Yvals{ii}>YLIMITS(1))&(FIGinfo.HorizontalLabels.Yvals{ii}<YLIMITS(2))
         text(1.005*FIGinfo.XLIMITS(2),FIGinfo.HorizontalLabels.Yvals{ii},FIGinfo.HorizontalLabels.Text{ii}, ...
            'Units','data','Vert','Middle','Horiz','left','FontSize',FIG.HORLABfontSize,'Color',FIGinfo.HorizontalLabels.color{ii})
         if ~isempty(FIGinfo.HorizontalLabels.linestyle{ii})
            plot(FIGinfo.XLIMITS,FIGinfo.HorizontalLabels.Yvals{ii}*ones(1,2),'Color',FIGinfo.HorizontalLabels.color{ii}, ...
               'LineStyle',FIGinfo.HorizontalLabels.linestyle{ii},'LineWidth',FIGinfo.HorizontalLabels.linewidth{ii})
         end
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ADD Other Lines
if isfield(FIGinfo,'OtherLines')
   for ii=1:length(FIGinfo.OtherLines.Xvals)
      plot(FIGinfo.OtherLines.Xvals{ii},FIGinfo.OtherLines.Yvals{ii},'Color',FIGinfo.OtherLines.color{ii}, ...
         'LineStyle',FIGinfo.OtherLines.linestyle{ii},'LineWidth',FIGinfo.OtherLines.linewidth{ii})
   end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% ADD Other Markers
% if isfield(FIGinfo,'OtherMarkers')
%    for ii=1:length(FIGinfo.OtherLines.Xvals)
% % 		plot(SACSCCmetrics.CDscc_usec/1000,SACSCCfunctions.SCC_AB_avg(find(SACSCCfunctions.delays_usec==SACSCCmetrics.CDscc_usec)),'bx','MarkerSize',10,'LineWidth',2)
% 		plot(FIGinfo.OtherMarkers.Xvals{ii},FIGinfo.OtherMarkers.Yvals{ii},'Color',FIGinfo.OtherMarkers.color{ii}, ...
% 			'Marker',FIGinfo.OtherMarkers.marker{ii},'MarkerSize',FIGinfo.OtherMarkers.markerSize{ii}, ...
% 			'LineWidth',FIGinfo.OtherLines.linewidth{ii})
%    end
% end
   
set(gca,'FontSize',FIG.ANNOTfontSize)


hold off

return;


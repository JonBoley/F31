function [Thresh_dBSPL_ret,BF_kHz_ret,Q10_ret] = plotTCs(PIClist,CALIBpic,PLOTyes)
% FILE: plotTCs
% Modified from : verifyBFQ10.m
% Usgae: [Thresh_dBSPL_ret,BF_kHz_ret,Q10_ret] =
% plotTCs(PIClist,CALIBpic,PLOTyes)
% Just a simple way to plot TCs from a given list of TC pics
%
% Modified on: 10May2007 M. Heinz for SAC_XAC Analysis
%
% Modified 11Feb2005 M. Heinz for NOHR Data
%
% Created 7/31/02: for choosing BF and verofying Q10
%
% 1) Picks BF & Threshold from actual data points
% 2) Generates a smoothed TC (without bias at BF) and saves as unit.tcdata(:,3)
% 3) Finds Q10 based on smoothed TC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('PIClist')
   % EXP:
   PIClist = [63];
   CALIBpic = 2;

   % EXP: 110306
   PIClist = [29 39 40 41 42 43 44 45 49];
   CALIBpic = 2;
   
   
end

if ~exist('PLOTyes')
   PLOTyes=1;
end

%%% READ in Calib Data
xCAL=loadPic(CALIBpic);
CalibData = xCAL.CalibData(:,1:2);

numTCs=length(PIClist);
TrackUnit=getTrackUnit(getFileName(PIClist(1)));
TRACK=TrackUnit(1);   UNIT=TrackUnit(2);

TFiltWidthTC=5;

if PLOTyes
   %set(0,'DefaultTextInterpreter','none');
   set(0,'DefaultTextUnits','data')

   colors = {'y','b','r','k','g','m','c'};

   h14=figure(14); clf;
   set(h14,'Position',[485   612   914   358])
 
   TextFontSize=5;
   DataMarkerSize=12;
   DataMarkStyle='.';
   DataFitStyle='-';

   xmin=0.03; xmax=39; ymin=-20; ymax=110;
   load normema
   legtext='';
end

for ind=1:numTCs
   PICind=PIClist(ind);
   x{ind}=loadPic(PICind);
   TCdata{ind}=x{ind}.TcData;
   TCdata{ind}=TCdata{ind}(find(TCdata{ind}(:,1)),:);  % Get rid of all 0 freqs
   TCdata{ind}=TCdata{ind}(find(TCdata{ind}(:,2)~=x{ind}.Stimuli.file_attlo),:);  % Get rid of all 'upper atten limit points'
   %% TCdata: 
   %     col 1: freq; 
   %     col 2: raw ATTENS; 
   %     col 3: raw dB SPL; 
   %     col 4: smoothed SPLS
   for i=1:size(TCdata{ind},1)
      TCdata{ind}(i,3)=CalibInterp(TCdata{ind}(i,1),CalibData)-TCdata{ind}(i,2);
   end
   TCdata{ind}(:,4)=trifilt(TCdata{ind}(:,3)',TFiltWidthTC)';
   
	%    % Set unit BF/Threshold to picked BF/Threshold
	%    [Thresh_dBSPL{ind},loc{ind}]=min(abs(TCdata{ind}(:,4)));
	%    BF_kHz{ind}=TCdata{ind}(loc{ind},1);

	% Set unit BF/Threshold to picked (during EXPERIMENT) BF/Threshold
	BF_kHz{ind}=x{ind}.Thresh.BF;
	Thresh_dBSPL{ind}=TCdata{ind}(find(TCdata{ind}(:,1)==BF_kHz{ind}),4);
	
	% Also store unit BF/Threshold based on min value
	[minThresh_dBSPL{ind},loc{ind}]=min(abs(TCdata{ind}(:,4)));
	minBF_kHz{ind}=TCdata{ind}(loc{ind},1);

	% % %% Generate smoothed TC, but avoiding upward bias at BF (tip)
   % % % Fits each side separately, and then sets equal to actual data point at BF
   % % % i.e., smoothes sides (e.g., 10 dB up) without biasing threshold at BF upward
   % % TCdata(1:loc,4)=trifilt(TCdata(1:loc,3)',TFiltWidthTC)';
   % % TCdata(loc:end,4)=trifilt(TCdata(loc:end,3)',TFiltWidthTC)';
   % % TCdata(loc,4)=TCdata(loc,3);

   % pass smoothed tcdata for q10 calculation (based on actual data point at BF, and smoothed TC otherwise
   % This avoids the bias in smoothing at the tip, i.e., raising threshold at BF
   [Q10{ind},Q10fhi{ind},Q10flo{ind},Q10lev{ind}] = findQ10(TCdata{ind}(:,1),TCdata{ind}(:,4),BF_kHz{ind});
	figure(14)

   if PLOTyes
      %%%%% PLOT TUNING CURVE
      h_line1{ind} = semilogx(TCdata{ind}(:,1),TCdata{ind}(:,3),DataMarkStyle,'MarkerSize',DataMarkerSize,'Color',colors{mod(ind,7)+1});
      hold on
      h_line2{ind} = semilogx(TCdata{ind}(:,1),TCdata{ind}(:,4),DataFitStyle,'Color',colors{mod(ind,7)+1});
      if ind == 1
         semilogx(normt(1,:),normt(2,:),'k')
         ylabel('dB SPL'); xlabel('Frequency (kHz)');
         axis([xmin xmax ymin ymax]);
         set(gca,'YTick',[0 20 40 60 80 100])
         set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
         title(sprintf('Unit: %d.%d; (Cal: P%d)',TRACK,UNIT,CALIBpic))
         if geomean(TCdata{ind}(:,1)) < 1
            Xtext=.55;
         else
            Xtext=.05;
         end
      end

      semilogx([Q10flo{ind} Q10fhi{ind}],Q10lev{ind}*ones(1,2),'-','linewidth',2,'Color',colors{mod(ind,7)+1});
      legtext{ind} = sprintf('P%d:  BF=%.3f kHz; Thr=%.1f dB SPL; Q10=%.1f',PICind,BF_kHz{ind},Thresh_dBSPL{ind},Q10{ind});

      text(Xtext,.95-.05*(ind-1),legtext{ind},'Units','norm','Color',colors{mod(ind,7)+1})
   end

   
   Thresh_dBSPL_ret(ind) = Thresh_dBSPL{ind};
   BF_kHz_ret(ind) = BF_kHz{ind};
   Q10_ret(ind) = Q10{ind};

end
hold off

   
return;

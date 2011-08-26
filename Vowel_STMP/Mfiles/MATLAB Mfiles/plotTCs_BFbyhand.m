function [Thresh_dBSPL_ret,BF_kHz_ret,Q10_ret,TCdata] = plotTCs_BFbyhand(PIClist,CALIBpic,PLOTyes,SR_sps)
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

if ~exist('PLOTyes')
   PLOTyes=1;
end
if ~exist('SR_sps')
   SR_sps=NaN;
end

%%% READ in Calib Data
xCAL=loadPic(CALIBpic);
CalibData = xCAL.CalibData(:,1:2);

numTCs=length(PIClist);
TrackUnit=getTrackUnit(getFileName(PIClist(1)));
TRACK=TrackUnit(1);   UNIT=TrackUnit(2);
[pathstr, EXPname]=fileparts(pwd);

TFiltWidthTC=5;

if PLOTyes
   %set(0,'DefaultTextInterpreter','none');
   set(0,'DefaultTextUnits','data')

   colors = {'y','b','r','k','g','m','c'};

   h14=figure(14); clf;
	set(h14,'units','norm')
	%     set(h14,'Position',[485   612   914   358])
	set(h14,'Position',[0 0 1 1])
 
   TextFontSize=8;
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
   %     col 4: smoothed SPLs
   for i=1:size(TCdata{ind},1)
      TCdata{ind}(i,3)=CalibInterp(TCdata{ind}(i,1),CalibData)-TCdata{ind}(i,2);
   end
   TCdata{ind}(:,4)=trifilt(TCdata{ind}(:,3)',TFiltWidthTC)';  % smoothed TC
   
	% Set unit BF/Threshold to picked (during EXPERIMENT) BF/Threshold
	BF_kHz{ind}=x{ind}.Thresh.BF;
	Thresh_dBSPL{ind}=TCdata{ind}(find(TCdata{ind}(:,1)==BF_kHz{ind}),4);
	
	% Also store unit BF/Threshold based on min value
	[minThresh_dBSPL{ind},loc{ind}]=min(abs(TCdata{ind}(:,4)));
	minBF_kHz{ind}=TCdata{ind}(loc{ind},1);
	
   % pass smoothed tcdata for q10 calculation
   [Q10{ind},Q10fhi{ind},Q10flo{ind},Q10lev{ind}] = findQ10(TCdata{ind}(:,1),TCdata{ind}(:,4),BF_kHz{ind});

   if PLOTyes
		%%%%% PLOT TUNING CURVE
		h_line1{ind} = semilogx(TCdata{ind}(:,1),TCdata{ind}(:,3),DataMarkStyle,'MarkerSize',DataMarkerSize,'Color',colors{mod(ind,7)+1});
		hold on
		grid on
		h_line2{ind} = semilogx(TCdata{ind}(:,1),TCdata{ind}(:,4),DataFitStyle,'Color',colors{mod(ind,7)+1},'LineWidth',2);
		h_lineQ10{ind}=semilogx([Q10flo{ind} Q10fhi{ind}],Q10lev{ind}*ones(1,2),'-','linewidth',2,'Color',colors{mod(ind,7)+1});
		if ind == 1
			semilogx(normt(1,:),normt(2,:),'k')
			ylabel('dB SPL'); xlabel('Frequency (kHz)');
			axis([xmin xmax ymin ymax]);
			set(gca,'YTick',[0 20 40 60 80 100])
			set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
			title(sprintf('Exp: %s;    Unit: %d.%d; (Cal: P%d)',EXPname,TRACK,UNIT,CALIBpic),'INterp','none')
			if geomean(TCdata{ind}(:,1)) < 1
				Xtext=.45;
			else
				Xtext=.35;
			end
		end
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Allow user to verify BF
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% mark min and book BFs
		h_minTC = semilogx(minBF_kHz{ind}*ones(1,2),minThresh_dBSPL{ind}+[5 30],'m:','LineWidth',2);
		% 		text(.1,.2,'minTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','m','units','data')
		text(minBF_kHz{ind},minThresh_dBSPL{ind}+30,'minTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center', ...
			'VerticalAlignment','bottom','Color','m','units','data')
		h_bookTC = semilogx(BF_kHz{ind}*ones(1,2),Thresh_dBSPL{ind}+[-5 -30],'c:','LineWidth',2);
		% 		text(.1,.1,'bookTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','VerticalAlignment','bottom','Color','c','units','data')
		text(BF_kHz{ind},Thresh_dBSPL{ind}-30,'bookTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center', ...
			'VerticalAlignment','top','Color','c','units','data')
		
		% setup for user to pick
		pickThr = Thresh_dBSPL{ind};
		pickBF = BF_kHz{ind};
		
		text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
		h_textBF= text(.05,.98,text_str,'Units','norm','FontSize',TextFontSize,'VerticalAlignment','top');
		h_pickBF=text(pickBF,pickThr,'\uparrow','Interpreter','tex','FontSize',16, ...
			'VerticalAlignment','top','HorizontalAlignment','center');
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%% Loop to wait for verifying BF
		x15=get(h14,'Position');  % These commands insure mouse is over Figure so that cursor keys work
		set(0,'PointerLocation',x15(1:2)+[10 x15(4)/2])
		set(gcf,'CurrentCharacter','x')
		loc = find(TCdata{ind}(:,1)==pickBF);
		if isempty(loc)
			[yy,loc]=min(abs(TCdata{ind}(:,1)-pickBF));
			disp('unit.BF NOT chosen from data points')
			pickThr = TCdata{ind}(loc,4);
			pickBF = TCdata{ind}(loc,1);
			set(h_pickBF,'Position',[pickBF pickThr]);
			text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
			set(h_textBF,'String',text_str);
		end
		
		while 1==1      % Wait for verifying BF to complete
			pause(.1)
			w = waitforbuttonpress;
			if w ~= 0   %%%%   ('Mouse-Button press')
				keypress=get(gcf,'CurrentCharacter');
				
				switch double(keypress)
					case 13  %%% 'RETURN'
						break;
					case 28  %%% 'LEFT cursor'
						loc=min(length(TCdata{ind}(:,1)),loc+1);
					case 29  %%% 'RIGHT cursor'
						loc=max(1,loc-1);
				end
			end
			
			pickThr = TCdata{ind}(loc,4);
			pickBF = TCdata{ind}(loc,1);
			% recompute Q10 from new threshold
			[Q10{ind},Q10fhi{ind},Q10flo{ind},Q10lev{ind}] = findQ10(TCdata{ind}(:,1),TCdata{ind}(:,4),pickBF);

			set(h_pickBF,'Position',[pickBF pickThr]);
			text_str = sprintf('%s %6.3f %s\n%s %4.2f %s\n%s %.1f','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL','Q10:',Q10{ind});
			set(h_textBF,'String',text_str);
% 			h_lineQ10{ind}=semilogx([Q10flo{ind} Q10fhi{ind}],Q10lev{ind}*ones(1,2),'-','linewidth',2,'Color',colors{mod(ind,7)+1});
			
			
			set(h_lineQ10{ind},'XData',[Q10flo{ind} Q10fhi{ind}],'YData',Q10lev{ind}*ones(1,2));

		end  % End wait for verifying BF
		
		% Set unit BF/Threshold to picked BF/Threshold
		BF_kHz{ind}=pickBF;
		Thresh_dBSPL{ind}=pickThr;
		
		legtext{ind} = sprintf('P%d:  BF=%.3f kHz; Thr=%.1f dB SPL; Q10=%.1f; SR=%.1f sps',PICind,BF_kHz{ind},Thresh_dBSPL{ind},Q10{ind},SR_sps);
		text(Xtext,.95-.05*(ind-1),legtext{ind},'Units','norm','Color',colors{mod(ind,7)+1})
	end
	
	
	
	
	Thresh_dBSPL_ret(ind) = Thresh_dBSPL{ind};
   BF_kHz_ret(ind) = BF_kHz{ind};
   Q10_ret(ind) = Q10{ind};

	if exist('h14','var')
		% 	set(h14,'Position',[485   612   914   358],'units','pixels')
		set(h14,'units','norm')
		set(h14,'Position',[0.2593    0.5800    0.6529    0.3410])
	end
end
if exist('h14','var')
	hold off
	orient landscape
end

   
return;






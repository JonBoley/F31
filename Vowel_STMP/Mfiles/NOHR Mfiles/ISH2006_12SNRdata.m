% File: ISH2006_2_04.m
% Modified from UnitLook_EHIN_CoincDet2.m (Apr 29, 2006) 
% function UnitLook_EHIN_CoincDet2(UnitName)
%
% cleans up figures for ISH paper.
%
% Verion 2: 21Apr2006 - improved for ISH paper.
%
% Date: 31Dec2005 - Modified to analyze SCC in terms of Coincidence
% detectors, rather than NSCCs.
%
% Date: 23Dec2005 (M. Heinz) (Modified from UnitPlot_SpatTemp_1SCC.m)
% For: R03 Experiments (ISH abstract - EHIN)
%
% This is specifically designed for EHINvNrBFi template used on 041805 and
% the conditions used there (4 noise attens, F1 and T1, 10 BF shifts)
%
% UnitName: '1.29' (converted later)
%
% Uses: setMonitor - if figures don't look right, re-run 'setMonitor'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% PRINTyes=input('Do you want to print figures automatically ([0]: no; 1: yes)?? ');
% if isempty(PRINTyes)
PRINTyes=0;
% end

doSCC=1;


global R03_dir MONITOR MONITORtext
global SavedPICS SavedPICnums SavedPICSuse
global FeaturesText FormsAtHarmonicsText InvertPolarityText
SavedPICSuse=1; SavedPICS=[]; SavedPICnums=[];
NOHR_dir=R03_dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ASSUME 04_18_05 is Experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExpDate='041805';  % HARDCODE
eval(['cd ''' fullfile(R03_dir,filesep,'ExpData',filesep,'MH-2005_04_18-ANnorm') ''''])
%%%% Find the full Experiment Name
[p,ExpName,e,v]=fileparts(pwd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir=fullfile(NOHR_dir,'ExpData',ExpName);
unitdata_dir=fullfile(data_dir,'UNITSdata');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
ISHfigs_dir=fullfile(data_dir,'UNITSdata','ISHfigs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HARD CODE 
UnitName='4.05';
% UnitName='1.04';
% beep;  
disp(sprintf('\n    **** HARD CODING Unit: %s ****',UnitName));

%%%% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));


% LOAD CALCS if available to speed up figure design
SAVECALCSfilename=sprintf('UnitLook_EHIN.%d.%02d.mat',TrackNum,UnitNum);
if exist(fullfile(unitdata_dir,SAVECALCSfilename),'file')
	% 	loadBOOL=input('Do you want to load the existing calculations for plotting (0: no; [1]: yes)?? ');
	loadBOOL=1;
	if loadBOOL
		disp(sprintf('Loading ''%s'' with all calculations already done!!',SAVECALCSfilename))
		thisPRINTyes=PRINTyes;
		eval(['load ''' fullfile(unitdata_dir,SAVECALCSfilename) ''''])
		PRINTyes=thisPRINTyes;
		setMonitor;  % Make sure you get the current monitor
	end
else
	error('DATA not stored, need to run UnitLook_EHIN_CoincDet2.m')
end

if ~loadBOOL
	%%% Load DataList for 041805
	global DataList
	if isempty(DataList)
		disp(' ... Loading DataList for 04_18-05');
		load DataList_2005_04_18
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%%%% Verify that there is data for this unit
	if isempty(DataList.Units{TrackNum,UnitNum})
		error('NO DATA FOR THIS UNIT!!');
	end
	Datafields=fieldnames(DataList.Units{TrackNum,UnitNum});
	if ~sum(strcmp(Datafields,{'EHrlv'}))
		error('NO ''EHrlv'' DATA FOR THIS UNIT, THEREFORE STOPPING!!');
	end

end
	
cd(data_dir)
disp(sprintf('Looking at Basic EHIN for Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify UNITSdata/unit analyses are all done ahead of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~loadBOOL
	%%%% Load unit structure for this unit
	UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
	eval(['ddd=dir(''' fullfile(unitdata_dir,UnitFileName) ''');'])
	% If UNITSdata file does not exist, load create from DataList, ow/ load it
	if isempty(ddd)
		unit=DataList.Units{TrackNum,UnitNum};
		eval(['save ''' fullfile(unitdata_dir,UnitFileName) ''' unit'])
	else
		eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
	end

	% If EHINvNreBFi_simFF analysis is not completed, run here (this will also
	% check EHINvNreBFi analysis and run if needed!)
	if ~isfield(unit,'EHvN_reBF_simFF')
		UnitAnal_EHvNrBF_simFF(ExpDate,UnitName,0);
		eval(['load ''' fullfile(unitdata_dir,UnitFileName) ''''])
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show EHrlfs for this unit
EHrlf_piclist=findPics('EHrlv',[TrackNum UnitNum]);
% quick_EHrlfs(EHrlf_piclist,CalibPic)
ISHplot_EHrlfs(EHrlf_piclist,CalibPic)

EHrlf_FIG=gcf;
set(EHrlf_FIG,'Name',sprintf('EHrlfs for Unit: %s',UnitName))
set(EHrlf_FIG,'NumberTitle','off')
YLIMITS_EHrlf=ylim;
% OLDtitle=get(get(gca,'Title'),'String');
% title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s', ...
% 	ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,OLDtitle), ...
% 	'units','norm')

% Label L_EH - tone level for rest of data

TextFontSize=18;
AnnotFonotSize=16;
LineWidth=3;

EHINrlf_piclist=findPics('EHINrlv',[TrackNum UnitNum]);
if ~isempty(EHINrlf_piclist)
	x=loadPic(EHINrlf_piclist(1));
	VOWELlevel_dBSPL=x.Stimuli.Condition.Level_dBSPL;
	plot(VOWELlevel_dBSPL*ones(2,1),[0 1000],'k:','LineWidth',LineWidth)
% 	text(VOWELlevel_dBSPL,0,'L_{EH}','HorizontalAlignment','center','VerticalAlignment','bottom')
	text(VOWELlevel_dBSPL+3,50,sprintf('Vowel Level for STMP = %.f dB SPL',VOWELlevel_dBSPL),'HorizontalAlignment','left','FontSize',TextFontSize,'units','data')
end
if strcmp(MONITOR,'laptop')
	set(EHrlf_FIG,'pos',[0   415   280   265],'Resize','off')
elseif strcmp(MONITOR,'docking')
	set(EHrlf_FIG,'pos',[0   420   263   265],'Resize','off')
else
	beep
	disp('Monitor unknown, custom figure placement not done (look at setMonitor.m)!')
end
clear x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show EHINrlfs for this unit
% quick_EHINrlfs(EHINrlf_piclist,CalibPic)
dBAtt_2_SNR = ISHplot_EHINrlfs(EHINrlf_piclist,CalibPic);

EHINrlf_FIG=gcf;
set(EHINrlf_FIG,'Name',sprintf('EHINrlfs for Unit: %s',UnitName))
set(EHINrlf_FIG,'NumberTitle','off')
YLIMITS_EHINrlf=ylim;
% OLDtitle=get(get(gca,'Title'),'String');
% title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n%s', ...
% 	ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10,OLDtitle), ...
% 	'units','norm')

% Label A_N - noise Atten for rest of data
EHvNreBFi_piclist=findPics('EHvNrBFi',[TrackNum UnitNum]);
if ~isempty(EHvNreBFi_piclist)
	x=loadPic(EHvNreBFi_piclist(1));
	SNR_dB=x.Stimuli.Condition.NoiseAttens_dB(3)+dBAtt_2_SNR;
	plot(SNR_dB*ones(2,1),[0 1000],'k:','LineWidth',LineWidth)
	plot(SNR_dB*ones(2,1)-10,[0 1000],'k:','LineWidth',LineWidth)
	plot(SNR_dB*ones(2,1)+10,[0 1000],'k:','LineWidth',LineWidth)
% 	text(SNR_dB,0,'SNR_{N}','HorizontalAlignment','center','VerticalAlignment','bottom')
	text(.05,.15,sprintf('Signal-to-Noise Ratios for STMP = %.f, %.f, %.f dB',SNR_dB+10,SNR_dB,SNR_dB-10),'HorizontalAlignment','left','units','norm','FontSize',TextFontSize)
end
if strcmp(MONITOR,'laptop')
	set(EHINrlf_FIG,'pos',[0   71   280   265],'Resize','off')
elseif strcmp(MONITOR,'docking')
	set(EHINrlf_FIG,'pos',[0   68   263   265],'Resize','off')
else
	beep
	disp('Monitor unknown, custom figure placement not done (look at setMonitor.m)!')
end
ymax=max([YLIMITS_EHrlf(2) YLIMITS_EHINrlf(2)]);
set(0,'CurrentFigure',EHrlf_FIG)
ylim([0 ymax])
set(0,'CurrentFigure',EHINrlf_FIG)
ylim([0 ymax])

if isempty(EHvNreBFi_piclist)
	return;
end
clear x




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL PERhist PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot data
ATTENcolors={'b','k','g','k','c','y'};
FEATmarkers={'.','x','s','^','*','<','^','>'};
FIG.FontSize=8;

TextFontSize=18;
AnnotFontSize=14;
LetterFontSize=20;
BASELineWidth=2;
NoiseLineStyle='-';
QuietLineStyle='--';

FeatureColors={'r','g'};
%%% Show multiple periods of tone PERhist
NperiodsTONE=5;
if isfield(unit,'EHvN_reBF_simFF')
	PERhist_XMAX=1/F0min*1000;  % take the PERhist_Xlims on vowels (in ms)
	XLIMITS_perhist=[0 PERhist_XMAX];
else
	% Need to find PERhist_XMAX by hand
	PERhist_XMAX=0;
	for ROWind=1:NUMrows
		for ATTind=1:length(Nattens_dB)
			for BFind=1:length(PERhistXs_sec{1,1})
				if max(PERhistXs_sec{ROWind,ATTind}{BFind})>PERhist_XMAX
					PERhist_XMAX=max(PERhistXs_sec{ROWind,ATTind}{BFind});
				end
			end
		end
	end
	PERhist_XMAX=PERhist_XMAX*1000; % in ms
	XLIMITS_perhist=[0 PERhist_XMAX*NperiodsTONE];
end
XLIMITS_rate=[0 300];
XLIMITS_synch=[0 1];
% XLIMITS_phase=[-pi pi];
%% Find PERhist_YMAX
PERhistGAIN=2.0; % # of channels covered by max PERhist
PERhists_logCHwidth=log10(highBF/lowBF)/(PERhistsYCHANS-1);  % log10 channel width
PERhist_YMIN=lowBF;
PERhist_YMAX=10^((PERhistsYCHANS-1+PERhistGAIN)*PERhists_logCHwidth)*lowBF;    % sets an extra (GAIN-1) log channel widths
YLIMITS=[PERhist_YMIN PERhist_YMAX];  % Used for all plots
%% This  is ALL needed to get the right LOG Yticks!!
YLIMunit=10^floor(log10(lowBF));
YLIMS=floor(lowBF/YLIMunit)*YLIMunit*[1 100]; % Do two decades to be sure we get all the ticks
YTICKS=[YLIMS(1):YLIMunit:YLIMS(1)*10 YLIMS(1)*20:YLIMunit*10:YLIMS(end)];
BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave

% % Setup parameters for title
% if isempty(unit.Info.Threshold_dBSPL)
% 	unit.Info.Threshold_dBSPL=NaN;
% end
% if isempty(unit.Info.SR_sps)
% 	unit.Info.SR_sps=NaN;
% end
% if isempty(unit.Info.Q10)
% 	unit.Info.Q10=NaN;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADD EXTRA PERhist plot at the end with all ATTENs on top of one anothe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1000); clf
PERhist_FIG=gcf;
set(PERhist_FIG,'Name',sprintf('PERhists for Unit: %s',UnitName))
set(PERhist_FIG,'NumberTitle','off')

if strcmp(MONITOR,'laptop')
	set(gcf,'pos',[287   730   550   240],'Resize','off')
elseif strcmp(MONITOR,'docking')
	set(gcf,'pos',[270   709   472   235],'Resize','off')   
else
	beep
	disp('Monitor unknown, custom figure placement not done (look at setMonitor.m)!')
end
ROWind=0;
ALLlevelsTriFilt=3;

ATTEN=max(Nattens_dB);

%%%%%%%%%%%%%%%%%%%% EH_reBF Plots
if isfield(unit,'EHvN_reBF_simFF')
	for FeatIND=FeatINDs
		ROWind=ROWind+1;
		eval(['yTEMP=unit.EHvN_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
		if ~isempty(yTEMP)
			%%%% EH_reBF plots
			ATTind=find(yTEMP.Nattens_dB==ATTEN);

			%%%% Spatio-Temporal Plots
			PLOTnum=(ROWind-1)*NUMcols+1;
			eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
			LEGtext='';
			for BFind=1:length(BFs_kHz{ROWind,ATTind})
				if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
					LINEwidth=2*BASELineWidth;
				else
					LINEwidth=BASELineWidth;
				end
				% This normalization plots each signal the same size on a log scale
				if ~isempty(PERhistXs_sec{ROWind,ATTind}{BFind})
					NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,ATTind}(BFind)/PERhistsMAX;
					semilogy(PERhistXs_sec{ROWind,ATTind}{BFind}*1000, ...
						trifilt(PERhists_Smoothed{ROWind,ATTind}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,ATTind}(BFind), ...
						'LineWidth',LINEwidth,'Color',ATTENcolors{ATTind},'LineStyle',QuietLineStyle)
					hold on
					if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
						LEGtext{length(LEGtext)+1}='In Quiet';
					end
% 					for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
					for ATTind2=2
						if ~isempty(PERhistXs_sec{ROWind,ATTind2}{BFind})
							semilogy(PERhistXs_sec{ROWind,ATTind2}{BFind}*1000, ...
								trifilt(PERhists_Smoothed{ROWind,ATTind2}{BFind},ALLlevelsTriFilt)*NormFact+BFs_kHz{ROWind,ATTind2}(BFind), ...
								'LineWidth',LINEwidth,'Color',ATTENcolors{ATTind2},'LineStyle',NoiseLineStyle)
							if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,ATTind}/unit.Info.BF_kHz))<BFoctCRIT))
								LEGtext{length(LEGtext)+1}=sprintf('SNR = %.f dB',Nattens_dB(ATTind2)+dBAtt_2_SNR);
							end
						end
					end
					if strcmp(FeaturesText{FeatIND},'F1')
						hleg1000=legend(LEGtext,1);
						set(hleg1000,'FontSize',8)
						set(hleg1000,'pos',[0.4    0.87    0.0942    0.0473])
					end
				end
			end
			if ROWind==2
				xlabel('Time (ms)','FontSize',TextFontSize)
			end
			ylabel('Effective Best Frequency (kHz)','FontSize',TextFontSize)
			set(gca,'FontSize',AnnotFontSize)
			if ROWind==1
				text(.05,.9,'F1','units','norm','FontSize',LetterFontSize)
			else
				text(.05,.9,'T1','units','norm','FontSize',LetterFontSize)
			end
			xlim(XLIMITS_perhist)
			PLOThand=eval(['h' num2str(PLOTnum)]);
			set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
			ylim(YLIMITS)  % Same Ylimits for all plots
			%%%%%%%%%%%%%%%%%%%%%
			% Plot lines at all features
			for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
				if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
% 					semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
					text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000, ...
						sprintf('%s',FeaturesText{FeatINDPlot}), ...
						'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2},'FontSize',TextFontSize)
				end
% 				for BFind=1:length(BFs_kHz{ROWind,ATTind})
					if ~isempty(PERhistXs_sec{ROWind,ATTind}{BFind})
						if (FeatINDPlot<=FeatIND)
							if (FeatINDPlot>1)
								text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
									'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2},'FontSize',AnnotFontSize)
							else
								text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
									'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k','FontSize',AnnotFontSize)
							end
						end
					end
% 				end
			end
			hold off


			%%%% Rate Plot
			PLOTnum=(ROWind-1)*NUMcols+2;
			eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
			semilogy(Rates{ROWind,ATTind},BFs_kHz{ROWind,ATTind},'*-','Color',ATTENcolors{ATTind},'LineStyle',QuietLineStyle,'LineWidth',BASELineWidth)
			hold on
			%                semilogy(Nsps{ROWind,ATTind}/10,BFs_kHz{ROWind,ATTind},'m+','MarkerSize',4,'Color',ATTENcolors{ATTind})
% 			semilogy(ALSRs{ROWind,ATTind},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',ATTENcolors{ATTind})
% 			for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
			for ATTind2=2
				semilogy(Rates{ROWind,ATTind2},BFs_kHz{ROWind,ATTind2},'*-','Color',ATTENcolors{ATTind2},'LineStyle',NoiseLineStyle,'LineWidth',BASELineWidth)
				%                   semilogy(Nsps{ROWind,ATTind2}/10,BFs_kHz{ROWind,ATTind2},'m+','MarkerSize',4,'Color',ATTENcolors{ATTind2})
% 				semilogy(ALSRs{ROWind,ATTind2},unit.Info.BF_kHz,'go','MarkerSize',6,'Color',ATTENcolors{ATTind2})
			end
% 			semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
			if ROWind==2
				xlabel(sprintf('Rate (sp/sec)'),'FontSize',TextFontSize)
			end
			PLOThand=eval(['h' num2str(PLOTnum)]);
			xlim(XLIMITS_rate)
			set(gca,'FontSize',AnnotFontSize)
			set(PLOThand,'XDir','reverse')
			set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
			ylim(YLIMITS)  % Same Ylimits for all plots
			%%%%%%%%%%%%%%%%%%%%%
			% Plot lines at all features
			for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
				if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
					semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
				end
			end
			hold off

			%%%% Synch Plot
			PLOTnum=(ROWind-1)*NUMcols+3;
			eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
			semilogy(Synchs{ROWind,ATTind},BFs_kHz{ROWind,ATTind},'*-','Color',ATTENcolors{ATTind},'LineStyle',QuietLineStyle,'LineWidth',BASELineWidth)
			hold on
% 			for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
			for ATTind2=2
				semilogy(Synchs{ROWind,ATTind2},BFs_kHz{ROWind,ATTind2},'*-','Color',ATTENcolors{ATTind2},'LineStyle',NoiseLineStyle,'LineWidth',BASELineWidth)
			end
			semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
			if ROWind==2
				xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}),'FontSize',TextFontSize)
			end
			PLOThand=eval(['h' num2str(PLOTnum)]);
			xlim(XLIMITS_synch)
			set(PLOThand,'XDir','reverse')
			set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
			set(gca,'FontSize',AnnotFontSize)
			set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
			ylim(YLIMITS)  % Same Ylimits for all plots
			%%%%%%%%%%%%%%%%%%%%%
			% Plot lines at all features
			for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
				if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
					semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
				end
			end
			hold off

		end %End if data for this condition, plot
	end % End Feature
end


Xcorner=0.1;
Xwidth1=.5;
Xshift1=0.05;
Xwidth2=.1;
Xshift2=0.03;

Ycorner=0.075;
Yshift=0.07;
Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

TICKlength=0.02;

if NUMrows>4
	set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
	set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% 	set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>3
	set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
	set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% 	set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>2
	set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
	set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% 	set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

if NUMrows>1
	set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
	set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% 	set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
end

set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
% set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

orient landscape


if doSCC
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% DO ALL SCC PLOTS
	% 4/27/06 - Plot SCC on sp/sec scale, give up on POP plotting
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	XLIMITS_sccF1=1.5*1000/unit.EHvN_reBF_simFF.F1{1,1}.FeatureFreqs_Hz{1}(1)*[-1 1];
	XLIMITS_sccT1=1.5*1000/unit.EHvN_reBF_simFF.T1{1,1}.FeatureFreqs_Hz{1}(1)*[-1 1];
	YLIMITS_scc=[0 9];
	

	TextFontSize=18;
	AnnotFontSize=14;
	LetterFontSize=20;
	BASELineWidth=2;
	ZEROLineWidth=2;
	NoiseLineStyle='-'
	QuietLineStyle='--'


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%% EXTRA SCC PLOT WITH ALL ATTENS
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	ATTEN=max(Nattens_dB);
	figure(1002); clf
	SCC_fig=gcf;
	set(SCC_fig,'Name',sprintf('SCCs for Unit: %s',UnitName))
	set(SCC_fig,'NumberTitle','off')
	if strcmp(MONITOR,'laptop')
		set(gcf,'pos',[287    72   550   240],'Resize','off')
	elseif strcmp(MONITOR,'docking')
		set(gcf,'pos',[270    68   472   235],'Resize','off')
	else
		beep
		disp('Monitor unknown, custom figure placement not done (look at setMonitor.m)!')
	end
	ROWind=0;

	%%%%%%%%%%%%%%%%%%%% EH_reBF Plots
	if isfield(unit,'EHvN_reBF_simFF')
		for FeatIND=FeatINDs
			ROWind=ROWind+1;
			eval(['yTEMP=unit.EHvN_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
			if ~isempty(yTEMP)
				%%%% EH_reBF plots
% 				ATTind=find(yTEMP.Nattens_dB==ATTEN);     % IN QUIET
				ATTind=2;     % middle SNR

				%%%% Spatio-Temporal Plots
				PLOTnum=(ROWind-1)*NUMcols+1;
				eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
				LEGtext='';
				for SCCind=1:length(NSCCs{ROWind,ATTind})
					LINEwidth=.5;
					if ~isempty(NSCCs{ROWind,ATTind}{SCCind})
						plot(NSCC_delays_usec{ROWind,ATTind}{SCCind}/1000, ...
							NSCCs{ROWind,ATTind}{SCCind},'LineWidth',BASELineWidth,'Color',ATTENcolors{ATTind})
						hold on
						if SCCind==1
							LEGtext{length(LEGtext)+1}=sprintf('In Quiet',Nattens_dB(ATTind)+dBAtt_2_SNR);
						end
					end
% 					for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
					for ATTind2=[]
						if ~isempty(NSCCs{ROWind,ATTind2}{SCCind})
							plot(NSCC_delays_usec{ROWind,ATTind2}{SCCind}/1000, ...
								NSCCs{ROWind,ATTind2}{SCCind},'LineWidth',BASELineWidth,'Color',ATTENcolors{ATTind2})
							hold on
							if SCCind==1
								LEGtext{length(LEGtext)+1}=sprintf('SNR = %.f dB',Nattens_dB(ATTind2)+dBAtt_2_SNR);
							end
						end
					end
					if SCCind==1
						hleg1001=legend(LEGtext,1);
						set(hleg1001,'FontSize',8)
						set(hleg1001,'pos',[0.400    0.4    0.0942    0.0473])
					end
					if ~isempty(NSCCs{ROWind,ATTind}{SCCind})
% 						% Plot 0 and 1 lines
% 						semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,ATTind}{SCCind}), ...
% 							'k-','LineWidth',LINEwidth/2)
% 						semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,ATTind}{SCCind}), ...
% 							'k-','LineWidth',LINEwidth/2)
						plot(zeros(1,2),YLIMITS_scc,'k:','LineWidth',ZEROLineWidth)
						%                      % Plot values at 0 delay and CD
% 						plot(0,NSCC_0delay{ROWind,ATTind}{SCCind},'s','Color',ATTENcolors{ATTind})
% 						plot(NSCC_CDs_usec{ROWind,ATTind}{SCCind}/1000,NSCC_peaks{ROWind,ATTind}{SCCind}, ...
% 							'o','Color',ATTENcolors{ATTind})
					end
% 					for ATTind2=fliplr(find(Nattens_dB~=max(Nattens_dB)))
					for ATTind2=[]
% 						% Plot 0 and 1 lines
% 						semilogy(XLIMITS_scc,zeros(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,ATTind2}{SCCind}), ...
% 							'k-','LineWidth',LINEwidth/2)
% 						semilogy(XLIMITS_scc,ones(1,2)*NormFact+geomean(NSCC_BFskHz{ROWind,ATTind2}{SCCind}), ...
% 							'k-','LineWidth',LINEwidth/2)
% 						semilogy(zeros(1,2),YLIMITS,'k:','LineWidth',LINEwidth/2)
						%             % Plot values at 0 delay and CD
% 						plot(0,NSCC_0delay{ROWind,ATTind2}{SCCind},'s','Color',ATTENcolors{ATTind2})
% 						plot(NSCC_CDs_usec{ROWind,ATTind2}{SCCind}/1000,NSCC_peaks{ROWind,ATTind2}{SCCind}, ...
% 							'o','Color',ATTENcolors{ATTind2})
					end

				end
				xlabel('Delay (ms)','FontSize',TextFontSize)
				ylabel('Coincidence Rate (sp/sec)','FontSize',TextFontSize)
				if ROWind==1
					text(.05,.9,'F1','units','norm','FontSize',LetterFontSize)
				else
					text(.05,.9,'T1','units','norm','FontSize',LetterFontSize)
				end
				set(gca,'FontSize',AnnotFontSize)
				if ROWind==1
					xlim(XLIMITS_sccF1)
				else
					xlim(XLIMITS_sccT1)
				end
				PLOThand=eval(['h' num2str(PLOTnum)]);
				ylim(YLIMITS_scc)  % Same Ylimits for all plots
				%%%%%%%%%%%%%%%%%%%%%
				% Plot lines at all features
				for FeatINDPlot=find(~strcmp(FeaturesText,'TN'))
					%                      if (yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000<=YLIMITS(2))
					%                         semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
					%                         text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000, ...
					%                            sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{ATTind}(FeatINDPlot)/1000), ...
					%                            'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
					%                      end
					for SCCind=1:length(NSCCs{ROWind,ATTind})
						if ~isempty(NSCCs{ROWind,ATTind}{SCCind})
							if (FeatINDPlot<=FeatIND)
								if (FeatINDPlot>1)
									text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS_scc(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
										'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2},'FontSize',AnnotFontSize)
									plot(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot)*ones(1,2),[0 .75],'-k','Linewidth',2)
								else
									text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS_scc(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
										'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k','FontSize',AnnotFontSize)
									plot(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot)*ones(1,2),[0 .75],'-k','Linewidth',2)
								end
							end
						end
					end
					XLIMITS_scc=1.5*1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot);
				end
				hold off

			end %End if data for this condition, plot
		end % End Feature
	end


	Xcorner=0.05;
	Xwidth1=.5;
	Xshift1=0.05;
	Xwidth2=.1;
	Xshift2=0.03;

	Ycorner=0.05;
	Yshift=0.07;
	Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2

	TICKlength=0.02;

	if NUMrows>4
		set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	end

	if NUMrows>3
		set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	end

	if NUMrows>2
		set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	end

	if NUMrows>1
		set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
	end

	set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])

	orient landscape

end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO SMP PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SMPlegFontSize=7;
FIGnum=1;
figure(FIGnum); clf
SMP_fig=gcf;
set(SMP_fig,'Name',sprintf('SMP for Unit: %s',UnitName))
set(SMP_fig,'NumberTitle','off')

if strcmp(MONITOR,'laptop')
	set(gcf,'pos',[842    74   562   898],'Resize','off')
elseif strcmp(MONITOR,'docking')
	set(gcf,'pos',[749    69   528   873],'Resize','off')
else
	beep
	disp('Monitor unknown, custom figure placement not done (look at setMonitor.m)!')
end

TextFontSize=18;
AnnotFontSize=14;
LetterFontSize=20;
SMPslopeLineWidth=2;
ZEROLineWidth=2;
NoiseLineStyle='-'
QuietLineStyle='--'
SMPslopeMarkerSize=12;


%%%%%%%%%%%%%%%%%%%%% CALCULATE SLOPES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rate vs. Feature Level
SensitivitySlopes_rate=NaN+ones(size(Nattens_dB));
SensitivityIntercepts_rate=NaN+ones(size(Nattens_dB));
for ATTind=1:length(Nattens_dB)
	if sum(~isnan(SMP_rate{ATTind}))>1
		x=yTEMP.levels_dBSPL+FeatureLevels_dB(find(~isnan(SMP_rate{ATTind})));
		y=SMP_rate{ATTind}(find(~isnan(SMP_rate{ATTind})));
		[Cfit,MSE,fit]=fit1slope(x,y);
		SensitivitySlopes_rate(ATTind)=Cfit(1);
		SensitivityIntercepts_rate(ATTind)=Cfit(2);
	end
end

% ALSR vs. Feature Level
SensitivitySlopes_alsr=NaN+ones(size(Nattens_dB));
SensitivityIntercepts_alsr=NaN+ones(size(Nattens_dB));
for ATTind=1:length(Nattens_dB)
	if sum(~isnan(SMP_alsr{ATTind}))>1
		x=yTEMP.levels_dBSPL+FeatureLevels_dB(find(~isnan(SMP_alsr{ATTind})));
		y=SMP_alsr{ATTind}(find(~isnan(SMP_alsr{ATTind})));
		[Cfit,MSE,fit]=fit1slope(x,y);
		SensitivitySlopes_alsr(ATTind)=Cfit(1);
		SensitivityIntercepts_alsr(ATTind)=Cfit(2);
	end
end

% NSCC_peak(1) vs. Feature Level
SensitivitySlopes_nsccPEAK{1}=NaN+ones(size(Nattens_dB));
SensitivityIntercepts_nsccPEAK{1}=NaN+ones(size(Nattens_dB));
for ATTind=1:length(Nattens_dB)
	if sum(~isnan(SMP_NSCC_peak{1,ATTind}))>1
		x=yTEMP.levels_dBSPL+FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,ATTind})));
		y=SMP_NSCC_peak{1,ATTind}(find(~isnan(SMP_NSCC_peak{1,ATTind})));
		[Cfit,MSE,fit]=fit1slope(x,y);
		SensitivitySlopes_nsccPEAK{1}(ATTind)=Cfit(1);
		SensitivityIntercepts_nsccPEAK{1}(ATTind)=Cfit(2);
	end
end


%%%
SMPLineWidths=[1 2 3 4];
SMPsymbols = {'x','o','^','s'};
SMPMarkerSize=12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rate-FeatureLevels for each OAL
hSMP1=subplot(521);
LEGtext='';
for ATTind=1:length(Nattens_dB)
	%% Plot fitted lines
	if ~isnan(SensitivitySlopes_rate(ATTind))
		xdata=[min(FeatureLevels_dB(find(~isnan(SMP_rate{ATTind})))+yTEMP.levels_dBSPL) max(FeatureLevels_dB(find(~isnan(SMP_rate{ATTind})))+yTEMP.levels_dBSPL)];
		ydata=SensitivitySlopes_rate(ATTind)*xdata+SensitivityIntercepts_rate(ATTind);
		plot(xdata,ydata,'k-','LineWidth',SMPLineWidths(ATTind),'Marker',SMPsymbols{ATTind},'MarkerSize',SMPMarkerSize)
		hold on
% 		LEGtext{length(LEGtext)+1}=sprintf('%.f dB (%.2f)',Nattens_dB(ATTind)+dBAtt_2_SNR,SensitivitySlopes_rate(ATTind)); %% sp/sec/dB
	end
end
% for ATTind=1:length(Nattens_dB)
% 	for FeatIND=FeatINDs
% 		if ~isnan(SMP_rate{ATTind}(FeatIND))
% % 			plot(yTEMP.levels_dBSPL+FeatureLevels_dB(FeatIND),SMP_rate{ATTind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',ATTENcolors{ATTind},'LineStyle','none')
% 			hold on
% 		end
% 	end
% end
YLIMITS_SMPrate=[0 300];
XLIMITS_SMPrate=[15 55];
ylim(YLIMITS_SMPrate)  % Fixed ordinate for all plots
xlim(XLIMITS_SMPrate)
%%%%%%%%%%%%%%%%%%%%%
% Home-made Feature Symbol Legend
%%%%%%%%%%%%%%%%%%%%%
% LEGXleft=.8; LEGYbottom=0.15; LEGYstep=0.1; LEGXsymbOFFset=0.05;
% FeatNum=0;
% for FeatIND=fliplr(FeatINDs)
% 	FeatNum=FeatNum+1;
% 	plot(XLIMITS_SMPrate(1)+diff(XLIMITS_SMPrate)*(LEGXleft),YLIMITS_SMPrate(1)+diff(YLIMITS_SMPrate)*(LEGYbottom+(FeatNum-1)*LEGYstep),FEATmarkers{FeatIND}, ...
% 		'Color','k','MarkerSize',6)
% 	text(LEGXleft+LEGXsymbOFFset,LEGYbottom+(FeatNum-1)*LEGYstep,FeaturesText{FeatIND},'Units','norm','FontSize',10)
% end
ylabel('Rate (sp/sec)','FontSize',TextFontSize)
% xlabel('Feature Level (dB SPL)')
% hleg=legend(LEGtext,1);
% set(hleg,'FontSize',SMPlegFontSize)
hold off
set(gca,'FontSize',AnnotFontSize)
% set(gca,'XTickLabel','')
% title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
% 	ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10), ...
% 	'units','norm')


%%%% ALSR-FeatureLevels for each OAL
hSMP3=subplot(523);
LEGtext='';
for ATTind=1:length(Nattens_dB)
	%% Plot fitted lines
	if ~isnan(SensitivitySlopes_alsr(ATTind))
		xdata=[min(FeatureLevels_dB(find(~isnan(SMP_alsr{ATTind})))+yTEMP.levels_dBSPL) max(FeatureLevels_dB(find(~isnan(SMP_alsr{ATTind})))+yTEMP.levels_dBSPL)];
		ydata=SensitivitySlopes_alsr(ATTind)*xdata+SensitivityIntercepts_alsr(ATTind);
% 		plot(xdata,ydata,'Marker','none','Color',ATTENcolors{ATTind},'LineStyle','-')
		plot(xdata,ydata,'k-','LineWidth',SMPLineWidths(ATTind),'Marker',SMPsymbols{ATTind},'MarkerSize',SMPMarkerSize)
		hold on
		LEGtext{length(LEGtext)+1}=sprintf('%.f dB (%.2f)',Nattens_dB(ATTind)+dBAtt_2_SNR,SensitivitySlopes_alsr(ATTind));  % sp/sec/dB
	end
end
% for ATTind=1:length(Nattens_dB)
% 	for FeatIND=FeatINDs
% 		if ~isnan(SMP_alsr{ATTind}(FeatIND))
% 			plot(yTEMP.levels_dBSPL+FeatureLevels_dB(FeatIND),SMP_alsr{ATTind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',ATTENcolors{ATTind},'LineStyle','none')
% 			hold on
% 		end
% 	end
% end
YLIMITS_SMPalsr=[0 300];
XLIMITS_SMPalsr=XLIMITS_SMPrate;
ylim(YLIMITS_SMPalsr)  % Fixed ordinate for all plots
xlim(XLIMITS_SMPalsr)
ylabel('ALSR (sp/sec)','FontSize',TextFontSize)
% xlabel('Feature Level (dB SPL)')
% hleg=legend(LEGtext,1);
% set(hleg,'FontSize',SMPlegFontSize)
hold off
% set(gca,'XTickLabel','')
set(gca,'FontSize',AnnotFontSize)


%%%% NSCC_peak-FeatureLevels for each OAL
%%%% [This is the Deng and Geisler Analysis!!!]
hSMP5=subplot(525);
LEGtext='';
for ATTind=1:length(Nattens_dB)
	%% Plot fitted lines
	if ~isnan(SensitivitySlopes_nsccPEAK{1}(ATTind))
		xdata=[min(FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,ATTind})))+yTEMP.levels_dBSPL) max(FeatureLevels_dB(find(~isnan(SMP_NSCC_peak{1,ATTind})))+yTEMP.levels_dBSPL)];
		ydata=SensitivitySlopes_nsccPEAK{1}(ATTind)*xdata+SensitivityIntercepts_nsccPEAK{1}(ATTind);
		plot(xdata,ydata,'k-','LineWidth',SMPLineWidths(ATTind),'Marker',SMPsymbols{ATTind},'MarkerSize',SMPMarkerSize)
% 		plot(xdata,ydata,'Marker','none','Color',ATTENcolors{ATTind},'LineStyle','-')
		hold on
		if ATTind~=4
			LEGtext{length(LEGtext)+1}=sprintf('%.f dB',Nattens_dB(ATTind)+dBAtt_2_SNR);  % sp/sec/dB SCC_CD
		else
			LEGtext{length(LEGtext)+1}=sprintf('In Quiet');  
		end
	end
end
% for ATTind=1:length(Nattens_dB)
% 	for FeatIND=FeatINDs
% 		if ~isnan(SMP_NSCC_peak{1,ATTind}(FeatIND))
% 			plot(yTEMP.levels_dBSPL+FeatureLevels_dB(FeatIND),SMP_NSCC_peak{1,ATTind}(FeatIND),'Marker',FEATmarkers{FeatIND},'Color',ATTENcolors{ATTind},'LineStyle','none')
% 			hold on
% 		end
% 	end
% end
YLIMITS_SMPnsccPEAK=[0 10];
XLIMITS_SMPnsccPEAK=XLIMITS_SMPrate;
ylim(YLIMITS_SMPnsccPEAK)  % Fixed ordinate for all plots
xlim(XLIMITS_SMPnsccPEAK)
ylabel(sprintf('Max Coinc. Rate\n(sp/sec)'),'FontSize',TextFontSize)
xlabel('Vowel Feature Level (dB SPL)','FontSize',TextFontSize)
hleg=legend(LEGtext,4);
set(hleg,'FontSize',SMPlegFontSize)
hold off
set(gca,'FontSize',AnnotFontSize)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO SUMMARY SMP PLOT versus Natten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SMPlegFontSize=7;

NOnoiseIND=find(Nattens_dB==120);
noiseINDs=setdiff(1:length(Nattens_dB),NOnoiseIND);
LOWnoiseIND=find(Nattens_dB==max(Nattens_dB(noiseINDs)));
xmax=Nattens_dB(LOWnoiseIND)+dBAtt_2_SNR+8; xmin=min(Nattens_dB(noiseINDs)+dBAtt_2_SNR);
xmin=-8;
ymin=0; ymax=6;
yminFACT=-.4; ymaxFACT=1.7;

% RATE
hSMP2=subplot(522);
plot(Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_rate(noiseINDs),'b-x','LineWidth',SMPslopeLineWidth,'MarkerSize',SMPslopeMarkerSize)
xlim([xmin-2 xmax+2])
hold on
plot(20,SensitivitySlopes_rate(NOnoiseIND),'rs','MarkerFaceColor','r','MarkerSize',SMPslopeMarkerSize-2)

% plot([xmax Nattens_dB(LOWnoiseIND)+dBAtt_2_SNR],SensitivitySlopes_rate([NOnoiseIND LOWnoiseIND]),'r--')
plot([xmin-2 xmax+2],[0 0],'k--','LineWidth',ZEROLineWidth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustrate 1/2 point SNR for degraded spectral coding
plot([xmin-2 20],SensitivitySlopes_rate(NOnoiseIND)*[.5 .5],'k:','LineWidth',ZEROLineWidth)
HALFpt_rate_SNR = interp1(SensitivitySlopes_rate(noiseINDs),Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_rate(NOnoiseIND)*.5);
HALFpt_rate_NOTexact=0;
if isnan(HALFpt_rate_SNR)&sum(SensitivitySlopes_rate(noiseINDs)>(SensitivitySlopes_rate(NOnoiseIND)*.5))==3
	%% All above cutoff
	HALFpt_rate_SNR= Nattens_dB(1)+dBAtt_2_SNR;    % set to lowest SNR measured
	HALFpt_rate_NOTexact=1;
elseif isnan(HALFpt_rate_SNR)&sum(SensitivitySlopes_rate(noiseINDs)<(SensitivitySlopes_alsr(NOnoiseIND)*.5))==3
	%% All below cutoff
	HALFpt_rate_SNR= Nattens_dB(3)+dBAtt_2_SNR;  % set to highest SNR measured
	HALFpt_rate_NOTexact=-1;
end	
plot(HALFpt_rate_SNR*[1 1],SensitivitySlopes_rate(NOnoiseIND)*[-10 .5],'k:','LineWidth',ZEROLineWidth)
disp(sprintf('RATE: HALF-point SNR = %.1f dB (%d)',HALFpt_rate_SNR,HALFpt_rate_NOTexact))

hold off
YLIMITS_rA=[-1 6];
ylim(YLIMITS_rA)
% ylim([yminFACT ymaxFACT]*SensitivitySlopes_rate(NOnoiseIND))
text(20,0,'IQ','units','data','HorizontalAlignment','center','VerticalAlignment','bottom','color','red')
% title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\nRATE', ...
% 	ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10), ...
% 	'units','norm')
% title('RATE')
% ylabel('SMP Slope (sp/sec/dB)')
% xlabel('Noise Attenuation (dB)')
% set(gca,'XTickLabel','')
set(gca,'XDir','reverse')
set(gca,'FontSize',AnnotFontSize)
xlim([-20 30])

% ALSR
hSMP4=subplot(524);
plot(Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_alsr(noiseINDs),'b-x','LineWidth',SMPslopeLineWidth,'MarkerSize',SMPslopeMarkerSize)
xlim([xmin-2 xmax+2])
hold on
plot(20,SensitivitySlopes_alsr(NOnoiseIND),'rs','MarkerFaceColor','r','MarkerSize',SMPslopeMarkerSize-2)
% plot([xmax Nattens_dB(LOWnoiseIND)+dBAtt_2_SNR],SensitivitySlopes_alsr([NOnoiseIND LOWnoiseIND]),'r--')
plot([xmin-2 xmax+2],[0 0],'k--','LineWidth',ZEROLineWidth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustrate 1/2 point SNR for degraded spectral coding
plot([xmin-2 20],SensitivitySlopes_alsr(NOnoiseIND)*[.5 .5],'k:','LineWidth',ZEROLineWidth)
HALFpt_alsr_SNR = interp1(SensitivitySlopes_alsr(noiseINDs),Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_alsr(NOnoiseIND)*.5);
HALFpt_alsr_NOTexact=0;
if isnan(HALFpt_alsr_SNR)&sum(SensitivitySlopes_alsr(noiseINDs)>(SensitivitySlopes_alsr(NOnoiseIND)*.5))==3
	%% All above cutoff
	HALFpt_alsr_SNR= Nattens_dB(1)+dBAtt_2_SNR;    % set to lowest SNR measured
	HALFpt_alsr_NOTexact=1;
elseif isnan(HALFpt_alsr_SNR)&sum(SensitivitySlopes_alsr(noiseINDs)<(SensitivitySlopes_alsr(NOnoiseIND)*.5))==3
	%% All below cutoff
	HALFpt_alsr_SNR= Nattens_dB(3)+dBAtt_2_SNR;  % set to highest SNR measured
	HALFpt_alsr_NOTexact=-1;
end	
plot(HALFpt_alsr_SNR*[1 1],SensitivitySlopes_alsr(NOnoiseIND)*[-10 .5],'k:','LineWidth',ZEROLineWidth)
disp(sprintf('ALSR: HALF-point SNR = %.1f dB (%d)',HALFpt_alsr_SNR,HALFpt_alsr_NOTexact))


hold off
ylim([YLIMITS_rA])
% ylim([yminFACT ymaxFACT]*SensitivitySlopes_alsr(NOnoiseIND))
text(20,0,'IQ','units','data','HorizontalAlignment','center','VerticalAlignment','bottom','color','red')
% title('ALSR')
ylabel('SMP Slope (sp/sec/dB)','FontSize',TextFontSize)
% xlabel('Noise Attenuation (dB)')
% set(gca,'XTickLabel','')
set(gca,'XDir','reverse')
set(gca,'FontSize',AnnotFontSize)
xlim([-20 30])

if strcmp(UnitName,'4.06')
	SensitivitySlopes_nsccPEAK{1}=[.0063 .0033 .0223 .0297];
elseif strcmp(UnitName,'4.11')
	SensitivitySlopes_nsccPEAK{1}=[-.0017 .0053 .0187 .0370];
elseif strcmp(UnitName,'4.05')
	SensitivitySlopes_nsccPEAK{1}=[0.0057 .0460 .0517 .0577];
end


% Deng and Geisler 1987
% yminDG=-6e-1; ymaxDG=0;
hSMP6=subplot(526);
plot(Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_nsccPEAK{1}(noiseINDs),'b-x','LineWidth',SMPslopeLineWidth,'MarkerSize',SMPslopeMarkerSize)
xlim([xmin-2 xmax+2])
hold on
plot(20,SensitivitySlopes_nsccPEAK{1}(NOnoiseIND),'rs','MarkerFaceColor','r','MarkerSize',SMPslopeMarkerSize-2)
% plot([xmax Nattens_dB(LOWnoiseIND)+dBAtt_2_SNR],SensitivitySlopes_nsccPEAK{1}([NOnoiseIND LOWnoiseIND]),'r--')
plot([xmin-2 xmax+2],[0 0],'k--','LineWidth',ZEROLineWidth)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustrate 1/2 point SNR for degraded spectral coding

%%% NEEDED TO QUICKLY GET INTERPOLATION TO WORK RIGHT for this condition
if strcmp(UnitName,'2.04')
	SensitivitySlopes_nsccPEAK{1}(1)=-999;
elseif strcmp(UnitName,'4.06')
	SensitivitySlopes_nsccPEAK{1}(1)=-999;
end

plot([xmin-2 20],SensitivitySlopes_nsccPEAK{1}(NOnoiseIND)*[.5 .5],'k:','LineWidth',ZEROLineWidth)
HALFpt_DG_SNR = interp1(SensitivitySlopes_nsccPEAK{1}(noiseINDs),Nattens_dB(noiseINDs)+dBAtt_2_SNR,SensitivitySlopes_nsccPEAK{1}(NOnoiseIND)*.5);
HALFpt_DG_NOTexact=0;
if isnan(HALFpt_DG_SNR)&sum(SensitivitySlopes_nsccPEAK{1}(noiseINDs)>(SensitivitySlopes_nsccPEAK{1}(NOnoiseIND)*.5))==3
	%% All above cutoff
	HALFpt_DG_SNR= Nattens_dB(1)+dBAtt_2_SNR;    % set to lowest SNR measured
	HALFpt_DG_NOTexact=1;
elseif isnan(HALFpt_DG_SNR)&sum(SensitivitySlopes_nsccPEAK{1}(noiseINDs)<(SensitivitySlopes_nsccPEAK{1}(NOnoiseIND)*.5))==3
	%% All below cutoff
	HALFpt_DG_SNR= Nattens_dB(3)+dBAtt_2_SNR;  % set to highest SNR measured
	HALFpt_DG_NOTexact=-1;
end	
plot(HALFpt_DG_SNR*[1 1],SensitivitySlopes_nsccPEAK{1}(NOnoiseIND)*[-10 .5],'k:','LineWidth',ZEROLineWidth)
disp(sprintf('DG: HALF-point SNR = %.1f dB (%d)',HALFpt_DG_SNR,HALFpt_DG_NOTexact))


hold off
ylim([-0.015 0.08])
text(20,0,'IQ','units','data','HorizontalAlignment','center','VerticalAlignment','bottom','color','red')
% 	title('Peak Cross-BF Coincidence Detection - DENG AND GEISLER 1987')
% 	ylabel('SMP Slope (sp/sec/dB)')
xlabel('SNR (dB)','FontSize',TextFontSize)
set(gca,'XDir','reverse')
set(gca,'FontSize',AnnotFontSize)
xlim([-20 30])

Xcorner=0.13;
Xwidth=.35;
Xshift=0.15;

Ycorner=0.1;
Yshift=0.03;
Ywidth=.17;

TICKlength=0.03;

set(hSMP5,'Position',[Xcorner Ycorner Xwidth Ywidth],'TickLength',[TICKlength 0.025])
set(hSMP3,'Position',[Xcorner Ycorner+(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
set(hSMP1,'Position',[Xcorner Ycorner+2*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])

set(hSMP6,'Position',[Xcorner+(Xwidth+Xshift) Ycorner Xwidth Ywidth],'TickLength',[TICKlength 0.025])
set(hSMP4,'Position',[Xcorner+(Xwidth+Xshift) Ycorner+(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
set(hSMP2,'Position',[Xcorner+(Xwidth+Xshift) Ycorner+2*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])


orient tall
%% END SMP figure





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unlock all figures
figlist=get(0,'Children');
for FIGind=figlist
	set(FIGind,'Resize','on')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print all important figures
if PRINTyes
% 	figPRINTlist=[14 100 101 1000 1002 1 2000];
	figPRINTlist=[101 1000 1002 1 2000];
	for FIGind=figPRINTlist
		print(FIGind,'-PChesapeake')
% 		print(FIGind,'-PChoptank')
	end
% 	figure(1001)  % DFT figure won't print right for some reason???
% 	print -P'Adobe PDF'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(EHrlf_FIG)
% EHIQ_figname=['RLIQ_' num2str(TrackNum) '_0' num2str(UnitNum) '.eps'];
% eval(['print -deps ''' fullfile(ISHfigs_dir,EHIQ_figname) ''''])
% 
% figure(EHINrlf_FIG)
% EHIN_figname=['RLIN_' num2str(TrackNum) '_0' num2str(UnitNum) '.eps'];
% eval(['print -deps ''' fullfile(ISHfigs_dir,EHIN_figname) ''''])

% figure(PERhist_FIG)
% PERhist_figname=['PERhist_' num2str(TrackNum) '_0' num2str(UnitNum) '.eps'];
% eval(['print -deps ''' fullfile(ISHfigs_dir,PERhist_figname) ''''])
% 
% figure(SCC_fig)
% SCC_figname=['SCC_' num2str(TrackNum) '_0' num2str(UnitNum) '.eps'];
% eval(['print -deps ''' fullfile(ISHfigs_dir,SCC_figname) ''''])

figure(SMP_fig)
SMP_figname=['SMP_' num2str(TrackNum) '_0' num2str(UnitNum) '.eps'];
eval(['print -deps ''' fullfile(ISHfigs_dir,SMP_figname) ''''])


% Turn off saved PICS feature
SavedPICS=[]; SavedPICnums=[];
SavedPICSuse=0;

% % For ARO2005: ANmodel figure
% save SMP_FIG

return;

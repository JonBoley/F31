% FILE: fig_ThrQ10.m
%     : STIMind: 1-8, or 0 for all units
%
% For: RLFpaper
% From: ISH2003paper
% From: IHCON2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
global Recruit_dir RLFpaper_dir

STIMind=0;
global SEVimpExps MILDimpExps 

eval(['cd ''' fullfile(Recruit_dir,'POPData2') ''''])

MAXQUAL=1;
NUM_SRgroups=3;
SRbounds=[0.5 18];  % N-1 Boundaries

STIMlabel={'T2','T1','TB','FN','SPv','B2','B1','T05','SPf'};
FitLabels={'Good','No response','Masked','Plateau (<Rsat)','Non-monotonic','Decreasing','C2','Questionable'};
FitLabelsShort={'G','NR','M','P','NM','D','C2','Q'};

load POPnorm
POP_N=POP;
POPunitsN=POPunits;
levelsN=levels;
load POPimp
POP_I=POP;
POPunitsI=POPunits;
levelsI=levels;
clear POP POPunits levels

%%%%%%% FIG
set(0,'DefaultTextUnits','Normalized')

DataMarkerSize=4;

TextFontSize=9;
AnnotFontSize=9;
LabelFontSize=9;
SlopeFontSize=9;

TICKlength=0.025;

Datacolor='r';
Fitcolor='k';
SRcolors={'k','k','k'};
SRmarkers={'^','s','x'};

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'inches', [7.5 0.1 7 4]};  %J.NeuroPhys: 7in width
figure(10+STIMind); clf
set(gcf,figure_prop_name,figure_prop_val);


xmin=.1; xmax=20;
%%%%%%%%%%%%%%% Threshold/ Q10
load normema
%%%%%%% Add lines from Miller et al 1997 Normal data
QlowM97(:,1)=[.27 10]';  QhighM97(:,1)=QlowM97(:,1);
% 5th and 95th percentiles from Miller et al, 1997
QlowM97(:,2)=10^.2779*QlowM97(:,1).^.4708;
QhighM97(:,2)=10^.6474*QlowM97(:,1).^.4708;

for SRind=1:NUM_SRgroups
   if STIMind==0             % Take all units  (USE THIS ONE IN THE END?)
      Ndata=POPunitsN{SRind};
      Idata=POPunitsI{SRind};
   else                      % Take all units for a given stimulus
      Ndata=POP_N{STIMind,SRind};
      Idata=POP_I{STIMind,SRind};
   end
   
   QUALinds_N{SRind}=find(Ndata.QUALITY<=MAXQUAL);
   BF_N{SRind}=Ndata.BF(QUALinds_N{SRind});
   Thr_N{SRind}=Ndata.Thresh(QUALinds_N{SRind});
   Q10_N{SRind}=Ndata.Q10(QUALinds_N{SRind});
   
   QUALinds_Isev{SRind}=[];
	for EXPind=1:length(SEVimpExps)
		QUALinds_Isev{SRind}=[QUALinds_Isev{SRind} find((Idata.QUALITY<=MAXQUAL)&(strcmp(Idata.ExpDate,SEVimpExps{EXPind})))];
	end
	BF_Isev{SRind}=Idata.BF(QUALinds_Isev{SRind});
	Thr_Isev{SRind}=Idata.Thresh(QUALinds_Isev{SRind});
	Q10_Isev{SRind}=Idata.Q10(QUALinds_Isev{SRind});
   
   QUALinds_Imild{SRind}=[];
	for EXPind=1:length(MILDimpExps)
		QUALinds_Imild{SRind}=[QUALinds_Imild{SRind} find((Idata.QUALITY<=MAXQUAL)&(strcmp(Idata.ExpDate,MILDimpExps{EXPind})))];
	end
	BF_Imild{SRind}=Idata.BF(QUALinds_Imild{SRind});
	Thr_Imild{SRind}=Idata.Thresh(QUALinds_Imild{SRind});
	Q10_Imild{SRind}=Idata.Q10(QUALinds_Imild{SRind});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find % of Units in each SR group 
BFrange=[0.5,4];  % within impaired BF region [0.5,4] kHz

for SRind=1:3
   Nunits_N(SRind)=length(find((BF_N{SRind}>=BFrange(1))&(BF_N{SRind}<=BFrange(2))));
   Nunits_Imild(SRind)=length(find((BF_Imild{SRind}>=BFrange(1))&(BF_Imild{SRind}<=BFrange(2))));
   Nunits_Isev(SRind)=length(find((BF_Isev{SRind}>=BFrange(1))&(BF_Isev{SRind}<=BFrange(2))));
end
PERCunitsLMH_N54=Nunits_N/sum(Nunits_N);
PERCunitsLMH_N=[length(BF_N{1}) length(BF_N{2}) length(BF_N{3})];
PERCunitsLMH_N=PERCunitsLMH_N/sum(PERCunitsLMH_N);
PERCunitsLMH_Imild54=Nunits_Imild/sum(Nunits_Imild);
PERCunitsLMH_Imild=[length(BF_Imild{1}) length(BF_Imild{2}) length(BF_Imild{3})];
PERCunitsLMH_Imild=PERCunitsLMH_Imild/sum(PERCunitsLMH_Imild);
PERCunitsLMH_Isev54=Nunits_Isev/sum(Nunits_Isev);
PERCunitsLMH_Isev=[length(BF_Isev{1}) length(BF_Isev{2}) length(BF_Isev{3})];
PERCunitsLMH_Isev=PERCunitsLMH_Isev/sum(PERCunitsLMH_Isev);

disp(sprintf('Percent units in SR groups [L, M, H] for all BFs:\nNormal:\t\t%.2f, %.2f, %.2f\nMild:\t\t%.2f, %.2f, %.2f\nMod/Severe: %.2f, %.2f, %.2f\n', ...
   PERCunitsLMH_N,PERCunitsLMH_Imild,PERCunitsLMH_Isev))
disp(sprintf('\nPercent units in SR groups [L, M, H] for BFs in [%.1f,%.1f] kHz:\nNormal:\t\t%.2f, %.2f, %.2f\nMild:\t\t%.2f, %.2f, %.2f\nMod/Severe: %.2f, %.2f, %.2f\n', ...
   BFrange,PERCunitsLMH_N54,PERCunitsLMH_Imild54,PERCunitsLMH_Isev54))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%% Thresholds
ymin=-10; ymax=110;
h1=subplot(231);     % Normal
semilogx(normt(1,:),normt(2,:),'k')
hold on

for SRind=NUM_SRgroups:-1:1
	semilogx(BF_N{SRind},Thr_N{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
end
hleg1=legend('Miller et al. (1997)',sprintf('high (SR>%.1f)', ...
	SRbounds(2)),sprintf('med (%.1f<SR<=%.1f)', ...
	SRbounds(1),SRbounds(2)),sprintf('low (SR<=%.1f)',SRbounds(1)));
axis([xmin xmax ymin ymax])
title('Normal','FontSize',LabelFontSize)
ylabel('Threshold (dB SPL)','FontSize',LabelFontSize)
text(.08,.15,'NBTC','FontSize',TextFontSize,'VerticalAlignment','top','HorizontalAlignment','left');

h2=subplot(232);     % Mild
for SRind=NUM_SRgroups:-1:1
   semilogx(BF_Imild{SRind},Thr_Imild{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
   hold on
end
semilogx(normt(1,:),normt(2,:),'k')
axis([xmin xmax ymin ymax])
title('Mild Loss','FontSize',LabelFontSize)

h3=subplot(233);     % Severe
for SRind=NUM_SRgroups:-1:1
	semilogx(BF_Isev{SRind},Thr_Isev{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
semilogx(normt(1,:),normt(2,:),'k')
axis([xmin xmax ymin ymax])
title('Moderate/Severe Loss','FontSize',LabelFontSize)

%%%%%%%%%%% Q10
ymin=.3; ymax=20;
h4=subplot(234);     % Normal
for SRind=1:NUM_SRgroups
	loglog(BF_N{SRind},Q10_N{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_N{1}) length(BF_N{2}) length(BF_N{3})])),'Units','norm','FontSize',TextFontSize)
ylabel('Q_{10}','FontSize',LabelFontSize,'Interpreter','tex');

h5=subplot(235);     % Mild
for SRind=1:NUM_SRgroups
	loglog(BF_Imild{SRind},Q10_Imild{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_Imild{1}) length(BF_Imild{2}) length(BF_Imild{3})])),'Units','norm','FontSize',TextFontSize)
xlabel('Best Frequency (kHz)','FontSize',LabelFontSize)

h6=subplot(236);     % Severe
for SRind=1:NUM_SRgroups
	loglog(BF_Isev{SRind},Q10_Isev{SRind},'Color',SRcolors{SRind},'Marker',SRmarkers{SRind},'MarkerSize',DataMarkerSize,'LineStyle','none')
	hold on
end
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
axis([xmin xmax ymin ymax])
text(.7,.125,sprintf('N=%d',sum([length(BF_Isev{1}) length(BF_Isev{2}) length(BF_Isev{3})])),'Units','norm','FontSize',TextFontSize)



Xwidth=.31;
Xcorner=.062;
Xshift=0;
Ywidth=.42;
Ycorner=0.523;
Yshift=0.015;

set(h1,'Position',[Xcorner Ycorner Xwidth Ywidth],'XTickLabel',[],'YTick',[0:20:100], ...\
   'FontSize',AnnotFontSize,'TickLength',[TICKlength 0.025])
set(hleg1,'Position',[Xcorner+0.016    0.82    0.3500    0.0987],'FontSize',TextFontSize)
set(h2,'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth],'XTickLabel',[], ...
   'FontSize',AnnotFontSize,'YTick',[0:20:100],'YTickLabel',[],'TickLength',[TICKlength 0.025])
set(h3,'Position',[Xcorner+2*Xwidth+Xshift Ycorner Xwidth Ywidth],'XTickLabel',[], ...
   'FontSize',AnnotFontSize,'YTick',[0:20:100],'YTickLabel',[],'TickLength',[TICKlength 0.025])

set(h4,'Position',[Xcorner Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...\
   'YTickLabel',[1 10],'FontSize',AnnotFontSize,'TickLength',[TICKlength 0.025])
set(h5,'Position',[Xcorner+Xwidth+Xshift Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...
   'FontSize',AnnotFontSize,'YTickLabel',[],'TickLength',[TICKlength 0.025])
set(h6,'Position',[Xcorner+2*Xwidth+Xshift Ycorner-(Ywidth+Yshift) Xwidth Ywidth],'XTickLabel',[0.1 1 10], ...
   'FontSize',AnnotFontSize,'YTickLabel',[],'TickLength',[TICKlength 0.025])

drawnow
hold off

eval(['cd ''' fullfile(RLFpaper_dir,'Mfiles','EPSfigs') ''''])

%orient landscape

%print -djpeg100 -r150 ThrQ10

print -depsc fig_ThrQ10


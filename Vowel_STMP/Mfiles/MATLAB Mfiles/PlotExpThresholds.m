function PlotExpThresholds(ExpDirName)
% File: PlotExpThresholds.m
%
% Plots TC Thresholds for a given experiment

close all

set(0,'DefaultTextInterpreter','none');
% ExpDataDir = 'C:\NEL\ExpData';
% ExpDataDir = 'C:\Documents and Settings\Mike\My Documents\Work\Research\R03 Experiments\ExpData';
ExpDataDir = 'C:\Documents and Settings\Mike\My Documents\Work\Research\ISH2009\ExpData';
cd(ExpDataDir)

if ~exist(ExpDirName,'dir')
   dir
   error('No experiment by that name')
else
   cd(fullfile(ExpDataDir,ExpDirName))
   dir
   CALIBpics=findPics('calib');
   TrackUnitList=getTrackUnitList;
   NUMunits=size(TrackUnitList,1);
   Thrs_dbSPL = -999+zeros(1,NUMunits);
   BFs_kHz = -999+zeros(1,NUMunits);
   Q10s = -999+zeros(1,NUMunits);
   for ind = 1:NUMunits
      % Plot TC
      TCpics=findPics('tc',TrackUnitList(ind,:));
      if ~isempty(TCpics)
         Thrs_dbSPL_TEMP=[];
         BFs_kHz_TEMP=[];
         Q10s_TEMP=[];
         for TCind=1:length(TCpics)
            CALIB_PIC = CALIBpics(max(find((CALIBpics<TCpics(TCind)))));
				% [Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs(TCpics(TCind),CALIB_PIC,0);
				[Thrs_dbSPL_TEMP(TCind) BFs_kHz_TEMP(TCind) Q10s_TEMP(TCind)] = plotTCs_BFbyhand(TCpics(TCind),CALIB_PIC,0);
			end
         [Thrs_dbSPL(ind) index_TC] = min(Thrs_dbSPL_TEMP);
         BFs_kHz(ind) = BFs_kHz_TEMP(index_TC);
         Q10s(ind) = Q10s_TEMP(index_TC);
      end
   end
   
end

xmin=0.03; xmax=39; ymin=-20; ymax=110;
figure(3); clf
load normema

% %%%%%%% Add lines from Miller et al 1997 Normal data
% QlowM97(:,1)=[.27 10]';  QhighM97(:,1)=QlowM97(:,1);
% % 5th and 95th percentiles from Miller et al, 1997
% QlowM97(:,2)=10^.2779*QlowM97(:,1).^.4708;
% QhighM97(:,2)=10^.6474*QlowM97(:,1).^.4708;

%%%%%%% Add lines from Kale and Heinz (2009) normal-hearing CHIN data
QlowM97(:,1)=[.4 10]';  QhighM97(:,1)=QlowM97(:,1);
% 5th and 95th percentiles from Kale and Heinz (2009)
QlowM97(:,2)=10^.1572*QlowM97(:,1).^.1416;
QhighM97(:,2)=10^.7989*QlowM97(:,1).^.1416;

% FROM SUSHRUT
%  Equations for the Q10 lines are as follows:
% 
%         yhat=0.1416*xhat+0.4780; where xhat(i)=(BFhat(i)-mu(BFhat))/std(BFhat);
%                                        BFhat=log10(BF_normal_data); 
% 
% This equation fits the straight line to the normal data where xhat is the normalized BF vector i.e yhat is the predicted mean Q10 value for normal data.
% 
% Remaining equations are as follows:
% 
%       yhat=0.1416*xhat+0.7989; 95% conf. interval line
%       yhat=0.1416*xhat+0.1572;  5% conf. interval line
%       yhat=0.1416*xhat+0.9013; 99% conf. interval line
%       yhat=0.1416*xhat+0.0548;  1% conf. interval line


subplot(211)
semilogx(BFs_kHz,Thrs_dbSPL,'x')
hold on
semilogx(normt(1,:),normt(2,:),'k')
hold off
ylabel('dB SPL'); xlabel('Frequency (kHz)');
axis([xmin xmax ymin ymax]);
set(gca,'YTick',[0:10:100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
title(sprintf('Exp: %s',ExpDirName))
grid on

subplot(212)
loglog(BFs_kHz,Q10s,'x')
hold on
loglog(QlowM97(:,1),QlowM97(:,2),'k-')
loglog(QhighM97(:,1),QhighM97(:,2),'k-')
hold off
ylabel('Q10'); xlabel('Frequency (kHz)');
xlim([xmin xmax]);
ylim([.1 10]);
% set(gca,'YTick',[0 20 40 60 80 100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
grid on

orient tall

plotTCs(findPics('tc'),CALIB_PIC)
title(ExpDirName)
orient landscape

return;



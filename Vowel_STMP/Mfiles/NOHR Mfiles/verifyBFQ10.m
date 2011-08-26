function unit=verifyBFQ10(unit)
% FILE: verifyBFQ10.m
% Modified 11Feb2005 M. Heinz for NOHR Data
%
% Created 7/31/02: for choosing BF and verofying Q10
%
% 1) Picks BF & Threshold from actual data points
% 2) Generates a smoothed TC (without bias at BF) and saves as unit.tcdata(:,3)
% 3) Finds Q10 based on smoothed TC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set(0,'DefaultTextInterpreter','none');
set(0,'DefaultTextUnits','data')

%% Show related BF values, and choose BF by hand
% verifyBFQ10 will show minTC, book, TChand (if exists, DEFAULT), and unit.BF (if exists, i.e., re-check old data)

h14=figure(14); clf;  
set(h14,'Position',[107   340   914   358])

TFiltWidthTC=3;

TextFontSize=8;
DataMarkerSize=12;
DataMarkStyle='b.';
DataFitStyle='b-';

xmin=0.03; xmax=39; ymin=-20; ymax=110;
load normema

%%%%%% PLOT TUNING CURVE
h_line1 = semilogx(unit.TC.tcdata(:,1),unit.TC.tcdata(:,2),DataMarkStyle,'MarkerSize',DataMarkerSize);
hold on
h_line2 = semilogx(unit.TC.tcdata(:,1),trifilt(unit.TC.tcdata(:,2)',TFiltWidthTC)',DataFitStyle);
semilogx(normt(1,:),normt(2,:),'k')
ylabel('dB SPL'); xlabel('Frequency (kHz)');
axis([xmin xmax ymin ymax]);
set(gca,'YTick',[0 20 40 60 80 100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
if(isfield(unit,'calNUM'))
   title(sprintf('%s (Cal: P%d)',unit.titletext,unit.calNUM))
else
   title(sprintf('%s (Cal: P%d)',unit.titletext,unit.Calib.calPICT))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Allow user to verify BF 
h_minTC = semilogx(unit.TC.minTC_BF_kHz*ones(1,2),unit.TC.minTC_Thresh_dBSPL+[5 30],'m:');
text(.1,.2,'minTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','m')
if isfield(unit.TC,'handTC_BF_kHz')
   h_handTC = semilogx(unit.TC.handTC_BF_kHz*ones(1,2),unit.TC.handTC_Thresh_dBSPL+[-5 -30],'c:');
   text(.1,.1,'handTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','c')
end
if isfield(unit,'BF_kHz')
   %    unit.TC.oldTC_BF_kHz=unit.BF_kHz;
   %    unit.TC.oldTC_Thresh_dBSPL=unit.Thresh_dBSPL;
   h_unitTC = semilogx(unit.BF_kHz*ones(1,2),unit.Thresh_dBSPL+[-15 15],'k:');
   text(.1,.15,'unitTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','k')
end

if isfield(unit,'BF_kHz')
   pickThr = unit.Thresh_dBSPL;
   pickBF = unit.BF_kHz;
elseif isfield(unit.TC,'handTC_BF_kHz')
   pickThr = unit.TC.handTC_Thresh_dBSPL;
   pickBF = unit.TC.handTC_BF_kHz;
else
   pickThr = unit.TC.minTC_Thresh_dBSPL;
   pickBF = unit.TC.minTC_BF_kHz;
end

text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
h_textBF= text(.05,.4,text_str,'Units','norm');

h_pickBF=text(pickBF,pickThr,'\uparrow','Interpreter','tex','FontSize',16, ...
   'VerticalAlignment','top','HorizontalAlignment','center');

if ~isempty(unit.Calib.cal_freq{3})
   for i=1:length(unit.Calib.cal_freq{3})
      semilogx(unit.Calib.cal_freq{3}(i)/1000*ones(1,2),pickThr+[-7 7],'r-');
      text(.1*2^(i-1),.25,'TB','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','r')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop to wait for verifying BF 
x15=get(h14,'Position');  % These commands insure mouse is over Figure so that cursor keys work
set(0,'PointerLocation',x15(1:2)+[10 x15(4)/2])
set(gcf,'CurrentCharacter','x')
loc = find(unit.TC.tcdata(:,1)==pickBF);
if isempty(loc)
   [yy,loc]=min(abs(unit.TC.tcdata(:,1)-pickBF));
   disp('unit.BF NOT chosen from data points')
   pickThr = unit.TC.tcdata(loc,2);
   pickBF = unit.TC.tcdata(loc,1);
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
         loc=min(length(unit.TC.tcdata(:,1)),loc+1);
      case 29  %%% 'RIGHT cursor'
         loc=max(1,loc-1);
      end
   end
	
   pickThr = unit.TC.tcdata(loc,2);
   pickBF = unit.TC.tcdata(loc,1);
   set(h_pickBF,'Position',[pickBF pickThr]);
   text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
   set(h_textBF,'String',text_str);
end  % End wait for verifying BF

% Set unit BF/Threshold to picked BF/Threshold
unit.BF_kHz=pickBF;
unit.Thresh_dBSPL=pickThr;

%% Generate smoothed TC, but avoiding upward bias at BF (tip)
% Fits each side separately, and then sets equal to actual data point at BF
% i.e., smoothes sides (e.g., 10 dB up) without biasing threshold at BF upward
unit.TC.tcdata(1:loc,3)=trifilt(unit.TC.tcdata(1:loc,2)',TFiltWidthTC)';
unit.TC.tcdata(loc:end,3)=trifilt(unit.TC.tcdata(loc:end,2)',TFiltWidthTC)';
unit.TC.tcdata(loc,3)=unit.TC.tcdata(loc,2);
set(h_line2,'YData',unit.TC.tcdata(:,3));

% pass smoothed tcdata for q10 calculation (based on actual data point at BF, and smoothed TC otherwise 
% This avoids the bias in smoothing at the tip, i.e., raising threshold at BF
[unit.Q10,unit.TC.Q10fhi,unit.TC.Q10flo,unit.TC.Q10lev]= ...
   findQ10(unit.TC.tcdata(:,1),unit.TC.tcdata(:,3),unit.BF_kHz);

semilogx([unit.TC.Q10flo unit.TC.Q10fhi],unit.TC.Q10lev*ones(1,2),'k-');
text_str = sprintf('%s %6.3f %s\n%s %4.2f %s\n%s %4.1f','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL','Q10: ',unit.Q10);
set(h_textBF,'String',text_str);

temp=input('Verify Q10 is OK [Return to accept; anything else to mark as QUESTIONABLE]: ');
if ~isempty(temp)
   unit.Q10bad=1;
   unit.Q10=NaN;
end


hold off
return;


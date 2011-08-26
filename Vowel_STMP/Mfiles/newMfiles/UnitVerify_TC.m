function UnitVerify_TC(ExpDate,UnitName,justPLOT)
% FILE: UnitVerify_TC.m
% Modified 11Feb2005 M. Heinz for NOHR Data
%
% Created 7/31/02: for choosing BF and verofying Q10
%
% 1) Picks BF & Threshold from actual data points
% 2) Generates a smoothed TC (without bias at BF) and saves as unit.tcdata(:,3)
% 3) Finds Q10 based on smoothed TC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global ROOT_dir ExpList


%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ExpDate','var')
   ExpDate=0;
   %%% HARD CODE FOR NOW
   ExpDate='111804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
   %%% HARD CODE FOR NOW
   UnitName='1.28'
   while ~ischar(UnitName)
      UnitName=input('Enter Unit Name (e.g., ''3.07''): ');
   end
end
if ~exist('justPLOT','var')
   justPLOT=0;
end

% Find the full Experiment Name 
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(ExpList)
   if ~isempty(strfind(ExpList{i},ExpDateText))
      ExpName=ExpList{i};
      break;
   end
end
if ~exist('ExpName','var')
   disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
   disp(strvcat(ExpList))
   beep
   error('STOPPED');
end

% Parse out the Track and Unit Number 
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(ROOT_dir,'ExpData');
anal_dir=fullfile(ROOT_dir,'Data Analysis');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

disp(sprintf('Processing (verify_TC) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

% Verify Unit Number is legitimate
UnitPicList=findPics('*',[TrackNum,UnitNum]);
if isempty(UnitPicList)
   error(sprintf('Unit Number: %d.%02d does not exist!!!!!',TrackNum,UnitNum))
end

%%%%%%%% Create unit structure for this unit
% LOAD UNITSdata file, if it exists, otherwise take from DataList
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
%%%%%%%%%%%%%%%%
%% load unit data if it exists,
if ~isempty(dir(fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum))))
   eval(['load ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''''])
   % Make sure the new TC is here, rather than the old (atten, not dB SPL)
   if ~isfield(unit.TC{1},'tcdata')
      disp(['   *** Loading file: "' FileName '" because UNITSdata/unit has OLD TCdata'])
      eval(['load ' FileName])
      unit.TC=DataList.Units{TrackNum,UnitNum}.TC;
      unit.Info.Threshold_dBSPL=DataList.Units{TrackNum,UnitNum}.Info.Threshold_dBSPL;
      unit.Info=rmfield(unit.Info,'Threshold_dBatten');
      unit.Info.Q10=[];
   end
else
   disp(['   *** Loading file: "' FileName '" because UNITSdata/unit does not exist yet'])
   eval(['load ' FileName])
   unit=DataList.Units{TrackNum,UnitNum};
   unit.IgnorePicNums=intersect(UnitPicList,DataList.IgnorePicNums);
end


% See if Q10 already calculated
CALCtc=1;
if ~isempty(unit.Info.Q10)
   if ~justPLOT
      beep
      temp=input(sprintf('--- Q10 field has already been calculated for this unit, do you want to recalculate (''y''/[''n''])? '));
      if isempty(temp)
         CALCtc=0;
      elseif ~strcmp(upper(temp(1)),'Y')
         CALCtc=0;
      end
   else
      CALCtc=0;
   end
end

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
h_line1 = semilogx(unit.TC{unit.Info.TCindToUse}.tcdata(:,1),unit.TC{unit.Info.TCindToUse}.tcdata(:,2),DataMarkStyle,'MarkerSize',DataMarkerSize);
hold on
h_line2 = semilogx(unit.TC{unit.Info.TCindToUse}.tcdata(:,1),trifilt(unit.TC{unit.Info.TCindToUse}.tcdata(:,2)',TFiltWidthTC)',DataFitStyle);
semilogx(normt(1,:),normt(2,:),'k')
ylabel('dB SPL'); xlabel('Frequency (kHz)');
title(sprintf('Exp%s, Unit: %s',ExpDate,UnitName))
axis([xmin xmax ymin ymax]);
set(gca,'YTick',[0 20 40 60 80 100])
set(gca,'XTick',[.1 1 10],'XTickLabel',[.1 1 10])
% if(isfield(unit,'calNUM'))
%    title(sprintf('%s (Cal: P%d)',unit.titletext,unit.calNUM))
% else
%    title(sprintf('%s (Cal: P%d)',unit.titletext,unit.Calib.calPICT))
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Allow user to verify BF 
h_minTC = semilogx(unit.TC{unit.Info.TCindToUse}.minTC_BF_kHz*ones(1,2),unit.TC{unit.Info.TCindToUse}.minTC_Thresh_dBSPL+[5 30],'m:');
text(.1,.2,'minTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','m')
if isfield(unit.TC{unit.Info.TCindToUse},'book_BF_kHz')
   h_handTC = semilogx(unit.TC{unit.Info.TCindToUse}.book_BF_kHz*ones(1,2),unit.TC{unit.Info.TCindToUse}.book_Thr_dBSPL+[-5 -30],'c:');
   text(.1,.1,'bookTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','c')
end
% if isfield(unit,'BF_kHz')
%    unit.TC{unit.Info.TCindToUse}.oldTC_BF_kHz=unit.BF_kHz;
%    unit.TC{unit.Info.TCindToUse}.oldTC_Thresh_dBSPL=unit.Thresh_dBSPL;
h_unitTC = semilogx(unit.Info.BF_kHz*ones(1,2),unit.Info.Threshold_dBSPL+[-15 15],'k--');
text(.1,.15,'unitTC','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','k')
% end

% if isfield(unit,'BF_kHz')
pickThr = unit.Info.Threshold_dBSPL;
pickBF = unit.Info.BF_kHz;
% elseif isfield(unit.TC{unit.Info.TCindToUse},'handTC_BF_kHz')
%    pickThr = unit.TC{unit.Info.TCindToUse}.handTC_Thresh_dBSPL;
%    pickBF = unit.TC{unit.Info.TCindToUse}.handTC_BF_kHz;
% else
%    pickThr = unit.TC{unit.Info.TCindToUse}.minTC_Thresh_dBSPL;
%    pickBF = unit.TC{unit.Info.TCindToUse}.minTC_BF_kHz;
% end

text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
h_textBF= text(.05,.98,text_str,'Units','norm','FontSize',TextFontSize,'VerticalAlignment','top');

h_pickBF=text(pickBF,pickThr,'\uparrow','Interpreter','tex','FontSize',16, ...
   'VerticalAlignment','top','HorizontalAlignment','center');

% if ~isempty(unit.Calib.cal_freq{3})
%    for i=1:length(unit.Calib.cal_freq{3})
%       semilogx(unit.Calib.cal_freq{3}(i)/1000*ones(1,2),pickThr+[-7 7],'r-');
%       text(.1*2^(i-1),.25,'TB','FontSize',TextFontSize,'Units','Norm','HorizontalAlignment','center','Color','r')
%    end
% end

if CALCtc
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Loop to wait for verifying BF 
   x15=get(h14,'Position');  % These commands insure mouse is over Figure so that cursor keys work
   set(0,'PointerLocation',x15(1:2)+[10 x15(4)/2])
   set(gcf,'CurrentCharacter','x')
   loc = find(unit.TC{unit.Info.TCindToUse}.tcdata(:,1)==pickBF);
   if isempty(loc)
      [yy,loc]=min(abs(unit.TC{unit.Info.TCindToUse}.tcdata(:,1)-pickBF));
      disp('unit.BF NOT chosen from data points')
      pickThr = unit.TC{unit.Info.TCindToUse}.tcdata(loc,2);
      pickBF = unit.TC{unit.Info.TCindToUse}.tcdata(loc,1);
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
            loc=min(length(unit.TC{unit.Info.TCindToUse}.tcdata(:,1)),loc+1);
         case 29  %%% 'RIGHT cursor'
            loc=max(1,loc-1);
         end
      end
      
      pickThr = unit.TC{unit.Info.TCindToUse}.tcdata(loc,2);
      pickBF = unit.TC{unit.Info.TCindToUse}.tcdata(loc,1);
      set(h_pickBF,'Position',[pickBF pickThr]);
      text_str = sprintf('%s %6.3f %s\n%s %4.2f %s','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL');
      set(h_textBF,'String',text_str);
   end  % End wait for verifying BF
   
   % Set unit BF/Threshold to picked BF/Threshold
   unit.Info.BF_kHz=pickBF;
   unit.Info.Threshold_dBSPL=pickThr;
   
   %% Generate smoothed TC, but avoiding upward bias at BF (tip)
   % Fits each side separately, and then sets equal to actual data point at BF
   % i.e., smoothes sides (e.g., 10 dB up) without biasing threshold at BF upward
   unit.TC{unit.Info.TCindToUse}.tcdata(1:loc,3)=trifilt(unit.TC{unit.Info.TCindToUse}.tcdata(1:loc,2)',TFiltWidthTC)';
   unit.TC{unit.Info.TCindToUse}.tcdata(loc:end,3)=trifilt(unit.TC{unit.Info.TCindToUse}.tcdata(loc:end,2)',TFiltWidthTC)';
   unit.TC{unit.Info.TCindToUse}.tcdata(loc,3)=unit.TC{unit.Info.TCindToUse}.tcdata(loc,2);
end
set(h_line2,'YData',unit.TC{unit.Info.TCindToUse}.tcdata(:,3));
   
% pass smoothed tcdata for q10 calculation (based on actual data point at BF, and smoothed TC otherwise 
% This avoids the bias in smoothing at the tip, i.e., raising threshold at BF
[unit.TC{unit.Info.TCindToUse}.Q10,unit.TC{unit.Info.TCindToUse}.Q10fhi,unit.TC{unit.Info.TCindToUse}.Q10flo,unit.TC{unit.Info.TCindToUse}.Q10lev]= ...
   findQ10(unit.TC{unit.Info.TCindToUse}.tcdata(:,1),unit.TC{unit.Info.TCindToUse}.tcdata(:,3),unit.Info.BF_kHz);

unit.Info.Q10=unit.TC{unit.Info.TCindToUse}.Q10;

semilogx([unit.TC{unit.Info.TCindToUse}.Q10flo unit.TC{unit.Info.TCindToUse}.Q10fhi],unit.TC{unit.Info.TCindToUse}.Q10lev*ones(1,2),'k-');
text_str = sprintf('%s %6.3f %s\n%s %4.2f %s\n%s %4.1f','BF:',pickBF,'kHz.','Thresh:',pickThr,'dB SPL','Q10: ',unit.Info.Q10);
set(h_textBF,'String',text_str);
   
if CALCtc
   temp=input('Verify Q10 is OK [Return to accept; anything else to mark as QUESTIONABLE]: ');
   if ~isempty(temp)
      unit.Info.Q10bad=1;
      unit.Info.Q10=NaN;
   end
   
   %%%%%%%%%%%% Save unit data only if updated
   disp(['   *** Saving unit file: "' fullfile(ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) '" ...']) 
   %    beep
   %    disp(' SAVING TURNED OFF')
   eval(['save ''' fullfile(data_dir,ExpName,'UNITSdata',sprintf('unit.%d.%02d.mat',TrackNum,UnitNum)) ''' unit'])
else
   disp('TC left as is, nothing saved')
end

hold off

return;


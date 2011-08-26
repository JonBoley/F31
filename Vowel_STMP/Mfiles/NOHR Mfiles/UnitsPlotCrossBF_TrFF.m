function UnitsPlotCrossBF_TrFF(ExpDate,Unit1Name,Unit2Name)
% File: UnitsPlotCrossBF_TrFF.m
% Date: 25Sep2004 (M. Heinz) (modified from UnitPlot_TrBF and UnitPLot_TrFF)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitNames: '3.07' (converted later)
%
% Plots Tone_reFF data for 2 units from a given experiment PLUS the same data from SCC Analysis of Cross-Unit Correlation.  
% Loads 'UNITSdata/unit.T.U.mat' file.
% UnitAnal_TrFF.m performs the relevant analysis.
%

global NOHR_dir NOHR_ExpList

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='090204'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''090204''): ');
   end
end
if ~exist('Unit1Name','var')
   Unit1Name=0;
   %%% HARD CODE FOR NOW
Unit1Name='1.07'
   while ~ischar(Unit1Name)
      Unit1Name=input('Enter Unit 1 Name (e.g., ''3.07''): ');
   end
end
if ~exist('Unit2Name','var')
   Unit2Name=0;
   %%% HARD CODE FOR NOW
Unit2Name='1.09'
   while ~ischar(Unit2Name)
      Unit2Name=input('Enter Unit 2 Name (e.g., ''3.07''): ');
   end
end
UnitNames={Unit1Name,Unit2Name};

%%%% Find the full Experiment Name 
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(NOHR_ExpList)
   if ~isempty(strfind(NOHR_ExpList{i},ExpDateText))
      ExpName=NOHR_ExpList{i};
      break;
   end
end
if ~exist('ExpName','var')
   disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
   disp(strvcat(NOHR_ExpList))
   beep
   break
end

%%%% Parse out the Track and Unit Numbers
TrackNums=cell(1,2);
UnitNums=cell(1,2);
for i=1:2
   TrackNums{i}=str2num(UnitNames{i}(1:strfind(UnitNames{i},'.')-1));
   UnitNums{i}=str2num(UnitNames{i}(strfind(UnitNames{i},'.')+1:end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
disp(sprintf('Plotting Cross_BF Analysis (TreFF) Experiment: ''%s''; Units: %d.%02d and %d.%02d',ExpName,TrackNums{1},UnitNums{1},TrackNums{2},UnitNums{2}))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
units=cell(1,2);
UnitFileNames=cell(1,2);
%%%% Load unit structures for these units
for i=1:2
   UnitFileNames{i}=sprintf('unit.%d.%02d.mat',TrackNums{i},UnitNums{i});
   eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileNames{i}) ''');'])
   % If UNITSdata file does not exist, run analysis
   if isempty(ddd)
      UnitAnal_TrFF(ExpDate,UnitNames{i},0);
   end
   eval(['load ''' fullfile('UNITSdata',UnitFileNames{i}) ''''])
   
   % Make sure there is Tone_reFF data
   if ~isfield(unit,'Tone_reFF')
      disp(sprintf('*** No Tone_reFF data for unit %s!!!',UnitNames{i}))
      beep
      return;
   end
   
   % If Tone_reFF analysis is not completed, run here
   if ~isfield(unit.Tone_reFF,'synch')
      UnitAnal_TrFF(ExpDate,UnitNames{i},0);
      eval(['load ''' fullfile('UNITSdata',UnitFileNames{i}) ''''])
   end
   units{i}=unit;
   clear unit   
end

%%%% Plot data
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FIG.FontSize=8;

levels_dBSPL=union(units{1}.Tone_reFF.levels_dBSPL,units{2}.Tone_reFF.levels_dBSPL);

%%%% Store data from 2 units and SCC analysis
BFs_kHz=NaN+ones(1,5);
rate=NaN+zeros(length(levels_dBSPL),5);  % NSCCrates = 3: NSCC_0, NSCC_HLCD, NSCC_max
synch=NaN+zeros(length(levels_dBSPL),3);
phase=NaN+zeros(length(levels_dBSPL),3);
rateSlopes=NaN+zeros(length(levels_dBSPL),5);
synchSlopes=NaN+zeros(length(levels_dBSPL),3);
phaseSlopes=NaN+zeros(length(levels_dBSPL),3);

% Store rate, synch, and phase data from 2 units
for i=1:2
   BFs_kHz(i)=units{i}.Info.BF_kHz;
   if i==1
      FF_kHz=units{i}.Tone_reFF.freqs_kHz;
   else
      if units{i}.Tone_reFF.freqs_kHz~=FF_kHz
         error('Difference in Fixed Frequency between 2 units')
      end
   end
   %%%% Take out (mark with NaNs) insignificant phase and synch values
   warning off % don't need to see "divide by zero" error over and over
   units{i}.Tone_reFF.synch=units{i}.Tone_reFF.synch.*units{i}.Tone_reFF.RaySig./units{i}.Tone_reFF.RaySig;
   units{i}.Tone_reFF.phase=units{i}.Tone_reFF.phase.*units{i}.Tone_reFF.RaySig./units{i}.Tone_reFF.RaySig;
   warning on
   
   for LEVind=1:length(units{i}.Tone_reFF.levels_dBSPL)
      rate(find(levels_dBSPL==units{i}.Tone_reFF.levels_dBSPL(LEVind)),i)=units{i}.Tone_reFF.rate(LEVind);
      synch(find(levels_dBSPL==units{i}.Tone_reFF.levels_dBSPL(LEVind)),i)=units{i}.Tone_reFF.synch(LEVind);
      phase(find(levels_dBSPL==units{i}.Tone_reFF.levels_dBSPL(LEVind)),i)=units{i}.Tone_reFF.phase(LEVind);
   end
end
BFs_kHz(3)=geomean(BFs_kHz(1:2))/1.05;   % NSCC_0
BFs_kHz(4)=geomean(BFs_kHz(1:2));        % NSCC_HLCD
BFs_kHz(5)=geomean(BFs_kHz(1:2))*1.05;   % NSCC_max

%%%%%%%%%%%%% 9/25/04  TO DO   %%%%%%%%%%%%%%%%%%%%
% *1) Basics for 2 units are steup for TrFF
% *2) Found 2 good units
% *3) NOW, need to do SCC analysis on these two units, and plot that data HERE!!!!
% *4) Do same for EHrFF
% *5) ONTO reBF!!!!!!!! & simulate 2 neurons

levels_dBSPL_SCC=sort(intersect(units{1}.Tone_reFF.levels_dBSPL,units{2}.Tone_reFF.levels_dBSPL));

calcSAC=0;
HLCD_us=0;

ADHOC_RATE=50;
for LEVind=length(levels_dBSPL_SCC):-1:1  % go from high to low, to get HLCD first
   % for LEVind=2:2
   for i=1:2
      LEVinds(i)=find(units{i}.Tone_reFF.levels_dBSPL==levels_dBSPL(LEVind));
      SCCtitles{i}=sprintf('UNIT %d: %s  [BF=%.4f kHz, Thr= %.f dB atten]', ...
         i,units{i}.Info.Unit,units{i}.Info.BF_kHz,units{i}.Info.Threshold_dBatten);
   end
   SCCtitles{3}=sprintf('Exp: %s;     [SCCdemo_TONEreFF]  (FF=%.2f kHz)   (%.f dB SPL)',units{1}.Info.ExpName,FF_kHz,levels_dBSPL_SCC(LEVind));
   disp(sprintf('%s\nProcessing: %.f dB SPL\n%s',repmat('*',1,21),levels_dBSPL_SCC(LEVind),repmat('*',1,21)))
   
   [CD_us,NSCC_0,NSCC_CD,NSCC_HLCD,NSCC_max]= ...
      SCCdemo_TrFF(units{1}.Tone_reFF.picNums{LEVinds(1)},units{2}.Tone_reFF.picNums{LEVinds(2)}, ...
      units{1}.Tone_reFF.excludeLines{LEVinds(1)},units{2}.Tone_reFF.excludeLines{LEVinds(2)},HLCD_us,SCCtitles,calcSAC);
   
   if LEVind==length(levels_dBSPL_SCC)
      HLCD_us=CD_us;  % Set HLCD to be used for rest
      NSCC_HLCD=NSCC_CD;
   end
   
   rate(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),3)= NSCC_0 * ADHOC_RATE;
   rate(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),4)= NSCC_HLCD * ADHOC_RATE;
   rate(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),5)= NSCC_max * ADHOC_RATE;
   
   %    input('Enter')
end
   


figure(12); clf
set(gcf,'pos',[420   199   976   761])
%%%% Rate
subplot(321)
clear LEGtext
for LEVind=1:length(levels_dBSPL)
   semilogx(BFs_kHz,rate(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
   hold on
   LEGtext{LEVind}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
end
xLIMS=[10^floor(log10(FF_kHz*.5)) 10^ceil(log10(FF_kHz*2))];
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
ylim([0 300])  % Fixed ordinate for all plots
xlim(xLIMS)
% set(gca,'XTickLabel',[.1 1 10])
ylabel('Rate (sp/sec)')
xlabel('[3 NSCC results: (L) 0-delay, (M) HL-CD, (R) max NSCC]')
ht1=title(sprintf('Exp: %s;     [UnitsPlotCrossBF_TONEreFF]  (FF=%.2f kHz)\nUNIT 1: %s  [BF=%.4f kHz, Thr= %.f dB atten] %40s UNIT 2: %s  [BF=%.4f kHz, Thr= %.f dB atten]', ...
   units{1}.Info.ExpName,FF_kHz,units{1}.Info.Unit,units{1}.Info.BF_kHz,units{1}.Info.Threshold_dBatten,' ',units{2}.Info.Unit,units{2}.Info.BF_kHz,units{2}.Info.Threshold_dBatten));
set(ht1,'units','norm','Interpreter','none','pos',[1.2 1.0292 0])
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Rate Change
subplot(322)
clear LEGtextSlope
for LEVind=2:length(levels_dBSPL)
   rateSlopes(LEVind,:)=(rate(LEVind,:)-rate(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
   semilogx(BFs_kHz,rateSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
   hold on
   LEGtextSlope{LEVind}=sprintf('%.f-%.f dB SPL',levels_dBSPL(LEVind-1),levels_dBSPL(LEVind));
end
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
semilogx(xLIMS,[0 0],'k-')
ylabel('Rate Slope (sp/sec/dB)')
ymax_DeltaRate=10;  % 0 to 300 in 30 dB
ylim([-1 1]*ymax_DeltaRate)
xlim(xLIMS)
% set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Synch
subplot(323)
for LEVind=1:length(levels_dBSPL)
   semilogx(BFs_kHz(1:3),synch(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
   hold on
end
hleg=legend(LEGtext,4);
set(hleg,'FontSize',6)
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
ylabel('Synch')
ylim([0 1])
xlim(xLIMS)
% set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Synch Change
subplot(324)
for LEVind=2:length(levels_dBSPL)
   synchSlopes(LEVind,:)=(synch(LEVind,:)-synch(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
   semilogx(BFs_kHz(1:3),synchSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
   hold on
end
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
semilogx(xLIMS,[0 0],'k-')
hlegSlope=legend(LEGtextSlope{2:length(levels_dBSPL)},4);
set(hlegSlope,'FontSize',6)
ylabel('Synch Slope (1/dB)')
ymax_DeltaSynch=1/30;  % 0 to 1 in 30 dB
ylim([-1 1]*ymax_DeltaSynch)
xlim(xLIMS)
% set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Phase
subplot(325)
for LEVind=1:length(levels_dBSPL)
   semilogx(BFs_kHz(1:3),phase(LEVind,:),'Marker',LEVELmarkers{LEVind}, ...
      'Color',LEVELcolors{LEVind},'LineStyle','none')
   hold on
end
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
ylabel('Phase (rad)')
xlabel('Frequency (kHz)')
ylim([-pi pi])
xlim(xLIMS)
% set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Phase Change
subplot(326)
for LEVind=2:length(levels_dBSPL)
   phaseLL=repmat(phase(LEVind-1,:),3,1);
   phaseLL(2,:)=phaseLL(1,:)+2*pi;
   phaseLL(3,:)=phaseLL(3,:)-2*pi;
   phaseHL=repmat(phase(LEVind,:),3,1);
   phaseDiffs=phaseHL-phaseLL;
   [phaseMinDiffs,phaseMinInds]=min(abs(phaseDiffs));
   for FREQind=1:size(phase,2)
      phaseSlopes(LEVind,FREQind)=phaseDiffs(phaseMinInds(FREQind),FREQind)/diff(levels_dBSPL(LEVind-1:LEVind));
   end
   semilogx(BFs_kHz(1:3),phaseSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
   hold on
end
semilogx(FF_kHz*ones(1,2),[-1e6 1e6],'k:')
semilogx(xLIMS,[0 0],'k-')
ylabel('Phase Slope (rad/dB)')
xlabel('Frequency (kHz)')
% ylim([-1 1]*max(ceil(max(abs(phaseSlopes))/.02)*.02))  
ymax_DeltaSynch=0.05;  % ~pi/2 in 30 dB
ylim([-1 1]*ymax_DeltaSynch)
xlim(xLIMS)
% set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)




return;

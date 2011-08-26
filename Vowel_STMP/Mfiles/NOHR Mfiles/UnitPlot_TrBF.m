function UnitPlot_TrBF(ExpDate,UnitName)
% File: UnitPlot_TrBF.m
% Date: 13Sep2004 (M. Heinz)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Plots Tone_reBF data for a given experiment and unit.  Loads 'UNITSdata/unit.T.U.mat' file.
% UnitAnal_TrBF.m performs the relevant analysis.
%

global NOHR_dir NOHR_ExpList

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='070804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
%%% HARD CODE FOR NOW
UnitName='3.07'
   while ~ischar(UnitName)
      UnitName=input('Enter Unit Name (e.g., ''3.07''): ');
   end
end

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

%%%% Parse out the Track and Unit Number 
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
disp(sprintf('Plotting  (TreBF) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
% If UNITSdata file does not exist, run analysis
if isempty(ddd)
   UnitAnal_TrBF(ExpDate,UnitName,0);
end
eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

% Make sure there is Tone_reBF data
if ~isfield(unit,'Tone_reBF')
   disp(sprintf('*** No Tone_reBF data for this unit!!!'))
   beep
   return;
end

% If Tone_reBF analysis is not completed, run here
if ~isfield(unit.Tone_reBF,'synch')
   UnitAnal_TrBF(ExpDate,UnitName,0);
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end
   
%%%% Plot data
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FIG.FontSize=8;

rateSlopes=NaN+zeros(size(unit.Tone_reBF.rate));
synchSlopes=NaN+zeros(size(unit.Tone_reBF.rate));
phaseSlopes=NaN+zeros(size(unit.Tone_reBF.rate));


%%%% Take out (mark with NaNs) insignificant phase and synch values
warning off % don't need to see "divide by zero" error over and over
unit.Tone_reBF.synch=unit.Tone_reBF.synch.*unit.Tone_reBF.RaySig./unit.Tone_reBF.RaySig;
unit.Tone_reBF.phase=unit.Tone_reBF.phase.*unit.Tone_reBF.RaySig./unit.Tone_reBF.RaySig;
warning on


figure(12); clf
set(gcf,'pos',[420   199   976   761])
%%%% Rate
subplot(321)
clear LEGtext
for LEVind=1:length(unit.Tone_reBF.levels_dBSPL)
   semilogx(unit.Tone_reBF.freqs_kHz,unit.Tone_reBF.rate(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle',LEVELlines{LEVind})
   hold on
   LEGtext{LEVind}=sprintf('%.f dB SPL',unit.Tone_reBF.levels_dBSPL(LEVind));
end
semilogx(unit.Info.BF_kHz*ones(1,2),[0 1e6],'k:')
% ylim([0 ceil(max(max(unit.Tone_reBF.rate))/50)*50])
ylim([0 300])  % Fixed ordinate for all plots
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
ylabel('Rate (sp/sec)')
ht1=title(sprintf('Exp: %s; Unit: %s  [BF=%.4f kHz, Thr= %.f dB atten]\n[UnitPlot_TONEreBF]',unit.Info.ExpName,unit.Info.Unit,unit.Info.BF_kHz,unit.Info.Threshold_dBatten));
set(ht1,'units','norm','Interpreter','none','pos',[1.2 1.0292 0])
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Rate Change
subplot(322)
clear LEGtextSlope
for LEVind=2:length(unit.Tone_reBF.levels_dBSPL)
   rateSlopes(LEVind,:)=(unit.Tone_reBF.rate(LEVind,:)-unit.Tone_reBF.rate(LEVind-1,:))/diff(unit.Tone_reBF.levels_dBSPL(LEVind-1:LEVind));
   semilogx(unit.Tone_reBF.freqs_kHz,rateSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle',LEVELlines{LEVind})
   hold on
   LEGtextSlope{LEVind}=sprintf('%.f-%.f dB SPL',unit.Tone_reBF.levels_dBSPL(LEVind-1),unit.Tone_reBF.levels_dBSPL(LEVind));
end
semilogx(unit.Info.BF_kHz*ones(1,2),[-1e6 1e6],'k:')
semilogx([unit.Tone_reBF.freqs_kHz(1) unit.Tone_reBF.freqs_kHz(end)],[0 0],'k-')
ylabel('Rate Slope (sp/sec/dB)')
% ylim([-1 1]*max(ceil(max(abs(rateSlopes))/.02)*.02))
ymax_DeltaRate=10;  % 0 to 300 in 30 dB
ylim([-1 1]*ymax_DeltaRate)
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Synch
subplot(323)
for LEVind=1:length(unit.Tone_reBF.levels_dBSPL)
   semilogx(unit.Tone_reBF.freqs_kHz,unit.Tone_reBF.synch(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle',LEVELlines{LEVind})
   hold on
end
hleg=legend(LEGtext,4);
set(hleg,'FontSize',6)
semilogx(unit.Info.BF_kHz*ones(1,2),[0 1],'k:')
ylabel('Synch')
ylim([0 1])
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Synch Change
subplot(324)
for LEVind=2:length(unit.Tone_reBF.levels_dBSPL)
   synchSlopes(LEVind,:)=(unit.Tone_reBF.synch(LEVind,:)-unit.Tone_reBF.synch(LEVind-1,:))/diff(unit.Tone_reBF.levels_dBSPL(LEVind-1:LEVind));
   semilogx(unit.Tone_reBF.freqs_kHz,synchSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle',LEVELlines{LEVind})
   hold on
end
semilogx(unit.Info.BF_kHz*ones(1,2),[-1e6 1e6],'k:')
semilogx([unit.Tone_reBF.freqs_kHz(1) unit.Tone_reBF.freqs_kHz(end)],[0 0],'k-')
hlegSlope=legend(LEGtextSlope{2:length(unit.Tone_reBF.levels_dBSPL)},4);
set(hlegSlope,'FontSize',6)
ylabel('Synch Slope (1/dB)')
% ylim([-1 1]*max(ceil(max(abs(synchSlopes))/.02)*.02))
ymax_DeltaSynch=1/30;  % 0 to 1 in 30 dB
ylim([-1 1]*ymax_DeltaSynch)
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Phase
subplot(325)
for LEVind=1:length(unit.Tone_reBF.levels_dBSPL)
   semilogx(unit.Tone_reBF.freqs_kHz,unit.Tone_reBF.phase(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle',LEVELlines{LEVind})
   hold on
end
semilogx(unit.Info.BF_kHz*ones(1,2),[-pi pi],'k:')
ylabel('Phase (rad)')
xlabel('Frequency (kHz)')
ylim([-pi pi])
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)

%%%% Phase Change
subplot(326)
for LEVind=2:length(unit.Tone_reBF.levels_dBSPL)
   phaseLL=repmat(unit.Tone_reBF.phase(LEVind-1,:),3,1);
   phaseLL(2,:)=phaseLL(1,:)+2*pi;
   phaseLL(3,:)=phaseLL(3,:)-2*pi;
   phaseHL=repmat(unit.Tone_reBF.phase(LEVind,:),3,1);
   phaseDiffs=phaseHL-phaseLL;
   [phaseMinDiffs,phaseMinInds]=min(abs(phaseDiffs));
   for FREQind=1:length(unit.Tone_reBF.freqs_kHz)
      phaseSlopes(LEVind,FREQind)=phaseDiffs(phaseMinInds(FREQind),FREQind)/diff(unit.Tone_reBF.levels_dBSPL(LEVind-1:LEVind));
   end
   semilogx(unit.Tone_reBF.freqs_kHz,phaseSlopes(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
      LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle',LEVELlines{LEVind})
   hold on
end
semilogx(unit.Info.BF_kHz*ones(1,2),[-pi pi],'k:')
semilogx([unit.Tone_reBF.freqs_kHz(1) unit.Tone_reBF.freqs_kHz(end)],[0 0],'k-')
ylabel('Phase Slope (rad/dB)')
xlabel('Frequency (kHz)')
% ylim([-1 1]*max(ceil(max(abs(phaseSlopes))/.02)*.02))  
ymax_DeltaSynch=0.05;  % ~pi/2 in 30 dB
ylim([-1 1]*ymax_DeltaSynch)
xlim(unit.Info.BF_kHz*[.5*.99 2/.99])
set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
hold off
set(gca,'FontSize',FIG.FontSize)




return;

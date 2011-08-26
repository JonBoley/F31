function ExpPlot_TrFF(ExpDate)
% File: ExpPlot_TrFF.m
% Date: 18Sep2004 (M. Heinz) (Modified from UnitPlot_TrBF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Plots Tone_reFF data for a given experiment.  Loads 'DataList_ExpName.m' file, and steps through each unit.
% Separate plots for each FF (fixed frequency).  Each plot has rate, synch, phase, and associated slopes vs BF for a given FF
% UnitAnal_TrFF.m performs the relevant analysis.
%

global NOHR_dir NOHR_ExpList

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='090204'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir=fullfile(NOHR_dir,'ExpData');
anal_dir=fullfile(NOHR_dir,'Data Analysis');
stim_dir=fullfile(NOHR_dir,'Stimuli');

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
disp(sprintf('Plotting  (TreFF) Experiment: ''%s''',ExpName))

eval(['cd ''' fullfile(data_dir,ExpName) ''''])
FileName=strcat('DataList_',ExpDateText,'.mat');

% Load DataList file with all picture numbers
if length(dir(FileName))
   disp(['   *** Loading file: "' FileName '"'])
   eval(['load ' FileName])
else
   error(sprintf('%s: FILE NOT FOUND!!!!!',FileName))
end
[NumTracks,NumUnits]=size(DataList.Units);


%%%% Go through DataList first time to find all FFreqs_kHz and levels_dBSPL
FFreqs_kHz=[];
levels_dBSPL=[];
for TrackIND=1:NumTracks
   for UnitIND=1:NumUnits
      
      % Does this unit exist?
      if ~isempty(DataList.Units{TrackIND,UnitIND})
         % Is there Tone_reFF data for this unit?
         if isfield(DataList.Units{TrackIND,UnitIND},'Tone_reFF')
            
            % Store FF for this unit
            FFreqs_kHz=union(FFreqs_kHz,DataList.Units{TrackIND,UnitIND}.Tone_reFF.freqs_kHz);
            
            % Store levels for this unit
            levels_dBSPL=union(levels_dBSPL,DataList.Units{TrackIND,UnitIND}.Tone_reFF.levels_dBSPL);
            
         end
      end
   end
end
FFreqs_kHz=sort(FFreqs_kHz);
levels_dBSPL=sort(levels_dBSPL);
NumFFs=length(FFreqs_kHz);
NumLs=length(levels_dBSPL);

%%%% Create data structures to store Experiment Data.  Data stored in cell arrays, with each cell for a different FFreq
rate=cell(1,NumFFs);
BFs_kHz=cell(1,NumFFs);
Units=cell(1,NumFFs);
synch=cell(1,NumFFs);
phase=cell(1,NumFFs);
RaySig=cell(1,NumFFs);
% Create NaN matrices before filling in data, the truncate at end
for FFind=1:NumFFs
   rate{FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
   BFs_kHz{FFind}=NaN*ones(1,DataList.General.TOTALunits);
   Units{FFind}=cell(1,DataList.General.TOTALunits);
   synch{FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
   phase{FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
   RaySig{FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
end

%%%% Go back through each unit and store data
NUMunits={0,0,0};
for TrackIND=1:NumTracks
   for UnitIND=1:NumUnits
      
      % Does this unit exist?
      if ~isempty(DataList.Units{TrackIND,UnitIND})
         % Is there Tone_reFF data for this unit?
         if isfield(DataList.Units{TrackIND,UnitIND},'Tone_reFF')
            
            %%%% Load or creat unit structure for this unit
            UnitFileName=sprintf('unit.%d.%02d.mat',TrackIND,UnitIND);
            eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
            % If UNITSdata file does not exist, run analysis
            if isempty(ddd)
               UnitAnal_TrFF(ExpDate,DataList.Units{TrackIND,UnitIND}.Info.Unit,0);
            end
            eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

            % If Tone_reFF analysis is not completed, run here
            if ~isfield(unit.Tone_reFF,'synch')
               UnitAnal_TrFF(ExpDate,unit.Info.Unit,0);
               eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
            end
            
            % Find FF for this unit
            FFind=find(FFreqs_kHz==unit.Tone_reFF.freqs_kHz);
            NUMunits{FFind}=NUMunits{FFind}+1;
            
            BFs_kHz{FFind}(NUMunits{FFind})=unit.Info.BF_kHz;
            Units{FFind}{NUMunits{FFind}}=unit.Info.Unit;
            for LEVind=1:length(unit.Tone_reFF.levels_dBSPL)
               % Find levelIND for this 
               levelIND=find(levels_dBSPL==unit.Tone_reFF.levels_dBSPL(LEVind));
               rate{FFind}(levelIND,NUMunits{FFind})=unit.Tone_reFF.rate(LEVind,1);
               synch{FFind}(levelIND,NUMunits{FFind})=unit.Tone_reFF.synch(LEVind,1);
               phase{FFind}(levelIND,NUMunits{FFind})=unit.Tone_reFF.phase(LEVind,1);
               RaySig{FFind}(levelIND,NUMunits{FFind})=unit.Tone_reFF.RaySig(LEVind,1);
            end
         end
      end
   end
end

%%%% Cleanup data structures, to trim down to actual data
for FFind=1:NumFFs
   rate{FFind}=rate{FFind}(:,1:NUMunits{FFind});
   BFs_kHz{FFind}=BFs_kHz{FFind}(1:NUMunits{FFind});
   Units{FFind}=Units{FFind}(1:NUMunits{FFind});
   synch{FFind}=synch{FFind}(:,1:NUMunits{FFind});
   phase{FFind}=phase{FFind}(:,1:NUMunits{FFind});
   RaySig{FFind}=RaySig{FFind}(:,1:NUMunits{FFind});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FIG.FontSize=8;

%%%% Create data structures to store Experiment Data.  Data stored in cell arrays, with each cell for a different FFreq
rateSlopes=cell(1,NumFFs);
synchSlopes=cell(1,NumFFs);
phaseSlopes=cell(1,NumFFs);
% Create NaN matrices before filling in data, the truncate at end
warning off % don't need to see "divide by zero" error over and over
for FFind=1:NumFFs
   rateSlopes{FFind}=NaN+zeros(size(rate{FFind}));
   synchSlopes{FFind}=NaN+zeros(size(synch{FFind}));
   phaseSlopes{FFind}=NaN+zeros(size(phase{FFind}));
   %%%% Take out (mark with NaNs) insignificant phase and synch values
   synch{FFind}=synch{FFind}.*RaySig{FFind}./RaySig{FFind};
   phase{FFind}=phase{FFind}.*RaySig{FFind}./RaySig{FFind};
end
warning on

%%%% Close all old Figures from EHFF plotting
OPEN_FIGS=get(0,'Children');
eval('close(intersect(1:3,OPEN_FIGS))','')
for FFind=1:NumFFs
   
   xLIMS=[.1 20];
   
   figure(FFind); clf
   set(gcf,'pos',[420   199   976   761])
   %%%% Rate
   subplot(321)
   clear LEGtext
   for LEVind=1:NumLs
      semilogx(BFs_kHz{FFind},rate{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
      hold on
      LEGtext{LEVind}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
   end
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   %    ylim([0 ceil(max(max(rate{FFind}))/50)*50])
   %    xlim(FFreqs_kHz(FFind)*[.5 2])
   ylim([0 300])  % Fixed ordinate for all plots
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   %    set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
   ylabel('Rate (sp/sec)')
   ht1=title(sprintf('Exp: %s\n[ExpPlot_TONEreFF]  (FF=%.2f kHz)',unit.Info.ExpName,FFreqs_kHz(FFind)));
   set(ht1,'units','norm','Interpreter','none','pos',[1.2 1.0292 0])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% Rate Change
   subplot(322)
   clear LEGtextSlope
   for LEVind=2:NumLs
      rateSlopes{FFind}(LEVind,:)=(rate{FFind}(LEVind,:)-rate{FFind}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
      semilogx(BFs_kHz{FFind},rateSlopes{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind}, ...
         'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
      hold on
      LEGtextSlope{LEVind}=sprintf('%.f-%.f dB SPL',levels_dBSPL(LEVind-1),levels_dBSPL(LEVind));
   end
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   semilogx(xLIMS,[0 0],'k-')
   ylabel('Rate Slope (sp/sec/dB)')
   % ylim([-1 1]*max(ceil(max(abs(rateSlopes))/.02)*.02))
   ymax_DeltaRate=10;  % 0 to 300 in 30 dB
   ylim([-1 1]*ymax_DeltaRate)
   %    xlim(unit.Info.BF_kHz*[.5 2])
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% Synch
   subplot(323)
   for LEVind=1:NumLs
      semilogx(BFs_kHz{FFind},synch{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
      hold on
   end
   hleg=legend(LEGtext,4);
   set(hleg,'FontSize',6)
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   ylabel('Synch')
   ylim([0 1])
   %    xlim(unit.Info.BF_kHz*[.5 2])
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% Synch Change
   subplot(324)
   for LEVind=2:NumLs
      synchSlopes{FFind}(LEVind,:)=(synch{FFind}(LEVind,:)-synch{FFind}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
      semilogx(BFs_kHz{FFind},synchSlopes{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind}, ...
         'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
      hold on
   end
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   semilogx(xLIMS,[0 0],'k-')
   hlegSlope=legend(LEGtextSlope{2:NumLs},4);
   set(hlegSlope,'FontSize',6)
   ylabel('Synch Slope (1/dB)')
   % ylim([-1 1]*max(ceil(max(abs(synchSlopes))/.02)*.02))
   ymax_DeltaSynch=1/30;  % 0 to 1 in 30 dB
   ylim([-1 1]*ymax_DeltaSynch)
   %    xlim(unit.Info.BF_kHz*[.5 2])
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% Phase
   subplot(325)
   for LEVind=1:NumLs
      semilogx(BFs_kHz{FFind},phase{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
      hold on
   end
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   ylabel('Phase (rad)')
   xlabel('Frequency (kHz)')
   ylim([-pi pi])
   %    xlim(unit.Info.BF_kHz*[.5 2])
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
   %%%% Phase Change
   subplot(326)
   for LEVind=2:NumLs
      phaseLL=repmat(phase{FFind}(LEVind-1,:),3,1);
      phaseLL(2,:)=phaseLL(1,:)+2*pi;
      phaseLL(3,:)=phaseLL(3,:)-2*pi;
      phaseHL=repmat(phase{FFind}(LEVind,:),3,1);
      phaseDiffs=phaseHL-phaseLL;
      [phaseMinDiffs,phaseMinInds]=min(abs(phaseDiffs));
      for FREQind=1:length(BFs_kHz{FFind})
         phaseSlopes{FFind}(LEVind,FREQind)=phaseDiffs(phaseMinInds(FREQind),FREQind)/diff(levels_dBSPL(LEVind-1:LEVind));
      end
      semilogx(BFs_kHz{FFind},phaseSlopes{FFind}(LEVind,:),'Marker',LEVELmarkers{LEVind}, ...
         'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
      hold on
   end
   semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
   semilogx(xLIMS,[0 0],'k-')
   ylabel('Phase Slope (rad/dB)')
   xlabel('Frequency (kHz)')
   % ylim([-1 1]*max(ceil(max(abs(phaseSlopes))/.02)*.02))  
   %    ymax_DeltaSynch=0.05;  % ~pi/2 in 30 dB
   ymax_DeltaSynch=0.1;  % ~pi in 30 dB
   ylim([-1 1]*ymax_DeltaSynch)
   %    xlim(unit.Info.BF_kHz*[.5 2])
   xlim(xLIMS)
   set(gca,'XTickLabel',[.1 1 10])
   hold off
   set(gca,'FontSize',FIG.FontSize)
   
end

return;

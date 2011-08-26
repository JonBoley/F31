function ExpPlot_EHrFF(ExpDate)
% File: ExpPlot_EHrFF.m
% Date: 19Sep2004 (M. Heinz) (Modified from ExpPlot_TrFF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
%
% Plots EH_reFF data for a given experiment.  Loads 'DataList_ExpName.m' file, and steps through each unit.
% Separate plots for each FF (fixed frequency), each feature, and each Harmonic/Polarity combo (up to 4).  
% Each plot has rate, synch, phase, and associated slopes vs BF for a given FF.
% UnitAnal_EHrFF.m performs the relevant analysis.
%

global NOHR_dir NOHR_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

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
disp(sprintf('Plotting  (EHreFF) Experiment: ''%s''',ExpName))

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Go through DataList first time to find all FFreqs_kHz and levels_dBSPL
FFreqs_kHz=[];
levels_dBSPL=[];
FeatINDs_ALL=[];
% HarmonicINDs=[];
% PolarityINDs=[];
for TrackIND=1:NumTracks
   for UnitIND=1:NumUnits
      
      % Does this unit exist?
      if ~isempty(DataList.Units{TrackIND,UnitIND})
         % Is there EH_reFF data for this unit?
         if isfield(DataList.Units{TrackIND,UnitIND},'EH_reFF')
            
            % Find all features
            UnitFeats=fieldnames(DataList.Units{TrackIND,UnitIND}.EH_reFF);
            clear FeatINDs
            for FeatIND=1:length(UnitFeats)
               FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
            end
            % Store Feature_Indices for this unit
            FeatINDs_ALL=union(FeatINDs_ALL,FeatINDs);
            
            for FeatIND=FeatINDs  % Step through each Feature we have data for
               for HarmonicsIND=1:2
                  for PolarityIND=1:2
                     eval(['yTEMP=DataList.Units{TrackIND,UnitIND}.EH_reFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
                     if ~isempty(yTEMP)
%                         % Store Harmonic_Indices for this unit
%                         HarmonicINDs=union(HarmonicINDs,HarmonicsIND);
%                         % Store Polarity_Indices for this unit
%                         PolarityINDs=union(PolarityINDs,PolarityIND);
                        % Store FF for this unit
                        FFreqs_kHz=union(FFreqs_kHz,yTEMP.freqs_kHz);
                        % Store levels for this unit
                        levels_dBSPL=union(levels_dBSPL,yTEMP.levels_dBSPL);
                     end
                  end
               end
            end
         end
      end
   end
end
FFreqs_kHz=sort(FFreqs_kHz);
levels_dBSPL=sort(levels_dBSPL);
NumFFs=length(FFreqs_kHz);
NumLs=length(levels_dBSPL);

FeatINDs_ALL=sort(FeatINDs_ALL);
TOTALfeats=length(FeaturesText);
FeatPlotIND=zeros(size(FeaturesText));
% HarmonicINDs=sort(HarmonicINDs);
% PolarityINDs=sort(PolarityINDs);
% NumFeats=length(FeatINDs_ALL);
% NumHarms=length(HarmonicINDs);
% NumPolarity=length(PolarityINDs);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 9/19/04
% ALL SET UP TO HERE, so we just need to store separate data for each Feature, Harmonic, Polarity 
% so, I guess we'll just make rate=cell(TOTALfeats,NumHarms=2,NumPolarity=2,NumFFs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Create data structures to store Experiment Data.  Data stored in cell arrays, with each cell for a different FFreq
rate=cell(TOTALfeats,2,2,NumFFs);
BFs_kHz=cell(TOTALfeats,2,2,NumFFs);
Units=cell(TOTALfeats,2,2,NumFFs);
synch=cell(TOTALfeats,2,2,NumFFs);
phase=cell(TOTALfeats,2,2,NumFFs);
NUMunits=cell(TOTALfeats,2,2,NumFFs);
% Create NaN matrices before filling in data, the truncate at end
for FeatIND=FeatINDs_ALL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:NumFFs
            rate{FeatIND,HarmonicsIND,PolarityIND,FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
            BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind}=NaN*ones(1,DataList.General.TOTALunits);
            Units{FeatIND,HarmonicsIND,PolarityIND,FFind}=cell(1,DataList.General.TOTALunits);
            synch{FeatIND,HarmonicsIND,PolarityIND,FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
            phase{FeatIND,HarmonicsIND,PolarityIND,FFind}=NaN*ones(NumLs,DataList.General.TOTALunits);
         end
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Go back through each unit and store ALL data vs BF
for FeatIND=FeatINDs_ALL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:NumFFs
            NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind}=0;
         end
      end
   end
end

for TrackIND=1:NumTracks
   for UnitIND=1:NumUnits
      
      % Does this unit exist?
      if ~isempty(DataList.Units{TrackIND,UnitIND})
         % Is there EH_reFF data for this unit?
         if isfield(DataList.Units{TrackIND,UnitIND},'EH_reFF')
            
            %%%% Load or creat unit structure for this unit
            UnitFileName=sprintf('unit.%d.%02d.mat',TrackIND,UnitIND);
            eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
            % If UNITSdata file does not exist, run analysis
            if isempty(ddd)
               UnitAnal_EHrFF(ExpDate,DataList.Units{TrackIND,UnitIND}.Info.Unit,0);
            end
            eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

            % If EH_reFF analysis is not completed, run here
            UnitFeats=fieldnames(unit.EH_reFF);
            clear processedBOOL
            for HarmonicsIND=1:2
               for PolarityIND=1:2
                  eval(['emptyBOOL=~isempty(unit.EH_reFF.' UnitFeats{1} '{HarmonicsIND,PolarityIND});'])
                  if emptyBOOL
                     eval(['processedBOOL=isfield(unit.EH_reFF.' UnitFeats{1} '{HarmonicsIND,PolarityIND},''synch'');'])
                     if ~processedBOOL
                        UnitAnal_EHrFF(ExpDate,unit.Info.Unit,0);
                        eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
                     end
                     break
                  end
               end
               if exist('processedBOOL','var')
                  break
               end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Only store Synch and Phase for Feature of interest (for now)
            clear FeatINDs
            for FeatIND=1:length(UnitFeats)
               FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
            end
            
            for FeatIND=FeatINDs
               for HarmonicsIND=1:2
                  for PolarityIND=1:2
                     
                     eval(['yTEMP=unit.EH_reFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
                     if ~isempty(yTEMP)
                        
                        %%%% Take out (mark with NaNs) insignificant phase and synch values
                        warning off % don't need to see "divide by zero" error over and over
                        for LEVind=1:length(yTEMP.levels_dBSPL)
                           yTEMP.synch{LEVind}=yTEMP.synch{LEVind}.*yTEMP.RaySig{LEVind}./yTEMP.RaySig{LEVind};
                           yTEMP.phase{LEVind}=yTEMP.phase{LEVind}.*yTEMP.RaySig{LEVind}./yTEMP.RaySig{LEVind};
                        end
                        warning on
                        
                        %%%% Later, Can choose feature to plot; FOR NOW, assume only plot Feature at FF
                        FEATtoPLOT=FeaturesText{FeatIND};
                        FeatPlotIND(FeatIND)=find(strcmp(FeaturesText,FEATtoPLOT));
                        
                        % Find FF for this unit
                        FFind=find(FFreqs_kHz==yTEMP.freqs_kHz);
                        NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind}=NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind}+1;
                        
                        BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind}(NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind})=unit.Info.BF_kHz;
                        Units{FeatIND,HarmonicsIND,PolarityIND,FFind}{NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind}}=unit.Info.Unit;
                        for LEVind=1:length(yTEMP.levels_dBSPL)
                           % Find levelIND for this 
                           levelIND=find(levels_dBSPL==yTEMP.levels_dBSPL(LEVind));
                           rate{FeatIND,HarmonicsIND,PolarityIND,FFind}(levelIND,NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind})= ...
                              yTEMP.rate(LEVind,1);
                           if ~isempty(yTEMP.synch{LEVind,1})
                              synch{FeatIND,HarmonicsIND,PolarityIND,FFind}(levelIND,NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind})= ...
                                 yTEMP.synch{LEVind,1}(FeatPlotIND(FeatIND));
                              phase{FeatIND,HarmonicsIND,PolarityIND,FFind}(levelIND,NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind})= ...
                                 yTEMP.phase{LEVind,1}(FeatPlotIND(FeatIND));
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
end

%%%% Cleanup data structures, to trim down to actual data
for FeatIND=FeatINDs_ALL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:NumFFs
            rate{FeatIND,HarmonicsIND,PolarityIND,FFind}=rate{FeatIND,HarmonicsIND,PolarityIND,FFind}(:,1:NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind});
            BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind}=BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind}(1:NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind});
            Units{FeatIND,HarmonicsIND,PolarityIND,FFind}=Units{FeatIND,HarmonicsIND,PolarityIND,FFind}(1:NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind});
            synch{FeatIND,HarmonicsIND,PolarityIND,FFind}=synch{FeatIND,HarmonicsIND,PolarityIND,FFind}(:,1:NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind});
            phase{FeatIND,HarmonicsIND,PolarityIND,FFind}=phase{FeatIND,HarmonicsIND,PolarityIND,FFind}(:,1:NUMunits{FeatIND,HarmonicsIND,PolarityIND,FFind});
         end
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FIG.FontSize=8;

%%%% Create data structures to store Experiment Data.  Data stored in cell arrays, with each cell for a different FFreq
rateSlopes=cell(TOTALfeats,2,2,NumFFs);
synchSlopes=cell(TOTALfeats,2,2,NumFFs);
phaseSlopes=cell(TOTALfeats,2,2,NumFFs);
% Create NaN matrices before filling in data, the truncate at end
for FeatIND=FeatINDs_ALL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:NumFFs
            rateSlopes{FFind}=NaN+zeros(size(rate{FFind}));
            synchSlopes{FFind}=NaN+zeros(size(synch{FFind}));
            phaseSlopes{FFind}=NaN+zeros(size(phase{FFind}));
         end
      end
   end
end

%%%% Close all old Figures from EHFF plotting
ALL_EHfigNUMS=[];
for FeatIND=1:length(FeaturesText)
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:3
            ALL_EHfigNUMS=[ALL_EHfigNUMS str2num(strcat(num2str(FFind),num2str(FeatIND),num2str(HarmonicsIND),num2str(PolarityIND)))];
         end
      end
   end
end
OPEN_FIGS=get(0,'Children');
eval('close(intersect(ALL_EHfigNUMS,OPEN_FIGS))','')

for FeatIND=FeatINDs_ALL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         for FFind=1:NumFFs
            
            if ~isempty(rate{FeatIND,HarmonicsIND,PolarityIND,FFind})
               FIGnum=str2num(strcat(num2str(FFind),num2str(FeatIND),num2str(HarmonicsIND),num2str(PolarityIND)));
               
               xLIMS=[.1 20];
               
               figure(FIGnum); clf
               set(gcf,'pos',[420   199   976   761])
               %%%% Rate
               subplot(321)
               clear LEGtext
               for LEVind=1:NumLs
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},rate{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
                  hold on
                  LEGtext{LEVind}=sprintf('%.f dB SPL',levels_dBSPL(LEVind));
               end
               semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
               
               if FIGnum==1321
                  disp('STOP HERE')
               end
               %                ylim([0 ceil(max(max(rate{FeatIND,HarmonicsIND,PolarityIND,FFind}))/50)*50])
               %    xlim(FFreqs_kHz(FFind)*[.5 2])
               ylim([0 300])  % Fixed ordinate for all plots
               xlim(xLIMS)
               set(gca,'XTickLabel',[.1 1 10])
               %    set(gca,'XTick',unit.Tone_reBF.freqs_kHz)
               ylabel('Overall Rate (sp/sec)')
               ht1=title(sprintf('Exp: %s\n[ExpPlot_EHreFF]  (FF=%.2f kHz);    Feature: %s (plotting: %s);    Formants at Harmonics: %s;    Invert Polarity: %s', ...
                  unit.Info.ExpName,FFreqs_kHz(FFind),FeaturesText{FeatIND}, ...
                  FeaturesText{FeatPlotIND(FeatIND)},upper(FormsAtHarmonicsText{HarmonicsIND}),upper(InvertPolarityText{PolarityIND})));
               set(ht1,'units','norm','Interpreter','none','pos',[1.2 1.0292 0])
               
               hold off
               set(gca,'FontSize',FIG.FontSize)
               
               %%%% Rate Change
               subplot(322)
               clear LEGtextSlope
               for LEVind=2:NumLs
                  rateSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:)=(rate{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:) ...
                     -rate{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},rateSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
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
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},synch{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
                  hold on
               end
               hleg=legend(LEGtext,4);
               set(hleg,'FontSize',6)
               semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
               ylabel(sprintf('Synch (to %s)',FeaturesText{FeatPlotIND(FeatIND)}))
               ylim([0 1])
               %    xlim(unit.Info.BF_kHz*[.5 2])
               xlim(xLIMS)
               set(gca,'XTickLabel',[.1 1 10])
               hold off
               set(gca,'FontSize',FIG.FontSize)
               
               %%%% Synch Change
               subplot(324)
               for LEVind=2:NumLs
                  synchSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:)=(synch{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:)- ...
                     synch{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},synchSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
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
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},phase{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
                  hold on
               end
               semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
               ylabel(sprintf('Phase (to %s) (rad)',FeaturesText{FeatPlotIND(FeatIND)}))
               xlabel('Best Frequency (kHz)')
               ylim([-pi pi])
               %    xlim(unit.Info.BF_kHz*[.5 2])
               xlim(xLIMS)
               set(gca,'XTickLabel',[.1 1 10])
               hold off
               set(gca,'FontSize',FIG.FontSize)
               
               %%%% Phase Change
               subplot(326)
               for LEVind=2:NumLs
                  phaseLL=repmat(phase{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind-1,:),3,1);
                  phaseLL(2,:)=phaseLL(1,:)+2*pi;
                  phaseLL(3,:)=phaseLL(3,:)-2*pi;
                  phaseHL=repmat(phase{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:),3,1);
                  phaseDiffs=phaseHL-phaseLL;
                  [phaseMinDiffs,phaseMinInds]=min(abs(phaseDiffs));
                  for FREQind=1:length(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind})
                     phaseSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,FREQind)=phaseDiffs(phaseMinInds(FREQind),FREQind)/ ...
                        diff(levels_dBSPL(LEVind-1:LEVind));
                  end
                  semilogx(BFs_kHz{FeatIND,HarmonicsIND,PolarityIND,FFind},phaseSlopes{FeatIND,HarmonicsIND,PolarityIND,FFind}(LEVind,:), ...
                     'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor',LEVELcolors{LEVind},'Color',LEVELcolors{LEVind-1},'LineStyle','none')
                  hold on
               end
               semilogx(FFreqs_kHz(FFind)*ones(1,2),[-1e6 1e6],'k:')
               semilogx(xLIMS,[0 0],'k-')
               ylabel('Phase Slope (rad/dB)')
               xlabel('Best Frequency (kHz)')
               % ylim([-1 1]*max(ceil(max(abs(phaseSlopes))/.02)*.02))  
               %    ymax_DeltaSynch=0.05;  % ~pi/2 in 30 dB
               ymax_DeltaSynch=0.1;  % ~pi in 30 dB
               ylim([-1 1]*ymax_DeltaSynch)
               %    xlim(unit.Info.BF_kHz*[.5 2])
               xlim(xLIMS)
               set(gca,'XTickLabel',[.1 1 10])
               hold off
               set(gca,'FontSize',FIG.FontSize)
                  
               % orient landscape
               % print -dwinc
               
            end  % If data, plot
            
         end
      end
   end
end


return;

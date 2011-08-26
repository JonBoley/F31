function UnitsPlotCrossBF_EHrFF(ExpDate,Unit1Name,Unit2Name)
% File: UnitsPlotCrossBF_EHrFF.m
% Date: 25Sep2004 (M. Heinz) (modified from UnitsPlotCrossBF_TrBF)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitNames: '3.07' (converted later)
%
% Plots EH_reFF data for 2 units from a given experiment PLUS the same data from SCC Analysis of Cross-Unit Correlation.  
% Loads 'UNITSdata/unit.T.U.mat' file.
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
      ExpDate=input('Enter Experiment Date (e.g., ''090204''): ');
   end
end
if ~exist('Unit1Name','var')
   Unit1Name=0;
   %%% HARD CODE FOR NOW
Unit1Name='1.18'
   while ~ischar(Unit1Name)
      Unit1Name=input('Enter Unit 1 Name (e.g., ''3.07''): ');
   end
end
if ~exist('Unit2Name','var')
   Unit2Name=0;
   %%% HARD CODE FOR NOW
Unit2Name='1.26'
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
disp(sprintf('Plotting Cross_BF Analysis (EHreFF) Experiment: ''%s''; Units: %d.%02d and %d.%02d',ExpName,TrackNums{1},UnitNums{1},TrackNums{2},UnitNums{2}))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
units=cell(1,2);
UnitFileNames=cell(1,2);
UnitFeats=cell(1,2);
%%%% Load unit structures for these units
for i=1:2
   UnitFileNames{i}=sprintf('unit.%d.%02d.mat',TrackNums{i},UnitNums{i});
   eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileNames{i}) ''');'])
   % If UNITSdata file does not exist, run analysis
   if isempty(ddd)
      UnitAnal_EHrFF(ExpDate,UnitNames{i},0);
   end
   eval(['load ''' fullfile('UNITSdata',UnitFileNames{i}) ''''])
   
   % Make sure there is EH_reFF data
   if ~isfield(unit,'EH_reFF')
      disp(sprintf('*** No EH_reFF data for unit %s!!!',UnitNames{i}))
      beep
      return;
   end
   
   % If EH_reFF analysis is not completed, run here
   UnitFeats{i}=fieldnames(unit.EH_reFF);
   clear processedBOOL
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         eval(['emptyBOOL=~isempty(unit.EH_reFF.' UnitFeats{i}{1} '{HarmonicsIND,PolarityIND});'])
         if emptyBOOL
            eval(['processedBOOL=isfield(unit.EH_reFF.' UnitFeats{i}{1} '{HarmonicsIND,PolarityIND},''synch'');'])         
            if ~processedBOOL
               UnitAnal_EHrFF(ExpDate,UnitNames{i},0);   
               eval(['load ''' fullfile('UNITSdata',UnitFileNames{i}) ''''])
            end
            break
         end
      end
      if exist('processedBOOL','var')
         break
      end
   end
   
   units{i}=unit;
   clear unit   
   
end

%%%% Plot data
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FIG.FontSize=8;

%%%%%%%%%%%%% 9/25/04  TO DO   %%%%%%%%%%%%%%%%%%%%
% *1) Basics for 2 units are steup for TrFF
% *2) Found 2 good units
% *3) NOW, need to do SCC analysis on these two units, and plot that data HERE!!!!
% *4) Do same for EHrFF
% *5) ONTO reBF!!!!!!!! & simulate 2 neurons



%%% Separate Plot for each Condition
%%% Only plot Synch and Phase for Feature of interest (for now)
clear FeatINDs
FeatINDs=cell(1,3);
FeatPlotIND=zeros(size(FeaturesText));
for i=1:2
   for FeatIND=1:length(UnitFeats{i})
      FeatINDs{i}(FeatIND)=find(strcmp(FeaturesText,UnitFeats{i}{FeatIND}));
   end
end
FeatINDs{3}=intersect(FeatINDs{1},FeatINDs{2});
TOTALfeats=length(FeaturesText);

%%%% Store data from 2 units and SCC analysis
BFs_kHz=NaN+ones(1,5);
rate=cell(TOTALfeats,2,2);
synch=cell(TOTALfeats,2,2);
phase=cell(TOTALfeats,2,2);
rateSlopes=cell(TOTALfeats,2,2);
synchSlopes=cell(TOTALfeats,2,2);
phaseSlopes=cell(TOTALfeats,2,2);

for FeatIND=FeatINDs{3}
   for HarmonicsIND=1:2
      for PolarityIND=1:2
         yTEMP=cell(1,2);
         
         eval(['yTEMP{1}=units{1}.EH_reFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         eval(['yTEMP{2}=units{2}.EH_reFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         if (~isempty(yTEMP{1})&~isempty(yTEMP{2}))
            
            levels_dBSPL=union(yTEMP{1}.levels_dBSPL,yTEMP{2}.levels_dBSPL);
            
            %%%% Later, Can choose feature to plot; FOR NOW, assume only plot Feature at FF
            FEATtoPLOT=FeaturesText{FeatIND};
            FeatPlotIND(FeatIND)=find(strcmp(FeaturesText,FEATtoPLOT));
            
            FIGnum=str2num(strcat(num2str(FeatIND),num2str(HarmonicsIND),num2str(PolarityIND)));
            
            %%%% Store data from 2 units and SCC analysis
            rate{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),5);
            synch{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),3);
            phase{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),3);
            rateSlopes{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),5);
            synchSlopes{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),3);
            phaseSlopes{FeatIND,HarmonicsIND,PolarityIND}=NaN+zeros(length(levels_dBSPL),3);
            
            % Store rate, synch, and phase data from 2 units
            for i=1:2
               BFs_kHz(i)=units{i}.Info.BF_kHz;
               if i==1
                  FF_kHz=yTEMP{i}.freqs_kHz;
               else
                  if yTEMP{i}.freqs_kHz~=FF_kHz
                     error('Difference in Fixed Frequency between 2 units')
                  end
               end
               
               %%%% Take out (mark with NaNs) insignificant phase and synch values
               warning off % don't need to see "divide by zero" error over and over
               for LEVind=1:length(yTEMP{i}.levels_dBSPL)
                  for FREQind=1:length(yTEMP{i}.freqs_kHz)
                     yTEMP{i}.synch{LEVind,FREQind}=yTEMP{i}.synch{LEVind,FREQind}.* ...
                        yTEMP{i}.RaySig{LEVind,FREQind}./yTEMP{i}.RaySig{LEVind,FREQind};
                     yTEMP{i}.phase{LEVind,FREQind}=yTEMP{i}.phase{LEVind,FREQind}.* ...
                        yTEMP{i}.RaySig{LEVind,FREQind}./yTEMP{i}.RaySig{LEVind,FREQind};
                  end
               end
               warning on
               
               for LEVind=1:length(yTEMP{i}.levels_dBSPL)
                  rate{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==yTEMP{i}.levels_dBSPL(LEVind)),i)= ...
                     yTEMP{i}.rate(LEVind,1);
                  if ~isempty(yTEMP{i}.synch{LEVind,1})
                     synch{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==yTEMP{i}.levels_dBSPL(LEVind)),i)= ...
                        yTEMP{i}.synch{LEVind,1}(FeatPlotIND(FeatIND));
                     phase{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==yTEMP{i}.levels_dBSPL(LEVind)),i)= ...
                        yTEMP{i}.phase{LEVind,1}(FeatPlotIND(FeatIND));
                  end
               end
            end
            BFs_kHz(3)=geomean(BFs_kHz(1:2))/1.05;   % NSCC_0
            BFs_kHz(4)=geomean(BFs_kHz(1:2));        % NSCC_HLCD
            BFs_kHz(5)=geomean(BFs_kHz(1:2))*1.05;   % NSCC_max
            
            
            %%%%%%%%%%%%%%%%%% Calculate SCC
            levels_dBSPL_SCC=intersect(yTEMP{1}.levels_dBSPL,yTEMP{2}.levels_dBSPL);
            
            calcSAC=0;

            beep;
            input(sprintf('**********\nNEED TO MODIFY HLCD TO BE CALCULATED FROM F1, and used for all features,\n   not calculated for each feature\n**********'));
            
            HLCD_us=0;
            ADHOC_RATE=50;
            for LEVind=length(levels_dBSPL_SCC):-1:1  % go from high to low, to get HLCD first
               %    for LEVind=1:1
               
               for i=1:2
                  LEVinds(i)=find(yTEMP{i}.levels_dBSPL==levels_dBSPL(LEVind));
                  SCCtitles{i}=sprintf('UNIT %d: %s  [BF=%.4f kHz, Thr= %.f dB atten]', ...
                     i,units{i}.Info.Unit,units{i}.Info.BF_kHz,units{i}.Info.Threshold_dBatten);
               end
               SCCtitles{3}=sprintf('Exp: %s;     [SCCdemo_EHreFF]  (FF=%.2f kHz)   (%.f dB SPL)\nFeature: %s (plot: %s);  Formants at Harms: %s;  Invert Polarity: %s', ...
                  units{1}.Info.ExpName,FF_kHz,levels_dBSPL_SCC(LEVind),FeaturesText{FeatIND},FeaturesText{FeatPlotIND(FeatIND)}, ...
                  upper(FormsAtHarmonicsText{HarmonicsIND}),upper(InvertPolarityText{PolarityIND}));
               
               if (~isempty(yTEMP{1}.picNums{LEVinds(1)}))&(~isempty(yTEMP{2}.picNums{LEVinds(2)}))
                  disp(sprintf('%s\nProcessing: %.f dB SPL\n%s',repmat('*',1,21),levels_dBSPL_SCC(LEVind),repmat('*',1,21)))
                  
                  
                  [CD_us,NSCC_0,NSCC_CD,NSCC_HLCD,NSCC_max]= ...
                     SCCdemo_EHrFF(yTEMP{1}.picNums{LEVinds(1)},yTEMP{2}.picNums{LEVinds(2)}, ...
                     yTEMP{1}.excludeLines{LEVinds(1)},yTEMP{2}.excludeLines{LEVinds(2)},HLCD_us,SCCtitles,calcSAC);
                  
                  if LEVind==length(levels_dBSPL_SCC)
                     HLCD_us=CD_us;  % Set HLCD to be used for rest
                     NSCC_HLCD=NSCC_CD;
                  end
                  
                  rate{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),3)=NSCC_0 * ADHOC_RATE;
                  rate{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),4)=NSCC_HLCD * ADHOC_RATE;
                  rate{FeatIND,HarmonicsIND,PolarityIND}(find(levels_dBSPL==levels_dBSPL_SCC(LEVind)),5)=NSCC_max * ADHOC_RATE;
                                    
                  input('Enter')
               end
            end
            
            figure(FIGnum); clf
            set(gcf,'pos',[420   199   976   761])
            %%%% Rate
            subplot(321)
            clear LEGtext
            for LEVind=1:length(levels_dBSPL)
               semilogx(BFs_kHz,rate{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
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
            ht1=title(sprintf('Exp: %s;  [UnitsPlotCrossBF_EHreFF] (FF=%.2f kHz)\nFeature: %s (plot: %s);  Formants at Harms: %s;  Invert Polarity: %s\nUNIT 1: %s  [BF=%.4f kHz, Thr= %.f dB atten] %40s UNIT 2: %s  [BF=%.4f kHz, Thr= %.f dB atten]', ...
               units{1}.Info.ExpName,FF_kHz,FeaturesText{FeatIND},FeaturesText{FeatPlotIND(FeatIND)}, ...
               upper(FormsAtHarmonicsText{HarmonicsIND}),upper(InvertPolarityText{PolarityIND}),units{1}.Info.Unit, ...
               units{1}.Info.BF_kHz,units{1}.Info.Threshold_dBatten,' ',units{2}.Info.Unit,units{2}.Info.BF_kHz, ...
               units{2}.Info.Threshold_dBatten));
            set(ht1,'units','norm','Interpreter','none','pos',[1.2 1.0292 0])
            hold off
            set(gca,'FontSize',FIG.FontSize)
            
            %%%% Rate Change
            subplot(322)
            clear LEGtextSlope
            for LEVind=2:length(levels_dBSPL)
               rateSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:)=(rate{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:)-rate{FeatIND,HarmonicsIND,PolarityIND}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
               semilogx(BFs_kHz,rateSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
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
               semilogx(BFs_kHz(1:3),synch{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind},'Color',LEVELcolors{LEVind},'LineStyle','none')
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
               synchSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:)=(synch{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:)-synch{FeatIND,HarmonicsIND,PolarityIND}(LEVind-1,:))/diff(levels_dBSPL(LEVind-1:LEVind));
               semilogx(BFs_kHz(1:3),synchSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
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
               semilogx(BFs_kHz(1:3),phase{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind}, ...
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
               phaseLL=repmat(phase{FeatIND,HarmonicsIND,PolarityIND}(LEVind-1,:),3,1);
               phaseLL(2,:)=phaseLL(1,:)+2*pi;
               phaseLL(3,:)=phaseLL(3,:)-2*pi;
               phaseHL=repmat(phase{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),3,1);
               phaseDiffs=phaseHL-phaseLL;
               [phaseMinDiffs,phaseMinInds]=min(abs(phaseDiffs));
               for FREQind=1:size(phase,2)
                  phaseSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,FREQind)=phaseDiffs(phaseMinInds(FREQind),FREQind)/diff(levels_dBSPL(LEVind-1:LEVind));
               end
               semilogx(BFs_kHz(1:3),phaseSlopes{FeatIND,HarmonicsIND,PolarityIND}(LEVind,:),'Marker',LEVELmarkers{LEVind},'MarkerEdgeColor', ...
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
            
            
         end
      end
   end
end

return;

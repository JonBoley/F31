function UnitPlot_SMP(ExpDate,UnitName)
% File: UnitPlot_SMP.m  
% Date: 24Sep2004 (M. Heinz) (Modified from UnitPlot_EHrBF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '070804' (converted later)
% UnitName: '3.07' (converted later)
%
% Plots SMP analysis of EH_reBF data for a given experiment and unit.  Loads 'UNITSdata/unit.T.U.mat' file.
% Only looks at at-BF data, replicates May et al analysis.
% UnitAnal_EHrBF.m performs the relevant analysis.
%

global NOHR_dir NOHR_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText
   
%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='080204'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''070804''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
%%% HARD CODE FOR NOW
UnitName='1.07'
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
disp(sprintf('Plotting (SMP) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
% If UNITSdata file does not exist, run analysis
if isempty(ddd)
   UnitAnal_EHrBF(ExpDate,UnitName,0);
end
eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

% Make sure there is EH_reBF data
if ~isfield(unit,'EH_reBF')
   disp(sprintf('*** No EH_reBF data for this unit!!!'))
   beep
   return;
end

% If EH_reBF analysis is not completed, run here
UnitFeats=fieldnames(unit.EH_reBF);
UnitFeats=UnitFeats(~strcmp(UnitFeats,'interleaved'));  %% 010705: M Heinz; takes out newly added "interleaved" field
clear processedBOOL
for HarmonicsIND=1:2
   for PolarityIND=1:2
      eval(['emptyBOOL=~isempty(unit.EH_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND});'])
      if emptyBOOL
         eval(['processedBOOL=isfield(unit.EH_reBF.' UnitFeats{1} '{HarmonicsIND,PolarityIND},''synch'');'])         
         if ~processedBOOL
            UnitAnal_EHrBF(ExpDate,UnitName,0);   
            eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
         end
         break
      end
   end
   if exist('processedBOOL','var')
      break
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot data
LEVELmarkers={'s','^','*'};
LEVELlines={'-','--',':'};
LEVELcolors={'b','r','g'};

FEATmarkers={'.','x','s','^','*','<','^','>'};
FEATlines={':','-.','-','--',':','-.','-','--'};
FEATcolors={'m','c','b','r','g','m','c','b'};

FIG.FontSize=8;

%%%% Close all old Figures from EHFF plotting
ALL_EHfigNUMS=[11 12 21 22];
OPEN_FIGS=get(0,'Children');
eval('close(intersect(ALL_EHfigNUMS,OPEN_FIGS))','')

clear FeatINDs
for FeatIND=1:length(UnitFeats)
   FeatINDs(FeatIND)=find(strcmp(FeaturesText,UnitFeats{FeatIND}));
end
F1IND=find(strcmp(FeaturesText,'F1'));

%%% Separate Plot for each FormantsAtHarmonics and InvertPolarity Condition
for HarmonicsIND=1:2
   for PolarityIND=1:2
      yTEMP=cell(1,length(FeaturesText));
      FeatureRatevsOAL=cell(1,length(FeaturesText));
      OALevels_dBSPL=[];
      
      ANYdata=0;
      for FeatIND=FeatINDs
         eval(['yTEMP{FeatIND}=unit.EH_reBF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
         if ~isempty(yTEMP{FeatIND})
            ANYdata=1;
            OALevels_dBSPL=union(OALevels_dBSPL,yTEMP{FeatIND}.levels_dBSPL);
         end
      end
      OALevels_dBSPL=sort(OALevels_dBSPL);
      
      if ANYdata  % For this HarmIND/PolIND condition
         
         EffectiveBFs=NaN+ones(size(FeaturesText));
         for FeatIND=FeatINDs
            if ~isempty(yTEMP{FeatIND})
               FeatureRatevsOAL{FeatIND}=NaN+ones(size(OALevels_dBSPL));
               BFind=find(round(yTEMP{FeatIND}.freqs_kHz*1000)==round(unit.Info.BF_kHz*1000)); % frequency corresponding to BF
               if isempty(BFind)
                  error(sprintf('NO DATA at BF for Unit: %s',unit.Info.Unit))
               end
               for LEVind=1:length(yTEMP{FeatIND}.levels_dBSPL)
                  Lind=find(OALevels_dBSPL==yTEMP{FeatIND}.levels_dBSPL(LEVind));
                  FeatureRatevsOAL{FeatIND}(Lind)=yTEMP{FeatIND}.rate(LEVind,BFind)';
               end
               if FeatIND==FeatINDs(1)
                  FeatureLevels=yTEMP{FeatIND}.FeatureLevels_dB;
%                   FeatureFreqs_kHz=yTEMP{FeatIND}.FeatureFreqs_Hz/1000;
               else
                  if sum(yTEMP{FeatIND}.FeatureLevels_dB(find(~isnan(FeatureLevels)))~=FeatureLevels(find(~isnan(FeatureLevels))))
                     error('FeatureLevels-dB MISMATCH across FeatINDs');
                  end
               end
               % Find Effective BF for this feature (assumes F1 is at BF)
               EffectiveBFs(FeatIND)=yTEMP{FeatIND}.FeatureFreqs_Hz{BFind}(FeatIND)*(unit.Info.BF_kHz/yTEMP{FeatIND}.FeatureFreqs_Hz{BFind}(F1IND));
            end
         end
         
         %%
         % Data for Plot 2: Rate vs. Feature Level
         RatevsFeatureLevel=cell(size(OALevels_dBSPL));
         SensitivitySlopes=NaN+ones(size(OALevels_dBSPL));
         for LEVind=1:length(OALevels_dBSPL)
            RatevsFeatureLevel{LEVind}=NaN+ones(size(FeaturesText));
            for FeatIND=FeatINDs
               if ~isempty(FeatureRatevsOAL{FeatIND})
                  RatevsFeatureLevel{LEVind}(FeatIND)=FeatureRatevsOAL{FeatIND}(LEVind);
               end
            end
            if sum(~isnan(RatevsFeatureLevel{LEVind}))>1
               x=OALevels_dBSPL(LEVind)+FeatureLevels(find(~isnan(RatevsFeatureLevel{LEVind})));
               y=RatevsFeatureLevel{LEVind}(find(~isnan(RatevsFeatureLevel{LEVind})));
               [Cfit,MSE,fit]=fit1slope(x,y);
               SensitivitySlopes(LEVind)=Cfit(1);
            end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % 9/24/04 TODO
         % *1) SMP: Organize data for plots 2 and 3, and plot
         %    *? WHY IS featureLevels all less than 0???  - OK, its F0 that has the peak, and is at 0 dB
         % 9/25/04 TODO
         % *2) Look at some units, to get basic feeling
         %    *- getting warning for 80 dB SPL, when not present for Pol=---
         %    *- look at a few examples to be sure we've got most of the bugs
         %          - see UnitSMP notes for list of good units (TRY 080204, 1.29 for design of analyses)
         %
         % 3) Setup SAC/SCC computations, FIGURE OUT HOW TO
         % 4) Try to show level-dependent cross-BF Corr and enhanced Spectral Coding from TreFF and EHrFF
         %
         % 5) Demo 2-neuron/same-stim simulation with PSTview and real data!, FIGURE OUT HOW TO
         % 6) DESIGN general way to organize data with shifted SMPs and SCC to look at basic questions!!!
         % 7) Look for enhanced spectral coding in individual units
         %
         % 8) FOR ARO:
         %    - method for simulating 2 neurons, with same stimulus
         %    - evidence for NL phase effects for vowels
         %    - general increased cross-BF correlation with Level
         %    - any evidence for enhance spectral coding based on X-BF correlation
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         
         FIGnum=str2num(strcat(num2str(HarmonicsIND),num2str(PolarityIND)));
         
         figure(FIGnum); clf
         set(gcf,'pos',[753   199   643   761])
         %%%% Rate-Levels for each Feature
         subplot(311)
         LEGtext='';
         for FeatIND=FeatINDs
            if ~isempty(FeatureRatevsOAL{FeatIND})
               plot(OALevels_dBSPL,FeatureRatevsOAL{FeatIND},'Marker',FEATmarkers{FeatIND},'Color',FEATcolors{FeatIND},'LineStyle',FEATlines{FeatIND})
               hold on
               LEGtext{length(LEGtext)+1}=sprintf('%s',FeaturesText{FeatIND});
            end
         end
         ylim([0 300])  % Fixed ordinate for all plots
         xlim([00 100])
         ylabel('Rate (sp/sec)')
         xlabel('Overall Vowel Level (dB SPL)')
         ht1=title(sprintf('Exp: %s; Unit: %s [BF=%.4f kHz, Thr= %.f dB atten]\n[UnitPlot_SMP]   Formants at Harmonics: %s;    Invert Polarity: %s', ...
            unit.Info.ExpName,unit.Info.Unit,unit.Info.BF_kHz,unit.Info.Threshold_dBatten, ...
            upper(FormsAtHarmonicsText{HarmonicsIND}),upper(InvertPolarityText{PolarityIND})));
         set(ht1,'units','norm','Interpreter','none')
         hleg=legend(LEGtext,4);
         set(hleg,'FontSize',6)
         hold off
         set(gca,'FontSize',FIG.FontSize)

         %%%% Rate-FeatureLevels for each OAL
         subplot(312)
         LEGtext='';
         for LEVind=1:length(OALevels_dBSPL)
            if sum(~isnan(RatevsFeatureLevel{LEVind}))
               plot(OALevels_dBSPL(LEVind)+FeatureLevels,RatevsFeatureLevel{LEVind},'Marker','none','Color',LEVELcolors{LEVind},'LineStyle',LEVELlines{LEVind})
               hold on
               LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL (slope=%.2f sp/sec/dB)',OALevels_dBSPL(LEVind),SensitivitySlopes(LEVind));
            end
         end
         for FeatIND=FeatINDs
            if ~isempty(FeatureRatevsOAL{FeatIND})
               plot(OALevels_dBSPL+FeatureLevels(FeatIND),FeatureRatevsOAL{FeatIND},'Marker',FEATmarkers{FeatIND},'Color',FEATcolors{FeatIND},'LineStyle','none')
            end
         end
         ylim([0 300])  % Fixed ordinate for all plots
         xlim([0 100])
         ylabel('Rate (sp/sec)')
         xlabel('Feature Level (dB SPL)')
         hleg=legend(LEGtext,4);
         set(hleg,'FontSize',8)
         hold off
         set(gca,'FontSize',FIG.FontSize)

         %%%% Effective Spectrum Plots for each OAL
         subplot(313)
         LEGtext='';
         for LEVind=1:length(OALevels_dBSPL)
            if sum(~isnan(RatevsFeatureLevel{LEVind}))
               semilogx(EffectiveBFs,RatevsFeatureLevel{LEVind},'Marker','none','Color',LEVELcolors{LEVind},'LineStyle',LEVELlines{LEVind})
               hold on
               LEGtext{length(LEGtext)+1}=sprintf('%.f dB SPL',OALevels_dBSPL(LEVind));
            end
         end
         for FeatIND=FeatINDs
            semilogx(EffectiveBFs(FeatIND)*ones(size(FeatureRatevsOAL{FeatIND})),FeatureRatevsOAL{FeatIND},'Marker',FEATmarkers{FeatIND},'Color',FEATcolors{FeatIND},'LineStyle','none')
         end
         ylim([0 300])  % Fixed ordinate for all plots
         xlim([.1 40])
         ylabel('Rate (sp/sec)')
         xlabel('Effective BF (kHz)')
         hleg=legend(LEGtext,3);
         set(hleg,'FontSize',8)
         hold off
         set(gca,'FontSize',FIG.FontSize,'XTickLabels',[.1 1 10])
       
         
         orient tall
         
      end %End if data for this condition, plot
      
      
   end % End Polarity
end % End Harmonics


return;

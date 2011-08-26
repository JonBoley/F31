function UnitPlot_SpatTemp(ExpDate,UnitName)
% File: UnitPlot_SpatTemp.m
% Date: 01Nov2004 (M. Heinz) (Modified from UnitPlot_EHrBF_simFF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '080204' (converted later)
% UnitName: '1.29' (converted later)
%
% Plots Simulated spatio-temporal response pattern (ala Shamma 1985) from EH_reBF_simFF and T_reBF_simFF data 
% for a given experiment and unit.  Loads 'UNITSdata/unit.T.U.mat' file.
% ???? UnitAnal_EHrBF_simFF.m performs the relevant analysis.
% will need both EHreBF and TrBF data!!!!!!!!!!!!!!!!!
%
% 11/2/04 TO DO
% 1) only plots Fn/Tn{1,1} for now, need to add rest later
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12/28/04 TODO
% - Plot EHrBFi data
% - verify consistent with Shamma I
% - setup adjacent SCCs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

global NOHR_dir NOHR_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText
   
%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
%%% HARD CODE FOR NOW
ExpDate='080204'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''080204''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
%%% HARD CODE FOR NOW
UnitName='1.29'
   while ~ischar(UnitName)
      UnitName=input('Enter Unit Name (e.g., ''1.29''): ');
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
disp(sprintf('Plotting (Spatio-Temporal) Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
% If UNITSdata file does not exist, run analysis
if isempty(ddd)
   UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);
end
eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

% Make sure there is Tone_reBF or EH_reBF data
if (~isfield(unit,'EH_reBF'))|(~isfield(unit,'Tone_reBF'))
   disp(sprintf('*** No Tone_reBF or EH_reBF data for this unit!!!'))
   beep
   return;
end

% If EH_reBF_simFF analysis is not completed, run here
if ~isfield(unit,'EH_reBF_simFF')
   UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);   
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end
         
% If Tone_reBF_simFF analysis is not completed, run here
if ~isfield(unit,'Tone_reBF_simFF')
   UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

%%%% Find number of features to plot (tone, F1, T1, ...)
%%%% Find all levels and BFs to plot
NUMcols=4;
NUMrows=0;
levels_dBSPL=[]; BFsTEMP_kHz=[];
if isfield(unit,'Tone_reBF_simFF')
   NUMrows=NUMrows+1;
   levels_dBSPL=union(levels_dBSPL,unit.Tone_reBF_simFF.levels_dBSPL);
   BFsTEMP_kHz=union(BFsTEMP_kHz,unit.Tone_reBF_simFF.BFs_kHz);
end
if isfield(unit,'EH_reBF_simFF')
   EHfeats=fieldnames(unit.EH_reBF_simFF);
   NUMrows=NUMrows+length(EHfeats);
   clear FeatINDs
   for FeatIND=1:length(EHfeats)
      FeatINDs(FeatIND)=find(strcmp(FeaturesText,EHfeats{FeatIND}));
   end
   
   for FeatIND=FeatINDs
      for HarmonicsIND=1:2
         for PolarityIND=1:2
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               levels_dBSPL=union(levels_dBSPL,yTEMP.levels_dBSPL);
               BFsTEMP_kHz=union(BFsTEMP_kHz,yTEMP.BFs_kHz);
            end
         end
      end
   end
end
lowBF=min(BFsTEMP_kHz);
YLIMocts=2.5;  % Number of octaves to show on oridinate
clear BFsTEMP_kHz;

BFs_kHz=cell(NUMrows,length(levels_dBSPL));
PERhists=cell(NUMrows,length(levels_dBSPL));
PERhistXs_sec=cell(NUMrows,length(levels_dBSPL));
Nsps=cell(NUMrows,length(levels_dBSPL));
Rates=cell(NUMrows,length(levels_dBSPL));
Synchs=cell(NUMrows,length(levels_dBSPL));
Phases=cell(NUMrows,length(levels_dBSPL));

for LEVEL=levels_dBSPL
   figure(round(LEVEL)); clf
   set(gcf,'pos',[420     4   977   976])
   
   %%%%%%%%%%%%%%%%%%%% Tone Plots
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      PERhists{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      PERhistXs_sec{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      PSTxmax_sec=-999;
      for BFind=1:length(unit.Tone_reBF_simFF.BFs_kHz)   
         if ~isempty(unit.Tone_reBF_simFF.picNums{LEVind,BFind})
            
            PIC=concatPICS_NOHR(unit.Tone_reBF_simFF.picNums{LEVind,BFind},unit.Tone_reBF_simFF.excludeLines{LEVind,BFind});
            % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
            PIC=simFF_PICshift(PIC);
            PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
            
            BFs_kHz{ROWind,LEVind}(BFind)=unit.Tone_reBF_simFF.BFs_kHz(BFind);
            PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
            PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;
            if max(PERhistXs_sec{ROWind,LEVind}{BFind})>PSTxmax_sec
               PSTxmax_sec=max(PERhistXs_sec{ROWind,LEVind}{BFind});
            end
            Nsps{ROWind,LEVind}(BFind)=PIC.PERhist.NumDrivenSpikes;
            Rates{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);
            if PIC.SynchRate_PERhist.FeatureRaySig
               Synchs{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeatureSynchs;
               Phases{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeaturePhases;
            else
               Synchs{ROWind,LEVind}(BFind)=NaN;
               Phases{ROWind,LEVind}(BFind)=NaN;
            end
         else
            BFs_kHz{ROWind,LEVind}(BFind)=unit.Tone_reBF_simFF.BFs_kHz(BFind);
            Nsps{ROWind,LEVind}(BFind)=NaN;
            Rates{ROWind,LEVind}(BFind)=NaN;
            Synchs{ROWind,LEVind}(BFind)=NaN;
            Phases{ROWind,LEVind}(BFind)=NaN;
         end
      end
      
      %%%% Spatio-Temporal Plots
      NORMrate=1000;
      TFiltWidth=31;
      %       figure(round(unit.Tone_reBF_simFF.levels_dBSPL(LEVind)))
      PLOTnum=(ROWind-1)*NUMcols+1;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave
      for BFind=1:length(BFs_kHz{ROWind,LEVind})   
         if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
            LINEwidth=2;
         else
            LINEwidth=.5;
         end
         % ?? Is it OK to plot this as log, or do we need to make it linear?  - should be OK visually
         semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000,trifilt(PERhists{ROWind,LEVind}{BFind}/NORMrate,TFiltWidth)+BFs_kHz{ROWind,LEVind}(BFind),'LineWidth',LINEwidth)
         hold on
      end
      xlabel('Time (ms)')
      ylabel('Best Frequency (kHz)')
      if isempty(unit.Info.Threshold_dBatten)
         unit.Info.Threshold_dBatten=NaN;
      end
      if isempty(unit.Info.SR_sps)
         unit.Info.SR_sps=NaN;
      end
      if isempty(unit.Info.Q10)
         unit.Info.Q10=NaN;
      end
      title(sprintf('     Exp%s, %s: BF=%.2f kHz, Thr=%.f dB Atten, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL', ...
         ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBatten,unit.Info.SR_sps,unit.Info.Q10,'TONE', ...
         unit.Tone_reBF_simFF.levels_dBSPL(LEVind)),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
      xlim([0 PSTxmax_sec]*1000)
      YLIMunit=10^floor(log10(lowBF));
      YLIMS=floor(lowBF/YLIMunit)*YLIMunit*[1 10];
      ylim(YLIMS)  % Show 1 decade of BFs
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
      hold off
      
      
      %%%% Rate Plot
      PLOTnum=(ROWind-1)*NUMcols+2;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Rate (sp/sec)')
      %       ylabel('Best Frequency (kHz)')
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
      xlim([0 300])
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
      set(PLOThand,'XDir','reverse')
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
      hold off
      
      %%%% Synch Plot
      PLOTnum=(ROWind-1)*NUMcols+3;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Synch Coef (to f)')
      %       ylabel('Best Frequency (kHz)')
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
      xlim([0 1])
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
      set(PLOThand,'XDir','reverse')
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
      hold off
      
      %%%% Phase Plot
      PLOTnum=(ROWind-1)*NUMcols+4;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Phase (cycles of f)')
      %       ylabel('Best Frequency (kHz)')
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
      xlim([-pi pi])
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
      set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
      ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
      hold off
      
   end   
   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Plots
   if isfield(unit,'EH_reBF_simFF')
      for FeatIND=FeatINDs
         ROWind=ROWind+1;
         for HarmonicsIND=1:1
            for PolarityIND=1:1
               eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
               if ~isempty(yTEMP)
                  %%%% EH_reBF plots
                  LEVind=find(yTEMP.levels_dBSPL==LEVEL);
                  PERhists{ROWind,LEVind}=cell(size(yTEMP.BFs_kHz));
                  PERhistXs_sec{ROWind,LEVind}=cell(size(yTEMP.BFs_kHz));
                  PSTxmax_sec=-999;
                  for BFind=1:length(yTEMP.BFs_kHz)
                     if ~isempty(yTEMP.picNums{LEVind,BFind})
                        PIC=concatPICS_NOHR(yTEMP.picNums{LEVind,BFind},yTEMP.excludeLines{LEVind,BFind});
                        % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                        PIC=simFF_PICshift(PIC);
                        PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
                        
                        BFs_kHz{ROWind,LEVind}(BFind)=yTEMP.BFs_kHz(BFind);
                        PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
                        PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;
                        if max(PERhistXs_sec{ROWind,LEVind}{BFind})>PSTxmax_sec
                           PSTxmax_sec=max(PERhistXs_sec{ROWind,LEVind}{BFind});
                        end                    
                        Nsps{ROWind,LEVind}(BFind)=PIC.PERhist.NumDrivenSpikes;
                        Rates{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.SynchRate_PERhist(1);
                        if PIC.SynchRate_PERhist.FeatureRaySig(FeatIND)
                           Synchs{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeatureSynchs(FeatIND);
                           Phases{ROWind,LEVind}(BFind)=PIC.SynchRate_PERhist.FeaturePhases(FeatIND);
                        else
                           Synchs{ROWind,LEVind}(BFind)=NaN;
                           Phases{ROWind,LEVind}(BFind)=NaN;
                        end
                     else
                        BFs_kHz{ROWind,LEVind}(BFind)=yTEMP.BFs_kHz(BFind);
                        Nsps{ROWind,LEVind}(BFind)=NaN;
                        Rates{ROWind,LEVind}(BFind)=NaN;
                        Synchs{ROWind,LEVind}(BFind)=NaN;
                        Phases{ROWind,LEVind}(BFind)=NaN;
                     end
                     
                     
                  end
                  
                  %%%% Spatio-Temporal Plots
                  %                      figure(round(yTEMP.levels_dBSPL(LEVind)))
                  PLOTnum=(ROWind-1)*NUMcols+1;
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  for BFind=1:length(BFs_kHz{ROWind,LEVind})
                     if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                        LINEwidth=2;
                     else
                        LINEwidth=.5;
                     end
                     % ?? Is it OK to plot this as log, or do we need to make it linear?  - should be OK visually
                     if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000,trifilt(PERhists{ROWind,LEVind}{BFind}/NORMrate,TFiltWidth)+BFs_kHz{ROWind,LEVind}(BFind),'LineWidth',LINEwidth)
                        hold on
                     end
                  end
                  xlabel('Time (ms)')
                  ylabel('Effective Best Frequency (kHz)')
                  if ROWind==2
                     title(sprintf('     %s @ %.f dB SPL,     Harm: %d, Polarity: %d',FeaturesText{FeatIND},yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  else
                     title(sprintf('%s',FeaturesText{FeatIND}),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  end
                  xlim([0 PSTxmax_sec]*1000)
                  YLIMunit=10^floor(log10(lowBF));
                  YLIMS=floor(lowBF/YLIMunit)*YLIMunit*[1 10];
                  ylim(YLIMS)  % Show 1 decade of BFs
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
                  hold off
                  
                  
                  %%%% Rate Plot
                  PLOTnum=(ROWind-1)*NUMcols+2;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel('Rate (sp/sec)')
                  %       ylabel('Best Frequency (kHz)')
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
                  xlim([0 300])
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
                  set(PLOThand,'XDir','reverse')
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
                  hold off
                  
                  %%%% Synch Plot
                  PLOTnum=(ROWind-1)*NUMcols+3;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
                  xlim([0 1])
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
                  set(PLOThand,'XDir','reverse')
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
                  hold off
                  
                  %%%% Phase Plot
                  PLOTnum=(ROWind-1)*NUMcols+4;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
                  %       ylabel('Best Frequency (kHz)')
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Same as Spat-Temp plots
                  xlim([-pi pi])
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YLIMS(1):YLIMunit:YLIMS(end),'YTickLabel',YLIMS(1):YLIMunit:YLIMS(end))
                  set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
                  ylim([YLIMS(1) YLIMS(1)*2^YLIMocts])  % Show 2 octaves of BFs, with semilog labels
                  hold off
                  
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end
   
   
   Xcorner=0.05;
   Xwidth1=.5;
   Xshift1=0.05;
   Xwidth2=.1;
   Xshift2=0.03;
   
   Ycorner=0.05;
   Yshift=0.07;
   Ywidth=(1-NUMrows*(Yshift+.01))/NUMrows;   %.26 for 3; .42 for 2
   
   TICKlength=0.02;
   
   if NUMrows>2
      set(h9,'Position',[Xcorner Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h10,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h11,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h12,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-3)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>1
      set(h5,'Position',[Xcorner Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h6,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h7,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h8,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-2)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   set(h1,'Position',[Xcorner Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
   set(h2,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h3,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   set(h4,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-1)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])

   orient landscape
end


return;

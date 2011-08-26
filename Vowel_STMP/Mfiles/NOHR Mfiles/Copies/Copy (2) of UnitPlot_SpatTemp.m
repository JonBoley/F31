function UnitPlot_SpatTemp(ExpDate,UnitName)
% File: UnitPlot_SpatTemp.m
%
% Modified Date: 07Jan2005 (M Heinz)
% Modified From:  UnitPlot_SpatTemp.m and UnitPlot_SpatTemp.m to use either Interleaved Data or regular in 1 file
% 
% Date: 01Nov2004 (M. Heinz) (Modified from UnitPlot_EHrBF_simFF.m)
% For: NOHR Experiments
%
% ExpDate: e.g., '080204' (converted later)
% UnitName: '1.29' (converted later)
%
% Plots Simulated spatio-temporal response pattern (ala Shamma 1985) from EH_reBF_simFF and/or T_reBF_simFF data 
% for a given experiment and unit.  Loads 'UNITSdata/unit.T.U.mat' file.
% UnitAnal_EHrBF_simFF.m and/or UnitAnal_TonerBF_simFF.m performs the relevant analysis.
%
% 11/2/04 TO DO
% 1) only plots Fn/Tn{1,1} for now, need to add rest later
% *2) Separate Calcs and Plots, to allow smart plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global NOHR_dir NOHR_ExpList
global SavedPICS SavedPICnums SavedPICSuse
global FeaturesText FormsAtHarmonicsText InvertPolarityText
SavedPICSuse=1; SavedPICS=[]; SavedPICnums=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Verify in passed parameters if needed
if ~exist('ExpDate','var')
   ExpDate=0;
   %%% HARD CODE FOR NOW
   ExpDate='080204'
   % ExpDate='111804'
   while ~ischar(ExpDate)
      ExpDate=input('Enter Experiment Date (e.g., ''080204''): ');
   end
end
if ~exist('UnitName','var')
   UnitName=0;
   
   %%% HARD CODE FOR NOW
   if strcmp(ExpDate,'080204')
      UnitName='1.29'
   elseif strcmp(ExpDate,'111804')
      UnitName='1.28'
   end
   
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
disp(sprintf('Plotting Spatio-Temporal Patterns for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify UNITSdata/unit analyses are all done ahead of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Load unit structure for this unit
UnitFileName=sprintf('unit.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile('UNITSdata',UnitFileName) ''');'])
% If UNITSdata file does not exist, load DataList to see if interleaved or not, and then run necessary analyses
if isempty(ddd)
   FileName=strcat('DataList_',ExpDateText,'.mat');
   disp(['   *** Loading file: "' FileName '" because UNITSdata/unit does not exist yet'])
   eval(['load ' FileName])
   
   if isfield(DataList.Units{TrackNum,UnitNum},'Tone_reBF')   
      UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   end
   if isfield(DataList.Units{TrackNum,UnitNum},'EH_reBF')   
      UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);
   end
   % Estimate SR
   UnitCalc_SR(ExpDate,UnitName);   
   
end
eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])

% Make sure there is either Tone_reBF or EH_reBF data
if (~isfield(unit,'EH_reBF'))&(~isfield(unit,'Tone_reBF'))
   disp(sprintf('*** No Tone_reBF or EH_reBF data for this unit!!!'))
   beep
   return;
end

% If EH_reBF_simFF analysis is not completed, run here
if isfield(unit,'EH_reBF')&~isfield(unit,'EH_reBF_simFF')
   UnitAnal_EHrBF_simFF(ExpDate,UnitName,0);   
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

% If Tone_reBF_simFF analysis is not completed, run here
if isfield(unit,'Tone_reBF')&~isfield(unit,'Tone_reBF_simFF')
   UnitAnal_TrBF_simFF(ExpDate,UnitName,0);
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

% Make sure SR is estimated
if isempty(unit.Info.SR_sps)
   UnitCalc_SR(ExpDate,UnitName);   
   eval(['load ''' fullfile('UNITSdata',UnitFileName) ''''])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this unit, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
   
   F0min=Inf;
   for FeatIND=FeatINDs
      for HarmonicsIND=1:2
         for PolarityIND=1:2
            eval(['yTEMP=unit.EH_reBF_simFF.' FeaturesText{FeatIND} '{HarmonicsIND,PolarityIND};'])
            if ~isempty(yTEMP)
               levels_dBSPL=union(levels_dBSPL,yTEMP.levels_dBSPL);
               BFsTEMP_kHz=union(BFsTEMP_kHz,yTEMP.BFs_kHz);
               % Find minimum F0 for PERhist XMAX
               for i=1:length(yTEMP.FeatureFreqs_Hz)
                  if yTEMP.FeatureFreqs_Hz{i}(1)<F0min
                     F0min=yTEMP.FeatureFreqs_Hz{i}(1);
                  end
               end
            end
         end
      end
   end
end
lowBF=min(BFsTEMP_kHz);
highBF=max(BFsTEMP_kHz);
clear BFsTEMP_kHz;
TFiltWidth=1;   % What is a good number here?? is 1 OK, ow, you get major smoothing for low F0, and not much smoothing for high F0s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL CALCS (BEFORE PLOTTING), e.g., PERhists, SCCs, ...
%%%%   - runs through all BFs, levels and saves ALL calcs prior to PLOTTING (allows amart plotting based on ALL data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BFs_kHz=cell(NUMrows,length(levels_dBSPL));
PERhists=cell(NUMrows,length(levels_dBSPL));
PERhists_Smoothed=cell(NUMrows,length(levels_dBSPL));
PERhistXs_sec=cell(NUMrows,length(levels_dBSPL));  % for plotting
Nsps=cell(NUMrows,length(levels_dBSPL));
Rates=cell(NUMrows,length(levels_dBSPL));
Synchs=cell(NUMrows,length(levels_dBSPL));
Phases=cell(NUMrows,length(levels_dBSPL));
PERhistsMAX=0;
PERhistsYCHANS=0;

for LEVEL=levels_dBSPL
% beep
% disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
% for LEVEL=levels_dBSPL(end)
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Calcs
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      PERhists{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      
      PERhistXs_sec{ROWind,LEVind}=cell(size(unit.Tone_reBF_simFF.BFs_kHz));
      for BFind=1:length(unit.Tone_reBF_simFF.BFs_kHz)
         if ~isempty(unit.Tone_reBF_simFF.picNums{LEVind,BFind})
            
            PIC=concatPICS_NOHR(unit.Tone_reBF_simFF.picNums{LEVind,BFind},unit.Tone_reBF_simFF.excludeLines{LEVind,BFind});
            % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
            PIC=simFF_PICshift(PIC);
            PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
            
            BFs_kHz{ROWind,LEVind}(BFind)=unit.Tone_reBF_simFF.BFs_kHz(BFind);
            PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
            
            % Filter PERhists with Triangular Filter: filter 3 reps, take middle rep to avoid edge effects and keep periodic
            N=length(PERhists{ROWind,LEVind}{BFind});
            SmoothedPERhist=trifilt(repmat(PERhists{ROWind,LEVind}{BFind},1,3),TFiltWidth);
            PERhists_Smoothed{ROWind,LEVind}{BFind}=SmoothedPERhist(N+1:2*N); % Take middle 3rd
            % Determine Maximum of ALL plotted PERhists (i.e., post-Smoothing)
            if max(PERhists_Smoothed{ROWind,LEVind}{BFind})>PERhistsMAX
               PERhistsMAX=max(PERhists_Smoothed{ROWind,LEVind}{BFind});
            end
            
            PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;
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
      % Determine how many CHANNELS to plot
      if length(unit.Tone_reBF_simFF.BFs_kHz)>PERhistsYCHANS
         PERhistsYCHANS=length(unit.Tone_reBF_simFF.BFs_kHz);
      end
   end   % if Tone data
   
   
   %%%%%%%%%%%%%%%%%%%% EH_reBF Calcs
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
                  for BFind=1:length(yTEMP.BFs_kHz)
                     if ~isempty(yTEMP.picNums{LEVind,BFind})
                        PIC=concatPICS_NOHR(yTEMP.picNums{LEVind,BFind},yTEMP.excludeLines{LEVind,BFind});
                        % Shift spikes and frequencies to simulate shifted-BF neuron with stimulus at nominal-BF
                        PIC=simFF_PICshift(PIC);
                        PIC=calcSynchRate_PERhist(PIC);  % Calculates PERIOD histogram as well
                        
                        BFs_kHz{ROWind,LEVind}(BFind)=yTEMP.BFs_kHz(BFind);
                        PERhists{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist;
                        
                        % Filter PERhists with Triangular Filter: filter 3 reps, take middle rep to avoid edge effects and keep periodic
                        N=length(PERhists{ROWind,LEVind}{BFind});
                        SmoothedPERhist=trifilt(repmat(PERhists{ROWind,LEVind}{BFind},1,3),TFiltWidth);
                        PERhists_Smoothed{ROWind,LEVind}{BFind}=SmoothedPERhist(N+1:2*N);
                        % Determine Maximum of ALL plotted PERhists (i.e., post-Smoothing)
                        if max(PERhists_Smoothed{ROWind,LEVind}{BFind})>PERhistsMAX
                           PERhistsMAX=max(PERhists_Smoothed{ROWind,LEVind}{BFind});
                        end
                        
                        PERhistXs_sec{ROWind,LEVind}{BFind}=PIC.PERhist.PERhist_X_sec;
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
                  % Determine how many CHANNELS to plot
                  if length(yTEMP.BFs_kHz)>PERhistsYCHANS
                     PERhistsYCHANS=length(yTEMP.BFs_kHz);
                  end
               end %End if data for this condition, plot
            end % End Polarity
         end % End Harmonics
      end % End Feature
   end % If EHrBF data
end % Levels


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DO ALL PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FeatureColors={'r','g'};
PERhist_XMAX=1/F0min*1000;
XLIMITS_perhist=[0 PERhist_XMAX];
XLIMITS_rate=[0 300];
XLIMITS_synch=[0 1];
XLIMITS_phase=[-pi pi];
%% Find PERhist_YMAX
PERhistGAIN=3.0; % # of channels covered by max PERhist
PERhists_logCHwidth=log10(highBF/lowBF)/(PERhistsYCHANS-1);  % log10 channel width
PERhist_YMIN=lowBF;
PERhist_YMAX=10^((PERhistsYCHANS-1+PERhistGAIN)*PERhists_logCHwidth)*lowBF;    % sets an extra (GAIN-1) log channel widths
YLIMITS=[PERhist_YMIN PERhist_YMAX];  % Used for all plots
%% This  is ALL needed to get the right LOG Yticks!!
YLIMunit=10^floor(log10(lowBF));
YLIMS=floor(lowBF/YLIMunit)*YLIMunit*[1 100]; % Do two decades to be sure we get all the ticks
YTICKS=[YLIMS(1):YLIMunit:YLIMS(1)*10 YLIMS(1)*20:YLIMunit*10:YLIMS(end)];
BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave

% Setup parameters for title
if isempty(unit.Info.Threshold_dBatten)
   unit.Info.Threshold_dBatten=NaN;
end
if isempty(unit.Info.SR_sps)
   unit.Info.SR_sps=NaN;
end
if isempty(unit.Info.Q10)
   unit.Info.Q10=NaN;
end

for LEVEL=levels_dBSPL
% beep
% disp('***** HARD CODED FOR ONLY 1 (highest) LEVEL *****')
% for LEVEL=levels_dBSPL(end)
   figure(round(LEVEL)); clf
   set(gcf,'pos',[420     4   977   976])
   ROWind=0;
   
   %%%%%%%%%%%%%%%%%%%% Tone Plots
   if isfield(unit,'Tone_reBF_simFF')
      ROWind=ROWind+1;
      
      %%%% Tone_reBF plots
      LEVind=find(unit.Tone_reBF_simFF.levels_dBSPL==LEVEL);
      
      %       figure(round(unit.Tone_reBF_simFF.levels_dBSPL(LEVind)))
      PLOTnum=(ROWind-1)*NUMcols+1;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      for BFind=1:length(BFs_kHz{ROWind,LEVind})   
         if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
            LINEwidth=2;
         else
            LINEwidth=.5;
         end
         NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;  % Normalizes so each plot is equal size on log y-axis
         semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000, ...
            trifilt(PERhists_Smoothed{ROWind,LEVind}{BFind}*NormFact,TFiltWidth)+BFs_kHz{ROWind,LEVind}(BFind), ...
            'LineWidth',LINEwidth)
         hold on
      end
      semilogy(XLIMITS_perhist,unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Time (ms)')
      ylabel('Effective Best Frequency (kHz)')
      title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB atten, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL', ...
         ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBatten,unit.Info.SR_sps,unit.Info.Q10,'TONE', ...
         unit.Tone_reBF_simFF.levels_dBSPL(LEVind)),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
      %       xlim([0 max(PERhistXs_sec{ROWind,LEVind}{BFind}*1000)]) % Different xlim for TONES ??? 
      xlim(XLIMITS_perhist)
      PLOThand=eval(['h' num2str(PLOTnum)]);
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      text(XLIMITS_perhist(2),YLIMITS(1),'1/f','units','data','HorizontalAlignment','center','VerticalAlignment','top')
      hold off
      
      %%%% Rate Plot
      PLOTnum=(ROWind-1)*NUMcols+2;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]'))
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_rate)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Synch Plot
      PLOTnum=(ROWind-1)*NUMcols+3;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Synch Coef (to f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_synch)
      set(PLOThand,'XDir','reverse')
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
      ylim(YLIMITS)  % Same Ylimits for all plots
      hold off
      
      %%%% Phase Plot
      PLOTnum=(ROWind-1)*NUMcols+4;   
      eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
      semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
      hold on
      semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
      xlabel('Phase (cycles of f)')
      PLOThand=eval(['h' num2str(PLOTnum)]);
      xlim(XLIMITS_phase)
      set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
      set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
      ylim(YLIMITS)  % Same Ylimits for all plots
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
                  
                  %%%% Spatio-Temporal Plots
                  PLOTnum=(ROWind-1)*NUMcols+1;
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  for BFind=1:length(BFs_kHz{ROWind,LEVind})
                     if ismember(BFind,find(abs(log2(BFs_kHz{ROWind,LEVind}/unit.Info.BF_kHz))<BFoctCRIT))
                        LINEwidth=2;
                     else
                        LINEwidth=.5;
                     end
                     % This normalization plots each signal the same size on a log scale
                     if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                        NormFact=(10^(PERhistGAIN*PERhists_logCHwidth)-1)*BFs_kHz{ROWind,LEVind}(BFind)/PERhistsMAX;
                        semilogy(PERhistXs_sec{ROWind,LEVind}{BFind}*1000, ...
                           trifilt(PERhists_Smoothed{ROWind,LEVind}{BFind}*NormFact,TFiltWidth)+BFs_kHz{ROWind,LEVind}(BFind), ...
                           'LineWidth',LINEwidth)
                        hold on
                     end
                  end
                  xlabel('Time (ms)')
                  ylabel('Effective Best Frequency (kHz)')                  
                  if ROWind==1
                     title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB atten, SR=%.1f sps, Q10=%.1f\n%s @ %.f dB SPL,   Harm: %d, Polarity: %d', ...
                        ExpDate,UnitName,unit.Info.BF_kHz,unit.Info.Threshold_dBatten,unit.Info.SR_sps,unit.Info.Q10,FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  else
                     title(sprintf('%s @ %.f dB SPL,   Harm: %d, Polarity: %d',FeaturesText{FeatIND}, ...
                        yTEMP.levels_dBSPL(LEVind),HarmonicsIND,PolarityIND),'units','norm','pos',[.1 1 0],'HorizontalAlignment','left')
                  end
                  xlim(XLIMITS_perhist)
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=1:length(FeaturesText)
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_perhist,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                        text(XLIMITS_perhist(2)*1.005,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000, ...
                           sprintf('%s (%.1f)',FeaturesText{FeatINDPlot},yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000), ...
                           'units','data','HorizontalAlignment','left','VerticalAlignment','middle','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                     for BFind=1:length(BFs_kHz{ROWind,LEVind})
                        if ~isempty(PERhistXs_sec{ROWind,LEVind}{BFind})
                           if (FeatINDPlot<=FeatIND)
                              if (FeatINDPlot>1)
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                              else
                                 text(1000/yTEMP.FeatureFreqs_Hz{1}(FeatINDPlot),YLIMITS(1),sprintf('1/%s',FeaturesText{FeatINDPlot}),'units','data', ...
                                    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',6,'Color','k')
                              end
                           end
                        end
                     end
                  end
                  hold off
                  
                  
                  %%%% Rate Plot
                  PLOTnum=(ROWind-1)*NUMcols+2;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Rates{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy(Nsps{ROWind,LEVind}/10,BFs_kHz{ROWind,LEVind},'m+','MarkerSize',4)
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Rate (sp/sec)\n[+: # of spikes/10]'))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_rate)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=1:length(FeaturesText)
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_rate,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Synch Plot
                  PLOTnum=(ROWind-1)*NUMcols+3;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Synchs{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Synch Coef (to %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_synch)
                  set(PLOThand,'XDir','reverse')
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=1:length(FeaturesText)
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_synch,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
                  hold off
                  
                  %%%% Phase Plot
                  PLOTnum=(ROWind-1)*NUMcols+4;   
                  eval(['h' num2str(PLOTnum) '=subplot(NUMrows,NUMcols,PLOTnum);'])
                  semilogy(Phases{ROWind,LEVind},BFs_kHz{ROWind,LEVind},'*-')
                  hold on
                  semilogy([-1000 1000],unit.Info.BF_kHz*[1 1],'k:')
                  xlabel(sprintf('Phase (cycles of %s)',FeaturesText{FeatIND}))
                  PLOThand=eval(['h' num2str(PLOTnum)]);
                  xlim(XLIMITS_phase)
                  set(PLOThand,'XDir','reverse','XTick',[-pi -pi/2 0 pi/2 pi],'XTickLabel',[-1 -1/2 0 1/2 1])
                  set(PLOThand,'YTick',YTICKS,'YTickLabel',YTICKS)
                  ylim(YLIMITS)  % Same Ylimits for all plots
                  %%%%%%%%%%%%%%%%%%%%%
                  % Plot lines at all features
                  for FeatINDPlot=1:length(FeaturesText)
                     if (yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000>=YLIMITS(1))&(yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000<=YLIMITS(2))
                        semilogy(XLIMITS_phase,yTEMP.FeatureFreqs_Hz{LEVind}(FeatINDPlot)/1000*[1 1],':','Color',FeatureColors{-rem(FeatINDPlot,2)+2})
                     end
                  end
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
   
   if NUMrows>4
      set(h17,'Position',[Xcorner Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h18,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h19,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h20,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-5)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
   if NUMrows>3
      set(h13,'Position',[Xcorner Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth1 Ywidth],'TickLength',[TICKlength 0.025])
      set(h14,'Position',[Xcorner+Xwidth1+Xshift1 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h15,'Position',[Xcorner+Xwidth1+Xshift1+Xwidth2+Xshift2 Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
      set(h16,'Position',[Xcorner+Xwidth1+Xshift1+2*(Xwidth2+Xshift2) Ycorner+(NUMrows-4)*(Ywidth+Yshift) Xwidth2 Ywidth],'TickLength',[TICKlength 0.025])
   end
   
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


% Turn off saved PICS feature
SavedPICS=[]; SavedPICnums=[];
SavedPICSuse=0;

return;

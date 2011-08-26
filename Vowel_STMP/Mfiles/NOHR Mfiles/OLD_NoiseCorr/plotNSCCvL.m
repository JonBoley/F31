% File: plotNSCvL.m
%
% For two given units and STIMind, calculates the NSCC (Normalized Shuffled Cross-Correlogram) as
% a function of level, as well as NSACs for both units
%
% 9/3/03: M. Heinz
% From plotNSCCvL.m
%


clear
global Recruit_dir

%%%%%%%%%%%%%%%%%%%%%%%%% Examples
ExpDates={'061202','041701','040501','013001','061202','012501','112901','041701','012501','112901','112901','112901'};
UnitNames={'2.08','7.10','3.08','2.04','4.01','2.07','6.03','7.05','2.05','2.04','2.01','3.02'};
EG_BFs=[
   1  .38;
   2  1.1;  %% LF1
   3 1.3;  %% LF2 (.40 oct diff from #2)
   4 2.42;
   5 2.7;
   6 3.4;
   7 4.3;
   8 4.4;  %% HF1
   9 5.4;
   10 6.2;
   11 6.5;  %% HF2 (.56 oct diff from #8)
   12 6.8];
EG_NUMS={10,12};


%%%%%%%%%%%%%%% Parameters
paramsIN.StartTime=.01;
paramsIN.EndTime=.21;
Duration=paramsIN.EndTime-paramsIN.StartTime;
DELAYbinwidth=50e-6;
STIMind=4;  % Noise

LRangedB=5;

%%%%%%%%%%%%%%%%%%%%%%%%% Load Units
disp(sprintf('\n'))
for UNITind=1:2
   ExpDate{UNITind}=ExpDates{EG_NUMS{UNITind}};
   UnitName{UNITind}=UnitNames{EG_NUMS{UNITind}};
   ExpName{UNITind}=strcat('Exp',ExpDate{UNITind});
   eval(['cd ''' fullfile(Recruit_dir,'NewData',ExpName{UNITind}) ''''])
   eval(['load UNITSdata/unit.' UnitName{UNITind} '.mat'])
   units{UNITind}=unit;
   clear unit
   disp(sprintf('Unit%d [%s, %s]: BF=%.2f kHz, Thr=%.1f dB SPL; SR=%.1f sp/sec, Q10=%.1f',UNITind,units{UNITind}.ExpDate,units{UNITind}.Unit,units{UNITind}.BF_kHz,units{UNITind}.Thresh_dBSPL,units{UNITind}.SR,units{UNITind}.Q10))
   plotAVGdata(units{UNITind});
   set(gcf,'Name',sprintf('Unit %d',UNITind'))
   if UNITind==1
      eval(['cd ''' fullfile(Recruit_dir,'Mfiles','MGHfunctions','NoiseCorr') ''''])
      saveas(gcf,'unit1FIG')
   else
      openfig('unit1FIG','reuse');
   end
   BF{UNITind}=units{UNITind}.BF_kHz*1000;
   
end
disp(sprintf('Octave Difference: %.2f',log2(BF{2}/BF{1})))

%%%%%%%%%% Get all spikes for this STIMind: All REPS, All LEVELS
CALCnew=1;
if CALCnew|~exist('EGspikes_NSCC.mat')
   for UNITind=1:2
      [RawSpikes{UNITind},levels{UNITind}]=collectSpikes(units{UNITind},STIMind,paramsIN);  %RawSpikes: NUMstimREPS x NUMlevels cell array; levels: NUMlevels vector
   end
   eval(['cd ''' fullfile(Recruit_dir,'Mfiles','MGHfunctions','NoiseCorr') ''''])
   save EGspikes_NSCC RawSpikes levels
else
   eval(['cd ''' fullfile(Recruit_dir,'Mfiles','MGHfunctions','NoiseCorr') ''''])
   load EGspikes_NSCC 
end
   
for UNITind=1:2
   NUMstimREPS{UNITind}=size(RawSpikes{UNITind},1);
end
%%% Write separate funct for AN model, then use same analyis for all!


%%%%%%%%%%%%%%%%%%%%%%%%% Setup Levels to plot
offset=0; 
for UNITind=1:2
   PLOTlevs{UNITind}=unique(round(levels{UNITind}/10)*10+offset);  % spectral levels in dB SPL
   PLOTlevs{UNITind}=PLOTlevs{UNITind}(find((PLOTlevs{UNITind}>=min(levels{UNITind}))&(PLOTlevs{UNITind}<=max(levels{UNITind}))));
   PLOTlevINDS{UNITind}=zeros(1,length(PLOTlevs{UNITind}));
   for LEVind=1:length(PLOTlevs{UNITind})
      PLOTlevINDS{UNITind}(LEVind)=find(round(levels{UNITind})==PLOTlevs{UNITind}(LEVind));
   end
   %%% Include LRangedB range to increase number of spikes
   NUMspikeREPS{UNITind}=LRangedB*NUMstimREPS{UNITind};
   dBoff=(LRangedB-1)/2;
   if(dBoff-round(dBoff))  error('Need odd LRangedB!!');  end
   % Take out any PLOTlevs for which we can't get all levels
   if PLOTlevINDS{UNITind}(1)-dBoff<1
      PLOTlevINDS{UNITind}=PLOTlevINDS{UNITind}(2:end);
      PLOTlevs{UNITind}=PLOTlevs{UNITind}(2:end);
   end
   if PLOTlevINDS{UNITind}(end)+dBoff>length(levels{UNITind})
      PLOTlevINDS{UNITind}=PLOTlevINDS{UNITind}(1:end-1);
      PLOTlevs{UNITind}=PLOTlevs{UNITind}(1:end-1);
   end
end
%% Take out non overlapping levels
PLOTlevsSAME=intersect(PLOTlevs{1},PLOTlevs{2});
for UNITind=1:2
   PLOTlevINDS{UNITind}=PLOTlevINDS{UNITind}(find(PLOTlevs{UNITind}==min(PLOTlevsSAME)):find(PLOTlevs{UNITind}==max(PLOTlevsSAME)));
end
PLOTlevs=PLOTlevsSAME; clear PLOTlevsSAME


%%%%%% Run through all levels
% LEVEL=10;
% PLOTINDS=find(round((PLOTlevs-LEVEL)/10)==0);
PLOTINDS=3:length(PLOTlevs);
for pLEVind=PLOTINDS

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% Create spike trains for Correlation Analysis
   SpikeTrains=cell(1,2);
   for UNITind=1:2
      SpikeTrains{UNITind}=cell(NUMspikeREPS{UNITind},1);
      spikeREPind=0;
      for LEVoffIND=(PLOTlevINDS{UNITind}(pLEVind)-dBoff):(PLOTlevINDS{UNITind}(pLEVind)+dBoff)
         for stimREPind=1:NUMstimREPS{UNITind}
            spikeREPind=spikeREPind+1;
            SpikeTrains{UNITind}{spikeREPind}=RawSpikes{UNITind}{stimREPind,LEVoffIND};
         end
      end
   end
   
   
   %%%%%%%%%%%%%%%% Compute Normalized Shuffled-AutoCorrelogram (SAC)
   Twin=[paramsIN.StartTime,paramsIN.EndTime];
   tic
   for UNITind=1:2
      [NSAC{UNITind}{pLEVind},NSACdelays{UNITind},AVGrate{UNITind}{pLEVind},TOTALspikes{UNITind}{pLEVind},ETIMEinner] = ShufAutoCorr(SpikeTrains{UNITind},DELAYbinwidth,Twin);
      disp(sprintf('NSAC %d: Level=%.f dB SPL;\tTOTAL spikes=%d;\telapsed time=%.2f sec (inner: %.2f sec)',UNITind,PLOTlevs(pLEVind),TOTALspikes{UNITind}{pLEVind},toc,ETIMEinner))
   end
   
   % Maybe eventually output NSAC also, but probably best to keep separate
   [NSCC{pLEVind},NSCCdelays,AVGrateNSCC{pLEVind},TOTALspikesNSCC{pLEVind},ETIMEinner] = ShufCrossCorr(SpikeTrains,DELAYbinwidth,Twin);
   disp(sprintf('NSCC: Level=%.f dB SPL;\tTOTAL spikes=[%d,%d];\telapsed time=%.2f sec (inner: %.2f sec)',PLOTlevs(pLEVind),TOTALspikes{1}{pLEVind},TOTALspikes{2}{pLEVind},toc,ETIMEinner))
   
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Plot results
%% Params for PLots
NUM_BFcycles=10;   % Ruggero showed decay within 3-7 cycles
XLIMIT=max([NUM_BFcycles/BF{1}/1e-6 NUM_BFcycles/BF{2}/1e-6]);
TFiltWidth=9;

if length(PLOTINDS)>1
   NUMcols=2;
   NUMrows=ceil(length(PLOTINDS)/2);
   subPLOTINDs=[1:2:length(PLOTINDS) 2:2:length(PLOTINDS)];
end

%%% Find YLIMITS
Ymax=0;  YRmax=0;
for UNITind=1:2
   for pLEVind=PLOTINDS
      plotNSAC=trifilt(NSAC{UNITind}{pLEVind},TFiltWidth);
      plotNSAC2=plotNSAC(find((NSACdelays{UNITind}>=-XLIMIT)&(NSACdelays{UNITind}<=XLIMIT)));
      if max(plotNSAC2)>Ymax
         Ymax=max(plotNSAC2);
      end
      if max(plotNSAC2*AVGrate{UNITind}{pLEVind})>YRmax
         YRmax=max(plotNSAC2*AVGrate{UNITind}{pLEVind});
      end
   end
end

for UNITind=1:2
   %%%%%%%%%%%% Plot SAC   
   figure(UNITind*10+11); clf
   set(gcf,'pos',[477     5   720   976])
   for pLEVind=PLOTINDS
      if length(PLOTINDS)>1
         NUMcols=2;
         NUMrows=ceil(length(PLOTINDS)/2);
         subPLOTind=subPLOTINDs(find(PLOTINDS==pLEVind));
         subplot(NUMrows,NUMcols,subPLOTind)
      else
         subPLOTind=-99;
      end
      
      plot(NSACdelays{UNITind},AVGrate{UNITind}{pLEVind}*trifilt(NSAC{UNITind}{pLEVind},TFiltWidth))
      hold on
      plot(0*ones(1,2),[0 YRmax],'r--')
      for CYCind=1:NUM_BFcycles
         plot(CYCind/BF{UNITind}/1e-6*ones(1,2),[0 YRmax],'k:')
         plot(-CYCind/BF{UNITind}/1e-6*ones(1,2),[0 YRmax],'k:')
      end
      plot([-XLIMIT XLIMIT],AVGrate{UNITind}{pLEVind}*ones(1,2),'r--')
      text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrate{UNITind}{pLEVind},TOTALspikes{UNITind}{pLEVind}), ...
         'units','norm','VerticalAlignment','top')
      
      if (subPLOTind==1)|(subPLOTind<0)
         title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
            ExpName{UNITind},UnitName{UNITind},NUMstimREPS{UNITind},units{UNITind}.BF_kHz,LRangedB,TFiltWidth))
      end
      if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
         xlabel('Delay (\mus)','Interpreter','tex')
         ylabel('Discharge Rate (sp/sec)')
      end
      hold off
      axis([-XLIMIT XLIMIT 0 YRmax])
      
      %   input('Enter')
   end
   
   
   %%%%%%%%%%%% Plot NSACs
   figure(UNITind*10+12); clf
   set(gcf,'pos',[677     5   720   976])
   for pLEVind=PLOTINDS
      if length(PLOTINDS)>1
         NUMcols=2;
         NUMrows=ceil(length(PLOTINDS)/2);
         subPLOTind=subPLOTINDs(find(PLOTINDS==pLEVind));
         subplot(NUMrows,NUMcols,subPLOTind)
      else
         subPLOTind=-99;
      end
      
      plot(NSACdelays{UNITind},trifilt(NSAC{UNITind}{pLEVind},TFiltWidth))
      hold on
      plot(0*ones(1,2),[0 Ymax],'r--')
      for CYCind=1:NUM_BFcycles
         plot(CYCind/BF{UNITind}/1e-6*ones(1,2),[0 Ymax],'k:')
         plot(-CYCind/BF{UNITind}/1e-6*ones(1,2),[0 Ymax],'k:')
      end
      plot([-XLIMIT XLIMIT],[1 1],'r--')
      text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrate{UNITind}{pLEVind},TOTALspikes{UNITind}{pLEVind}), ...
         'units','norm','VerticalAlignment','top')
      
      if (subPLOTind==1)|(subPLOTind<0)
         title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
            ExpName{UNITind},UnitName{UNITind},NUMstimREPS{UNITind},units{UNITind}.BF_kHz,LRangedB,TFiltWidth))
      end
      if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
         xlabel('Delay (\mus)','Interpreter','tex')
         ylabel('Normalized SAC')
      end
      hold off
      axis([-XLIMIT XLIMIT 0 Ymax])
   end  
end


%%%%%%%%%%%% Plot SCC   
figure(41); clf
set(gcf,'pos',[477     5   720   976])
for pLEVind=PLOTINDS
   if length(PLOTINDS)>1
      NUMcols=2;
      NUMrows=ceil(length(PLOTINDS)/2);
      subPLOTind=subPLOTINDs(find(PLOTINDS==pLEVind));
      subplot(NUMrows,NUMcols,subPLOTind)
   else
      subPLOTind=-99;
   end
   
   plot(NSCCdelays,AVGrateNSCC{pLEVind}{1}*trifilt(NSCC{pLEVind},TFiltWidth))
   hold on
   plot(0*ones(1,2),[0 YRmax],'r--')
   for CYCind=1:NUM_BFcycles
      plot(CYCind/BF{1}/1e-6*ones(1,2),[0 YRmax],'k:')
      plot(-CYCind/BF{1}/1e-6*ones(1,2),[0 YRmax],'k:')
      plot(CYCind/BF{2}/1e-6*ones(1,2),[0 YRmax],'g:')
      plot(-CYCind/BF{2}/1e-6*ones(1,2),[0 YRmax],'g:')
   end
   plot([-XLIMIT XLIMIT],AVGrateNSCC{pLEVind}{1}*ones(1,2),'r--')
   text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrateNSCC{pLEVind}{UNITind},TOTALspikesNSCC{pLEVind}{UNITind}), ...
      'units','norm','VerticalAlignment','top')
   
   if (subPLOTind==1)|(subPLOTind<0)
      title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
         ExpName{UNITind},UnitName{UNITind},NUMstimREPS{UNITind},units{UNITind}.BF_kHz,LRangedB,TFiltWidth))
   end
   if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('Discharge Rate (sp/sec)')
   end
   hold off
   axis([-XLIMIT XLIMIT 0 YRmax])
   
   %   input('Enter')
end

   

%%%%%%%%%%%% Plot NSCC   
figure(42); clf
set(gcf,'pos',[677     5   720   976])
for pLEVind=PLOTINDS
   if length(PLOTINDS)>1
      NUMcols=2;
      NUMrows=ceil(length(PLOTINDS)/2);
      subPLOTind=subPLOTINDs(find(PLOTINDS==pLEVind));
      subplot(NUMrows,NUMcols,subPLOTind)
   else
      subPLOTind=-99;
   end
   
   plot(NSCCdelays,trifilt(NSCC{pLEVind},TFiltWidth))
   hold on
   plot(0*ones(1,2),[0 Ymax],'r--')
   for CYCind=1:NUM_BFcycles
      plot(CYCind/BF{1}/1e-6*ones(1,2),[0 Ymax],'k:')
      plot(-CYCind/BF{1}/1e-6*ones(1,2),[0 Ymax],'k:')
      plot(CYCind/BF{2}/1e-6*ones(1,2),[0 Ymax],'g:')
      plot(-CYCind/BF{2}/1e-6*ones(1,2),[0 Ymax],'g:')
   end
   plot([-XLIMIT XLIMIT],ones(1,2),'r--')
   text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrateNSCC{pLEVind}{UNITind},TOTALspikesNSCC{pLEVind}{UNITind}), ...
      'units','norm','VerticalAlignment','top')
   
   if (subPLOTind==1)|(subPLOTind<0)
      title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
         ExpName{UNITind},UnitName{UNITind},NUMstimREPS{UNITind},units{UNITind}.BF_kHz,LRangedB,TFiltWidth))
   end
   if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('Discharge Rate (sp/sec)')
   end
   hold off
   axis([-XLIMIT XLIMIT 0 Ymax])
   
   %   input('Enter')
end

   


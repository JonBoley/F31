% File: plotNSACvL.m
%
% For a given unit and STIMind calculates the NSAC (Normalized Shuffled Auto-Correlogram) as
% a function of level
%
% 8/20/03: M. Heinz
% From plotNSCCvL.m
%

%%%%%%%%%% TO DO
% *0) Setup general NSAC function
% *1) Speed up SAC calculation
%    *LATER: Try Ccode into MEX to check speed
% *2) Fix binning to be consistent
% *3) Cleanup basic/general multi-level layout for NSAC
% *3.5) Need to verify AVGrate CALC when not all 5 REPS are used!!!!, probably will underestimate rate!!! (e.g 041701, 7.10 L=-60)
%
% 08Nov2003
% **3.a) Get SAC-MEX going
% 3.b) Look at some basic trends across Level, to verify calc is OK
%          Clean up general NSAC function, e.g., in terms of parameters passed
% 
% 4) Setup general NSCC function, SCC-MEX
% 5) Setup basic/general multi-level layout for NSAC/NSCC (both NSACs and NSCC)
% 6) Look at examples and basic trends, see if obvious, ow/ test with AN model
%
% 7) Clean up analysis, find decent params for analysis (LRangedB, binwidth)
%    make sure transients aren't getting in there, i.e., start late enough!!
% 8) Look at Spectra of SAC!!  SCCs!!


clear
global Recruit_dir

%%%%%%%%%%%%%%%%%%%%%%%%% Examples
ExpDates={'061202','041701','040501','013001','061202','012501','112901','041701','012501','112901','112901','112901'};
UnitNames={'2.08','7.10','3.08','2.04','4.01','2.07','6.03','7.05','2.05','2.04','2.01','3.02'};
EG_BFs=[
   1  .38;
   2  1.1;  %% LF1
   3 1.3;  %% LF2 (.25 oct diff from #2)
   4 2.42;
   5 2.7;
   6 3.4;
   7 4.3;
   8 4.4;  %% HF1
   9 5.4;
   10 6.2;
   11 6.5;  %% HF2 (.56 oct diff from #8)
   12 6.8];
EG_NUM=2;


%%%%%%%%%%%%%%% Parameters
paramsIN.StartTime=.01;
paramsIN.EndTime=.21;
Duration=paramsIN.EndTime-paramsIN.StartTime;

DELAYbinwidth=50e-6;
STIMind=4;  % Noise


%%%%%%%%%%%%%%%%%%%%%%%%% Load Unit
ExpDate=ExpDates{EG_NUM};
UnitName=UnitNames{EG_NUM};
ExpName=strcat('Exp',ExpDate);
eval(['cd ''' fullfile(Recruit_dir,'NewData',ExpName) ''''])
eval(['load UNITSdata/unit.' UnitName '.mat'])
disp(sprintf('Unit [%s, %s]: BF=%.2f kHz, Thr=%.1f dB SPL; SR=%.1f sp/sec, Q10=%.1f',unit.ExpDate,unit.Unit,unit.BF_kHz,unit.Thresh_dBSPL,unit.SR,unit.Q10))
plotAVGdata(unit);
BF=unit.BF_kHz*1000;


%%%%%%%%%% Get all spikes for this STIMind: All REPS, All LEVELS
CALCnew=0;
if CALCnew|~exist('EGspikes_NSAC.mat')
   [RawSpikes,levels]=collectSpikes(unit,STIMind,paramsIN);  %RawSpikes: NUMstimREPS x NUMlevels cell array; levels: NUMlevels vector
   eval(['cd ''' fullfile(Recruit_dir,'Mfiles','MGHfunctions','NoiseCorr') ''''])
   save EGspikes_NSAC RawSpikes levels
else
   eval(['cd ''' fullfile(Recruit_dir,'Mfiles','MGHfunctions','NoiseCorr') ''''])
   load EGspikes_NSAC 
end
   
NUMstimREPS=size(RawSpikes,1);
%%% Write separate funct for AN model, then use same analyis for all!


%%%%%%%%%%%%%%%%%%%%%%%%% Setup Levels to plot
offset=0; 
PLOTlevs=unique(round(levels/10)*10+offset);  % spectral levels in dB SPL
PLOTlevs=PLOTlevs(find((PLOTlevs>=min(levels))&(PLOTlevs<=max(levels))));
PLOTlevINDS=zeros(1,length(PLOTlevs));
for LEVind=1:length(PLOTlevs)
   PLOTlevINDS(LEVind)=find(round(levels)==PLOTlevs(LEVind));
end

%%% Include LRangedB range to increase number of spikes
LRangedB=5;
NUMspikeREPS=LRangedB*NUMstimREPS;
dBoff=(LRangedB-1)/2;
if(dBoff-round(dBoff))  error('Need odd LRangedB!!');  end
% Take out any PLOTlevs for which we can't get all levels
if PLOTlevINDS(1)-dBoff<1
   PLOTlevINDS=PLOTlevINDS(2:end);
   PLOTlevs=PLOTlevs(2:end);
end
if PLOTlevINDS(end)+dBoff>length(levels)
   PLOTlevINDS=PLOTlevINDS(1:end-1);
   PLOTlevs=PLOTlevs(1:end-1);
end

%%%%%% Run through all levels
% LEVEL=10;
% %LEVEL=-50;
% PLOTINDS=find(round((PLOTlevs-LEVEL)/10)==0);
PLOTINDS=1:length(PLOTlevs);
for pLEVind=PLOTINDS

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%% Create spike trains for Correlation Analysis
   
   spikeREPind=0;
   for LEVoffIND=(PLOTlevINDS(pLEVind)-dBoff):(PLOTlevINDS(pLEVind)+dBoff)
      for stimREPind=1:NUMstimREPS
         if ~isempty(RawSpikes{stimREPind,LEVoffIND})  % Take out all empty spike trains
            spikeREPind=spikeREPind+1;
            SpikeTrain1{spikeREPind,1}=RawSpikes{stimREPind,LEVoffIND};
         end
      end
   end

   %%%%%%%%%%%%%%%% Compute Normalized Shuffled-AutoCorrelogram (SAC)
   Twin=[paramsIN.StartTime,paramsIN.EndTime];
   tic
   [NSAC{pLEVind},NSACdelays,AVGrate{pLEVind},TOTALspikes{pLEVind},ETIMEinner] = ShufAutoCorr(SpikeTrain1,DELAYbinwidth,Twin);
   disp(sprintf('Level=%.f dB SPL;\tTOTAL spikes=%d;\telapsed time=%.2f sec (inner: %.2f sec)',PLOTlevs(pLEVind),TOTALspikes{pLEVind},toc,ETIMEinner))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Plot results
%% Params for PLots
NUM_BFcycles=10;   % Ruggero showed decay within 3-7 cycles
XLIMIT=NUM_BFcycles/BF/1e-6;
TFiltWidth=9;

if length(PLOTINDS)>1
   NUMcols=2;
   NUMrows=ceil(length(PLOTINDS)/2);
   subPLOTINDs=[1:2:length(PLOTINDS) 2:2:length(PLOTINDS)];
end


%%% Find YLIMITS
Ymax=0;  YRmax=0;
for pLEVind=PLOTINDS
   plotNSAC=trifilt(NSAC{pLEVind},TFiltWidth);
   plotNSAC2=plotNSAC(find((NSACdelays>=-XLIMIT)&(NSACdelays<=XLIMIT)));
   if max(plotNSAC2)>Ymax
      Ymax=max(plotNSAC2);
   end
   if max(plotNSAC2*AVGrate{pLEVind})>YRmax
      YRmax=max(plotNSAC2*AVGrate{pLEVind});
   end
end

%%%%%%%%%%%% Plot SAC

figure(11); clf
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
   
   plot(NSACdelays,AVGrate{pLEVind}*trifilt(NSAC{pLEVind},TFiltWidth))
   hold on
   plot(0*ones(1,2),[0 YRmax],'r--')
   for CYCind=1:NUM_BFcycles
      plot(CYCind/BF/1e-6*ones(1,2),[0 YRmax],'k:')
      plot(-CYCind/BF/1e-6*ones(1,2),[0 YRmax],'k:')
   end
   plot([-XLIMIT XLIMIT],AVGrate{pLEVind}*ones(1,2),'r--')
   text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrate{pLEVind},TOTALspikes{pLEVind}), ...
      'units','norm','VerticalAlignment','top')
   
   if (subPLOTind==1)|(subPLOTind<0)
      title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
         ExpName,UnitName,NUMstimREPS,unit.BF_kHz, LRangedB,TFiltWidth))
   end
   if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('Discharge Rate (sp/sec)')
   end
   hold off
   axis([-XLIMIT XLIMIT 0 YRmax])
   
   %   input('Enter')
   
end


%%%%%%%%%%%% Plot NSAC
figure(12); clf
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
      
   plot(NSACdelays,trifilt(NSAC{pLEVind},TFiltWidth))
   hold on
   plot(0*ones(1,2),[0 Ymax],'r--')
   for CYCind=1:NUM_BFcycles
      plot(CYCind/BF/1e-6*ones(1,2),[0 Ymax],'k:')
      plot(-CYCind/BF/1e-6*ones(1,2),[0 Ymax],'k:')
   end
   plot([-XLIMIT XLIMIT],[1 1],'r--')
   text(.05,.975,sprintf('L=%d dB SPL\nRate = %.1f sp/s\nNsps=%d',PLOTlevs(pLEVind),AVGrate{pLEVind},TOTALspikes{pLEVind}), ...
      'units','norm','VerticalAlignment','top')
   
   if (subPLOTind==1)|(subPLOTind<0)
      title(sprintf('%s, %s (%d reps): BF= %.2f kHz\n(%d dB range; TFilt=%d)', ...
         ExpName,UnitName,NUMstimREPS,unit.BF_kHz, LRangedB,TFiltWidth))
   end
   if (subPLOTind==length(PLOTINDS)-1)|(subPLOTind<0)
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('Normalized SAC')
   end
   hold off
   axis([-XLIMIT XLIMIT 0 Ymax])
end





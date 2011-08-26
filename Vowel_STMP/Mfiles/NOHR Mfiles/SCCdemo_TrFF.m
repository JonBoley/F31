function [CD_us,NSCC_0,NSCC_CD,NSCC_HLCD,NSCC_max]=SCCdemo_TrFF(picNums1,picNums2,excludeLines1,excludeLines2,HLCD_us,SCCtitles,calcSAC)
% M.Heinz 25Sep2004.  Modified from PSTview.m
% Computes and plots SAC for each unit, then computes and plots SCC and returns XBF-corr rate
%
% Allows multiple pst-style pics to be concatenated.  Also allows lines to be excluded for each picture.
% Assumes current directory is the data directory
% All calcs are done from external MFile calls, while all plotting is done here
%
% Usage: PSTview(picNums,excludeLines)
% picNums: vector of picture numbers
% excludeLines: cell array with vectors for each picture with any lines to be excluded
% SCCtitles: cell array with two titles for 2 columns (2 units) and overall title

disp(sprintf('Unit1: PICS: %s\nUnit2: PICS: %s\n',mat2str(picNums1),mat2str(picNums2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get Spike Trains and Parameters organized
picNums{1}=picNums1; picNums{2}=picNums2; 
clear picNums1 picNums2
for i=1:2
   NUMpics{i}=length(picNums{i});
end
if ~exist('excludeLines1','var')
   excludeLines{1}=cell(1,NUMpics{1});
else
   excludeLines{1}=excludeLines1;
   clear excludeLines1
end
if ~exist('excludeLines2','var')
   excludeLines{2}=cell(1,NUMpics{2});
else
   excludeLines{2}=excludeLines2;
   clear excludeLines2
end
if ~exist('SCCtitles','var')
   SCCtitles={'UNIT 1','UNIT 2','SCCdemo_TONEreFF'};
end
if ~exist('calcSAC','var')
   calcSAC=1;
end

%%%% Concatenate pictures for each unit
PIC=cell(1,2);
for i=1:2
   PIC{i}=concatPICS_NOHR(picNums{i},excludeLines{i});
end

%%%% Set parameters
clear paramsIN
paramsIN.SAC.StartTime_sec=.01;
if (PIC{1}.x.Hardware.Trigger.StmOn/1000==PIC{2}.x.Hardware.Trigger.StmOn/1000)
   paramsIN.SAC.EndTime_sec=PIC{1}.x.Hardware.Trigger.StmOn/1000;
else
   error(sprintf('Stimulus Durations do not match: unit 1 duration = %.3f ms; unit 2 duration = %.3f ms',PIC{1}.x.Hardware.Trigger.StmOn/1000,PIC{2}.x.Hardware.Trigger.StmOn/1000))
end
paramsIN.SAC.DELAYbinwidth_sec=50e-6;


%%%% Create cell arrays of Spike Trains
SpikeTrains=cell(1,2);
Twin_sec=[paramsIN.SAC.StartTime_sec paramsIN.SAC.EndTime_sec];
for i=1:2
   SpikeTrains{i}=getDrivenSpikeTrains(PIC{i}.x.spikes{1},[],Twin_sec);
end


%% 9/27/04 TODO
% *1) Verify SAC calculation, by testing vs slower Mfiles with few reps
% *2) Clean up and take out all timing stuff
% *3) Clean up SACdemo plots for Units 1 and 2
% *4) Setup SCC calculation with MEX files,
% 5) Setup all SCC analysis into UnitsPlot



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute SAC
Duration_sec=diff(Twin_sec);
if calcSAC
   NSAC=cell(1,2); NSACdelays=cell(1,2); AVGrate=cell(1,2); TOTALspikes=cell(1,2);
   for i=1:2
      disp(sprintf('Computing SAC for Unit %d ........',i))
      [NSAC{i},NSACdelays{i},AVGrate{i},TOTALspikes{i}] = ShufAutoCorr(SpikeTrains{i},paramsIN.SAC.DELAYbinwidth_sec,Duration_sec);
      disp(sprintf('TOTAL spikes=%d',TOTALspikes{i}))
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compute SCC
% NSAC=cell(1,2); NSACdelays=cell(1,2); AVGrate=cell(1,2); TOTALspikes=cell(1,2); ETIMEinner=cell(1,2);
disp(sprintf('Computing SCC  ........',i))
[NSCC,NSCCdelays,AVGratesNSCC,TOTALspikesNSCC] = ShufCrossCorr(SpikeTrains,paramsIN.SAC.DELAYbinwidth_sec,Duration_sec);
disp(sprintf('NSCC: TOTAL spikes=[%d,%d]',TOTALspikesNSCC{1},TOTALspikesNSCC{2}))


%%%%%%%%% TODO
% *1) Get Spike Trains orgainized
% *2) go through old SAC stuff and bring over and get going with SAC - 
% *3) plot SAC for each unit
% *4) setup SCC
% *5) plot SCC

%%%%%%%%%%% TODO 9/27/04
% *1) InnerSCCmex seems ALL SET
% *2) Clean up testing, timing
% *3) Plot NSCC (see plotNSCCvL)
% 4) Look at real data, and see what we've got!!!
% 5) on to EHreFF, 1 good example to show effect
% 6) on to EHrBF, simulate 2 neurons, to show effect



NSCC_0=NSCC(find(NSCCdelays==0));

%%%% find CD, i.e., local max closest to zero delay
diff_LtoR=diff(NSCC);
diff_RtoL=diff(fliplr(NSCC));
diffDelays_LtoR=NSCCdelays(1:end-1);
diffDelays_RtoL=fliplr(NSCCdelays(2:end));
LocalMaxDelays=intersect(diffDelays_LtoR(find(diff_LtoR<0)),diffDelays_RtoL(find(diff_RtoL<0)));
[y,i]=min(abs(LocalMaxDelays));
CD_us=LocalMaxDelays(i);

NSCC_CD=NSCC(find(NSCCdelays==CD_us));
NSCC_HLCD=NSCC(find(NSCCdelays==HLCD_us));
NSCC_max=max(NSCC);

% SCCrate_zeroDelay=sqrt(AVGratesNSCC{1}*AVGratesNSCC{2})*NSCC_zeroDelay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SAC/SCC PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PIC{1}.FixedFreq_Hz==PIC{2}.FixedFreq_Hz)
   FF_Hz=PIC{1}.FixedFreq_Hz;
else
   error(sprintf('Fixed Frequencies do not match: unit 1 FF = %.1f Hz; unit 2 FF = %.1f Hz',PIC{1}.FixedFreq_Hz,PIC{2}.FixedFreq_Hz))
end

if strcmp(SCCtitles{3},'SCCdemo_TONEreFF')
   SCCtitles{3}=strcat(SCCtitles{3},sprintf('     (FF = %.3f kHz)',FF_Hz/1000));
end

%% Params for Plots
NUM_cycles=10;   % Ruggero showed decay within 3-7 cycles
XLIMIT=NUM_cycles/FF_Hz/1e-6;
TFiltWidth=1;

%%%% Plot SAC/NSAC
%%%%%%%%%%%% Plot SAC
YRmax=500;
Ymax=YRmax/150;

if calcSAC
   figure(31); clf
   set(gcf,'pos',[208     5   1189   976])
   for i=1:2
      subplot(2,2,i)
      plot(NSACdelays{i},AVGrate{i}*trifilt(NSAC{i},TFiltWidth))
      hold on
      plot(0*ones(1,2),[0 YRmax],'r--')
      for CYCind=1:NUM_cycles
         plot(CYCind/FF_Hz/1e-6*ones(1,2),[0 YRmax],'k:')
         plot(-CYCind/FF_Hz/1e-6*ones(1,2),[0 YRmax],'k:')
      end
      plot([-XLIMIT XLIMIT],AVGrate{i}*ones(1,2),'r--')
      
      eval(['ht' num2str(i) '=title(SCCtitles{i});'])
      eval(['set(ht' num2str(i) ',''units'',''norm'',''Interpreter'',''none'')'])
      
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('SAC Discharge Rate (sp/sec)')
      
      hold off
      xlim([-XLIMIT XLIMIT])
      ylim([0 YRmax])
      
      
      %%%%%%%%%%%% Plot NSAC
      subplot(2,2,i+2)
      plot(NSACdelays{i},trifilt(NSAC{i},TFiltWidth))
      hold on
      plot(0*ones(1,2),[0 Ymax],'r--')
      for CYCind=1:NUM_cycles
         plot(CYCind/FF_Hz/1e-6*ones(1,2),[0 Ymax],'k:')
         plot(-CYCind/FF_Hz/1e-6*ones(1,2),[0 Ymax],'k:')
      end
      plot([-XLIMIT XLIMIT],[1 1],'r--')
      title(sprintf('%d reps;   TOTAL spikes = %d;   AVGrate = %.2f sp/sec',length(SpikeTrains{i}),TOTALspikes{i}, AVGrate{i}))   
      xlabel('Delay (\mus)','Interpreter','tex')
      ylabel('Normalized SAC')
      hold off
      xlim([-XLIMIT XLIMIT])
      ylim([0 Ymax])
   end
   subplot(221)
   text(1.175,1.12,SCCtitles{3},'Units','Norm','Interpreter','none','Horiz','center')
end

%%%% Plot SCC/NSCC
%%%%%%%%%%%% Plot SCC
figure(32); clf
set(gcf,'pos',[679     5   720   976])
subplot(211)
plot(NSCCdelays,sqrt(AVGratesNSCC{1}*AVGratesNSCC{2})*trifilt(NSCC,TFiltWidth))
hold on
plot(0*ones(1,2),[0 YRmax],'r--')
for CYCind=1:NUM_cycles
   plot(CYCind/FF_Hz/1e-6*ones(1,2),[0 YRmax],'k:')
   plot(-CYCind/FF_Hz/1e-6*ones(1,2),[0 YRmax],'k:')
   %    plot(CYCind/FF_Hz{2}/1e-6*ones(1,2),[0 YRmax],'g:')
   %    plot(-CYCind/FF_Hz{2}/1e-6*ones(1,2),[0 YRmax],'g:')
end
plot([-XLIMIT XLIMIT],sqrt(AVGratesNSCC{1}*AVGratesNSCC{2})*ones(1,2),'r--')

eval(['ht' num2str(i) '=title(sprintf(''%s\n%s %20s %s'',SCCtitles{3},SCCtitles{1},'' '',SCCtitles{2}));'])
eval(['set(ht' num2str(i) ',''units'',''norm'',''Interpreter'',''none'')'])

xlabel('Delay (\mus)','Interpreter','tex')
ylabel('SCC Discharge Rate (sp/sec)')
hold off
xlim([-XLIMIT XLIMIT])
ylim([0 YRmax])


%%%%%%%%%%%% Plot NSCC   
subplot(212)
plot(NSCCdelays,trifilt(NSCC,TFiltWidth))
hold on

plot(CD_us,NSCC_CD,'r^')
plot(HLCD_us,NSCC_HLCD,'ko')

plot(0*ones(1,2),[0 Ymax],'r--')
for CYCind=1:NUM_cycles
   plot(CYCind/FF_Hz/1e-6*ones(1,2),[0 Ymax],'k:')
   plot(-CYCind/FF_Hz/1e-6*ones(1,2),[0 Ymax],'k:')
%    plot(CYCind/BF{2}/1e-6*ones(1,2),[0 Ymax],'g:')
%    plot(-CYCind/BF{2}/1e-6*ones(1,2),[0 Ymax],'g:')
end
plot([-XLIMIT XLIMIT],ones(1,2),'r--')

title(sprintf('# reps: %d,%d;   TOTAL #spikes = %d,%d;   AVGrate = %.2f,%.2f sp/sec;   SCC_AVGrate = %.2f sp/sec', ...
   length(SpikeTrains{1}),length(SpikeTrains{2}),TOTALspikesNSCC{1},TOTALspikesNSCC{2},AVGratesNSCC{1},AVGratesNSCC{2}, ...
   sqrt(AVGratesNSCC{1}*AVGratesNSCC{2})),'Interpreter','none')

xlabel('Delay (\mus)','Interpreter','tex')
ylabel('Normalized SAC')
hold off
xlim([-XLIMIT XLIMIT])
ylim([0 Ymax])
   
return;

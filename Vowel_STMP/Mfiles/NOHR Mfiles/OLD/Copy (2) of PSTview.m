function PIC=PSTview(picNums, excludeLines, figNum)
% M.Heinz 29Jul2004.  Modified from pst.m by GE
% Plots raster,pst,rate vs. rep, and PerHist for a given condition
% Allows multiple pst-style pics to be concatenated.  Also allows lines to be excluded for each picture.
% Assumes current directory is the data directory
%
% Usage: PSTview(picNums,excludeLines)
% picNums: vector of picture numbers
% excludeLines: cell array with vectors for each picture with any lines to be excluded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7/29/04 TODO
% *1) set up rate plot and perHist
% *2) finish laying out raster and pst
% Later:3) Setup general pst,raster,rateplot,PerHist to be called from anywhere based on spikes{1}
% Later:4) Setup general synch,phase calcs from PerHist, and add to plot
% Later:5) Setup basic Anal for 1 unit, phase vs freq (L)
%
% *6) Setup for vowels
% *7) Make basic plot for tone & vowel, to verify basic results
% 8) PLAN: plots to show main points
% 9) PLAN: data to collect, priorities!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/20/04 TODO
% *- Find Period for stimulus, (tone: f; vowel: F0)
% *- Compute PerHistogram from Period
% *- show raster, PSTH, PerH, FFT, RLFs
%%%%%%%%%%%%%% TODO 8/20/04 %%%%%%%%%%%%%%%%%%%%%%%%
% PERhist is GOOD
% *1) Need to plot Synch Rate plot from PerHist
% *2) Verify we get basically the same info as from PSThist
% *3) Calculate Synch and Phase for each feature?
%%%%%%%%%%%%%%%%%%%%%%%%%% 8/29/04 TODO %%%%%%%%%%%%%%%%%%%%%%%
% *1) store BF and FeatureFreqs_Hz_Hz in PIC
% *2) re-work PST and PERH with PIC-saved values, clean up
% *3) Calculate Synch and Phase for each feature
%%%%%%%%%% TODO 8/28/04
% *1) For now, plot both SynchRateFFT from PST and PERH
% *2) Calculate Synch and Phase from PERHFFT for each of 7 features
% 3) Setup more general function to do these calcs, and incorporate into plots

%%%%%%%%%%% TODO 8/30/04
% *1) Add Condition Text to LLcorner
% *2) List pictures used
% *3) Add Synch 0-1 as right axes
% *4) Generalize for ALL conditions (tones, FF) [e.g., TEMPLATE='EHrBF']
% *5) Compute Synch/Phase from PST SynchR
%%%%%%%%%%%% 9/7/04 TODO
% *5.5) Compute and Store Synch/Phase at F0
% 5.75 ADD ExcludeLines
% 6) Find good examples of effect to show (tones and F1)
% 7) Setup general Analysis of single unit: 
%        - database:
%        - Summary Plot for tones: r,s,p vs L, f
%        - extend to vowels
%        - setup SCC
% ?LATER?) How to deal with Features not at Harmonics, in terms of computing Synch/Phase???

clear FIG PICS PIC
global FIG PICS PIC

NUMpics=length(picNums);
if ~exist('excludeLines','var')
   for PICind=1:NUMpics
      excludeLines{PICind}=[];
   end
end
if ~exist('figNum','var')
   figNum=100;
end

PICS=cell(1,NUMpics);
for PICind=1:NUMpics
   PICS{PICind}.num = picNums(PICind);
   PICS{PICind}.excludeLines = excludeLines{PICind};
   [PICS{PICind}.x,errorMSG]=loadPic(PICS{PICind}.num);
   if ~isempty(errorMSG)
      error(errorMSG);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%% CHECK for ERRORS in each picture
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
concatPICs;

PIC.TEMPLATE=PIC.x.Stimuli.short_description(1:end-1);
if sum(strcmp(PIC.TEMPLATE,{'EHrBF','EHrFF'}))
   % Troughs and Formants
   PIC.FeatureFreqs_Hz=PIC.x.Stimuli.BASELINE.FeatFreqs_Hz*PIC.x.Stimuli.updateRate_Hz/PIC.x.Stimuli.origstimUpdateRate_Hz;
elseif sum(strcmp(PIC.TEMPLATE,{'TrFF','TrBF'}))
   % tone frequency   
   PIC.FeatureFreqs_Hz=PIC.x.Stimuli.main.freqs;
end

% Find Fundamental Frequency for Stimulus
if sum(strcmp(PIC.TEMPLATE,{'TrFF','TrBF'}))  % Tones
   PIC.FundamentalFreq_Hz=PIC.x.Stimuli.main.freqs;
elseif sum(strcmp(PIC.TEMPLATE,{'EHrFF','EHrBF'}))  % Vowels
   PIC.FundamentalFreq_Hz=PIC.x.Stimuli.BASELINE.F0_Hz*PIC.x.Stimuli.updateRate_Hz/PIC.x.Stimuli.origstimUpdateRate_Hz;
   PIC.FeatureFreqs_Hz=[PIC.FundamentalFreq_Hz PIC.FeatureFreqs_Hz];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(strcmp(PIC.TEMPLATE(end-1:end),'BF'))  % BF=BaseFreq
   PIC.BF_Hz=PIC.x.Stimuli.Condition.BaseFrequency_kHz*1000;
elseif sum(strcmp(PIC.TEMPLATE(end-1:end),'FF'))  % BF=NEED TO GET FROM TC!!!!!!!!!!!!!!!!!!!!!
   PIC.BF_Hz=input('Enter BF (Hz): ');
   disp('***TODO***: NEED TO GET BF from TC, rather than from BaseFrequency: FOR FF templates')
   PIC.FixedFreq_Hz=PIC.x.Stimuli.Condition.BaseFrequency_kHz*1000;
end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if(~exist('figNum', 'var'))
   FIG.num = 100;   % default, if not specified.
else
   FIG.num = figNum;
end
FIG.fontSize = 8; % used for axis labelling.
FIG.FeatureColors={'g','r'};  % Troughs,Formants
FIG.FeatureLines={':',':'};  % Troughs,Formants

layoutFigure;
do_raster;
do_rate;
do_PST_histo;
do_PerHist;

return;


%%################################################################################################
function concatPICs()
global PICS PIC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LATER: Verify all PICs are the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PIC=PICS{1};  % For now, just assume so
PIC=rmfield(PIC,'excludeLines');
for PICind=1:length(PICS)
   PIC.nums(PICind)=PICS{PICind}.num;
   PIC.excludeLines{PICind}=PICS{PICind}.excludeLines;
end
PIC=rmfield(PIC,'num');

%%%%%%%%%%%% Concatenate ALL SPIKES
% verify all spikes in PICS{1} are from fully presented lines
LASTgoodSpike=max(find(PIC.x.spikes{1}(:,1)==PIC.x.Stimuli.fully_presented_lines));
if isempty(LASTgoodSpike)  % if last spike is an earlier line than all presented
   LASTgoodSpike=length(PIC.x.spikes{1}(:,1));
end
% SpikeINDs=1:LASTgoodSpike;
SpikeINDs=find(~ismember(PICS{1}.x.spikes{1}(1:LASTgoodSpike,1)',PICS{1}.excludeLines));
PIC.x.spikes{1}=PICS{1}.x.spikes{1}(SpikeINDs,:);  %Take spikes from 1st of PICS
NumLines=length(unique(PIC.x.spikes{1}(:,1)));

for PICind=2:length(PICS)
   x=PICS{PICind}.x;

   LASTgoodSpike=max(find(x.spikes{1}(:,1)==x.Stimuli.fully_presented_lines));
   x.spikes{1}=x.spikes{1}(1:LASTgoodSpike,:);
   % Adjust Line Numbers
   x.spikes{1}(1:end,1)=x.spikes{1}(1:end,1)+NumLines;
   PIC.x.spikes{1}=[PIC.x.spikes{1};x.spikes{1}];
   NumLines=PIC.x.spikes{1}(end,1);
end
PIC.x.Stimuli.fully_presented_lines=NumLines;  % update for concatenated PIC

return;

%%################################################################################################
function layoutFigure()
global FIG PIC

x = PIC.x;

figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'inches', [8.7083    2.20    5.8333    8]};
FIG.handles.main = figure(FIG.num); clf;
set(gcf,figure_prop_name,figure_prop_val);

NameText=sprintf('''pst'' analysis for picture(s): p%04d, ', PIC.nums(1));
for PICind=2:length(PIC.nums)
   NameText=[NameText sprintf('p%04d ', PIC.nums(PICind))];
end
NameText=NameText(1:end-1);
set(gcf, 'Name', NameText);

Yshift=0.05;
rasterXcorner=0.35; rasterYcorner=0.65; rasterWidth=0.6; rasterHeight=0.3;
rateXcorner=0.07; rateYcorner=rasterYcorner; rateWidth=0.20; rateHeight=rasterHeight;
psthXcorner=rasterXcorner; psthWidth=rasterWidth; psthHeight=0.15; psthYcorner=rasterYcorner-psthHeight-Yshift; 
perhXcorner=rateXcorner; perhYcorner=psthYcorner; perhWidth=rateWidth; perhHeight=psthHeight;
synchrXcorner=rasterXcorner; synchrWidth=rasterWidth; synchrHeight=psthHeight; synchrYcorner=psthYcorner-synchrHeight-Yshift; ; 
DFTXcorner=rasterXcorner; DFTWidth=rasterWidth; DFTHeight=psthHeight; DFTYcorner=synchrYcorner-DFTHeight-Yshift; 
FIG.handles.rate   = subplot('Position',[rateXcorner rateYcorner rateWidth rateHeight]);
FIG.handles.raster = subplot('Position',[rasterXcorner rasterYcorner rasterWidth rasterHeight]);
FIG.handles.psth   = subplot('Position',[psthXcorner psthYcorner psthWidth psthHeight]);
FIG.handles.perh   = subplot('Position',[perhXcorner perhYcorner perhWidth perhHeight]);
FIG.handles.SynchR = subplot('Position',[synchrXcorner synchrYcorner synchrWidth synchrHeight]);
FIG.handles.SynchR_right = subplot('Position',[synchrXcorner+synchrWidth synchrYcorner .001 synchrHeight]);
FIG.handles.DFT    = subplot('Position',[DFTXcorner DFTYcorner DFTWidth DFTHeight]);
FIG.handles.DFT_right    = subplot('Position',[DFTXcorner+DFTWidth DFTYcorner 0.001 DFTHeight]);

titleString = sprintf('picture(s) %s recorded on %s\n%s', ...
               mat2str(PIC.nums),...
               x.General.date, ...
               x.Stimuli.description{1});
titleString = strrep(titleString, '_', '\_');
FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 1, titleString, 'Units', 'normalized', 'FontSize', FIG.fontSize-1,...
   'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

OffsetText={'+','-'};
if strcmp(PIC.TEMPLATE(end-1:end),'BF')
   OffsetIND=strcmp(x.Stimuli.Condition.Offset_Direction,'below ')+1;
end

if strcmp(PIC.TEMPLATE,'EHrBF');
   if isnumeric(x.Stimuli.Condition.FreqOffset_octs)
      x.Stimuli.Condition.FreqOffset_octs=sprintf('%s%s',num2str(x.Stimuli.Condition.FreqOffset_octs),' ');
   end
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %s%soctaves\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      OffsetText{OffsetIND}, ...
      x.Stimuli.Condition.FreqOffset_octs, ...
      x.Stimuli.Condition.Feature, ...
      x.Stimuli.Condition.FormsAtHarmonics, ...
      x.Stimuli.Condition.InvertPolarity, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'TrBF')
   if isnumeric(x.Stimuli.Condition.FreqOffset_octs)
      x.Stimuli.Condition.FreqOffset_octs=sprintf('%s%s',num2str(x.Stimuli.Condition.FreqOffset_octs),' ');
   end
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nOffset = %s%soctaves\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      OffsetText{OffsetIND}, ...
      num2str(x.Stimuli.Condition.FreqOffset_octs), ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'TrFF')
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nFF = %.3f kHz\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      x.Stimuli.Condition.BaseFrequency_kHz, ...
      x.Stimuli.Condition.Level_dBSPL);
elseif strcmp(PIC.TEMPLATE,'EHrFF');
   ConditionString = sprintf('Template: %s\n\nBF = %.3f kHz\nFF = %.3f kHz\nFeature: %s\nForms at Harmonics: %s\nPolarity Inverted: %s\nLevel= %.f dB SPL\n', ...
      PIC.TEMPLATE, ...
      PIC.BF_Hz/1000, ...
      x.Stimuli.Condition.BaseFrequency_kHz, ...
      x.Stimuli.Condition.Feature, ...
      x.Stimuli.Condition.FormsAtHarmonics, ...
      x.Stimuli.Condition.InvertPolarity, ...
      x.Stimuli.Condition.Level_dBSPL);
end
text(0.03, synchrYcorner+synchrHeight-Yshift/3, ConditionString, 'Units', 'normalized', 'FontSize', FIG.fontSize,...
   'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','b');

   
   
return;

%%################################################################################################
function do_raster()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.raster);

plot(x.spikes{1}(:,2),x.spikes{1}(:,1), 'k.', 'MarkerSize', 4);

FIG.raster.xmax = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
FIG.raster.ymax = ceil(x.Stimuli.fully_presented_lines/10)*10;
xlim([0 FIG.raster.xmax]);
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
FIG.raster.hx=xlabel('time (sec)');
set(FIG.raster.hx,'units','norm','pos',[0.4926   -0.01         0])
set(gca,'YTick',0:10:FIG.raster.ymax)
set(gca, 'TickDir', 'out');
return;

%%################################################################################################
function do_rate()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.rate);

% compute driven and spontaneous rates for each line
spikeTimes = x.spikes{1};
driven_dur_sec = x.Hardware.Trigger.StmOn / 1000;
spont_dur_sec = x.Hardware.Trigger.StmOff / 1000;
for line = 1:x.Stimuli.fully_presented_lines
   spikeIndices = find( (spikeTimes(:,1) == line ) );
   driv(line) = length (find(spikeTimes(spikeIndices,2) <= driven_dur_sec));
   spont(line) = length(spikeIndices) - driv(line);
end
driv = driv / driven_dur_sec;
spont = spont / spont_dur_sec;

plot(driv,1:length(driv),'r','LineWidth',2)
hold on
plot(spont,1:length(spont),'g','LineWidth',2)
hold off
ylim([0 FIG.raster.ymax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca,'XDir','reverse')
set(gca,'YTick',0:10:FIG.raster.ymax)
YLabel('Rep Number');
FIG.rate.hx=xlabel('(sp/s)');
set(FIG.rate.hx,'units','norm','pos',[-0.2   -0.018         0])

PIC.rate.driv=driv;
PIC.rate.spont=spont;

return;


%%################################################################################################
function do_PST_histo()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.psth);

% FIG.psth.binWidth_sec = 100e-6; % matches Miller et al 1999a
FIG.psth.binWidth_sec = 20e-6; % matches Wong et al 1998
FIG.psth.StartTime_sec=20/1000;  % Skip onset transient
FIG.psth.EndTime_sec=PIC.x.Hardware.Trigger.StmOn/1000;

lastBin_sec = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
pst_X = [0:FIG.psth.binWidth_sec:lastBin_sec];
pst_Y=hist(x.spikes{1}(:,2), pst_X);
%%% Convert PST to spikes per second
pst_Y=pst_Y/PIC.x.Stimuli.fully_presented_lines/FIG.psth.binWidth_sec; % Convert to sp/sec

%%%%%%%%%%%%%%%%%%% Window for Calculating Synch Rate
drivenSpikeIndices = find((pst_X >= FIG.psth.StartTime_sec)&(pst_X <= FIG.psth.EndTime_sec));
NumDrivenSpikes=sum(pst_Y(drivenSpikeIndices)*PIC.x.Stimuli.fully_presented_lines*FIG.psth.binWidth_sec);

plot(pst_X, pst_Y, 'k');
hold on
plot(pst_X(drivenSpikeIndices), pst_Y(drivenSpikeIndices), 'b');
xlim([0 FIG.raster.xmax]);
set(gca, 'FontSize', FIG.fontSize);
set(gca, 'TickDir', 'out');
FIG.psth.hy=ylabel('(sp/s)');
set(FIG.psth.hy,'units','norm','pos',[-0.108    0.4874         0])
FIG.psth.hx=xlabel('time (sec)');
set(FIG.psth.hx,'units','norm','pos',[0.4926   -0.01         0])

PIC.pst.pst_X=pst_X;
PIC.pst.pst_Y=pst_Y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find DFT - to get Synch and phase at SynchFreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SynchRateCFFT=fft(pst_Y(drivenSpikeIndices))/length(drivenSpikeIndices);
FFTfreqs=(0:length(SynchRateCFFT)-1)*(1/FIG.psth.binWidth_sec)/length(SynchRateCFFT);

PIC.SynchR.NumDrivenSpikes=NumDrivenSpikes;
PIC.SynchR.SynchRateCFFT=SynchRateCFFT;
PIC.SynchR.FFTfreqs=FFTfreqs;

subplot(FIG.handles.SynchR);
plot(FFTfreqs,abs(SynchRateCFFT))
hold on
% ADD Synch Coef Labels to right Axis
set(FIG.handles.SynchR_right,'YAxisLocation','right','YTick',[0:.2:1],'FontSize',FIG.fontSize,'TickDir','out','TickLength',[.03 0.025])

% Plot BF 
plot(PIC.BF_Hz*ones(1,2),[0 1e4],'k--')
FIG.SynchR.ymax=max(abs(SynchRateCFFT));
text(PIC.BF_Hz,FIG.SynchR.ymax*1.06,'BF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')

% Label F0
text(PIC.FundamentalFreq_Hz,FIG.SynchR.ymax*.98,'F0','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center','Color',FIG.FeatureColors{2})

% Plot FF (Fixed Frequency, if used)
if strcmp(PIC.TEMPLATE(end-1:end),'FF')
   plot(PIC.FixedFreq_Hz*ones(1,2),[0 1e4],'k-.')
   text(PIC.FixedFreq_Hz,FIG.SynchR.ymax*1.06,'FF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
end

% Plot Features (e.g., Tone Freq, Formants and Troughs)
for i=1:length(PIC.FeatureFreqs_Hz)
   plot(PIC.FeatureFreqs_Hz(i)*ones(1,2),[0 1e4],':','Color',FIG.FeatureColors{mod(i,2)+1})
end

FIG.SynchR.xmax=min([ceil(max([PIC.FeatureFreqs_Hz PIC.BF_Hz])/5000)*5000 20000]);

xlim([0 FIG.SynchR.xmax])
ylim([0 FIG.SynchR.ymax])
set(gca, 'FontSize', FIG.fontSize);
% xlabel('Frequency (Hz)');
ylabel('Synchronized Rate (sp/sec)')
hold off


PIC.SynchR.Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
%%%%%%%%%%% Compute Synchs and Phases for each feature
FeatureINDs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureSynchs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeaturePhases=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureRaySig=zeros(size(PIC.FeatureFreqs_Hz));
for i=1:length(PIC.FeatureFreqs_Hz)
   %% Here, the FFTfreqs are not equal to the harmonics, and so we need to pick the closest one
   %    FEATind=find(round(FFTfreqs*1000)==round(PIC.FeatureFreqs_Hz(i)*1000));
   [yyy,FEATind]=min(abs(FFTfreqs-PIC.FeatureFreqs_Hz(i)));
   if ~isempty(FEATind)
      FeatureINDs(i)=FEATind;
      FeatureSynchs(i)=abs(SynchRateCFFT(FeatureINDs(i)))/SynchRateCFFT(1);
      FeaturePhases(i)=angle(SynchRateCFFT(FeatureINDs(i)));
      RayleighStat=2*NumDrivenSpikes*FeatureSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
      RayleighCRIT=chi2inv(1-PIC.SynchR.Rayleigh_P,2);
      if RayleighStat>RayleighCRIT
         FeatureRaySig(i)=1;
      end
   end
end

PIC.SynchR.FeatureINDs=FeatureINDs;
PIC.SynchR.FeatureSynchs=FeatureSynchs;
PIC.SynchR.FeaturePhases=FeaturePhases;
PIC.SynchR.FeatureRaySig=FeatureRaySig;


return;


%%################################################################################################
function do_PerHist()
global FIG PIC

x = PIC.x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DFT usage for getting Synch/Phase %%%%%%%%%%%
% K=16 used for all frequencies: estimates 1st 7 harmonics, which is plenty for us.
% Johnson (1980) used 30-200 bins/cycle for the PerHist;
% Anderson et al (1971) used 10 bins/cycle for PerHist.
% With K bins/cycle for all frequencies, Synch/Phase from PerHist DFT is identical to from PSTHist.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=256;  %% # bins/cycle
BINwidth=1/PIC.FundamentalFreq_Hz/K;
FIG.perh.binWidth_sec=BINwidth;
M=floor((FIG.psth.EndTime_sec-FIG.psth.StartTime_sec)*PIC.FundamentalFreq_Hz);  % Integer number of cycles to include in spike window
EndTime=FIG.psth.StartTime_sec+M/PIC.FundamentalFreq_Hz; % Limit to integer number of cycles

% Find driven spikes to use
spikeTimes = x.spikes{1};
drivenSpikeIndices = find( (spikeTimes(:,2) >= FIG.psth.StartTime_sec)&(spikeTimes(:,2) <= EndTime) );
drivenSpikes=spikeTimes(drivenSpikeIndices,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Period Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivenSpikes_BINS=rem(floor(drivenSpikes/BINwidth),K)+1; %Convert times to PerHist bins (1:K)
[PerHist,xxx]=hist(drivenSpikes_BINS,(1:K));  % Make Histogram from BINS (1:K)
NumDrivenSpikes=sum(PerHist);

%%% Convert PerHist to spikes per second
PerHist=PerHist/PIC.x.Stimuli.fully_presented_lines/M/FIG.perh.binWidth_sec; % Convert to sp/sec

PIC.PerH.PerHist=PerHist;

subplot(FIG.handles.perh);
bar((1:K)-.5,PerHist,1)
xlim([0 K])
set(gca, 'FontSize', FIG.fontSize);
text(.05,.92,sprintf('F0 = %.f (Hz)',PIC.FundamentalFreq_Hz),'units','norm', 'FontSize', FIG.fontSize)
hy=ylabel('(sp/s)');
set(hy,'units','norm')
set(hy,'pos',[-0.27 0.4914 0])
set(gca,'XTick',[0 .5 1]*K)
set(gca,'XTickLabel',[0 .5 1]*K*FIG.perh.binWidth_sec*1000)
set(gca, 'TickDir', 'out');
xlabel('time (msec)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find DFT - to get Synch and phase at PIC.FundamentalFreq_Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SynchRateDFT=fft(PerHist)/length(PerHist);
FFTfreqs=(0:length(SynchRateDFT)-1)*(1/FIG.perh.binWidth_sec)/length(SynchRateDFT);

PIC.DFT.NumDrivenSpikes=NumDrivenSpikes;
PIC.DFT.SynchRateDFT=SynchRateDFT;
PIC.DFT.FFTfreqs=FFTfreqs;

subplot(FIG.handles.DFT);
%%%%%%% Get both Synhchronized Rate and SynchCoef labels
plot(FFTfreqs,abs(SynchRateDFT),'b-x')
hold on
% Plot BF (Target Frequency)
FIG.DFT.ymax=max(abs(SynchRateDFT));
text(PIC.BF_Hz,FIG.DFT.ymax*1.06,'BF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
plot(PIC.BF_Hz*ones(1,2),[0 1e4],'k--')

% Label F0
text(PIC.FundamentalFreq_Hz,FIG.SynchR.ymax*.98,'F0','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center','Color',FIG.FeatureColors{2})

% Plot FF (Fixed Frequency, if used)
if strcmp(PIC.TEMPLATE(end-1:end),'FF')
   plot(PIC.FixedFreq_Hz*ones(1,2),[0 1e4],'k-.')
   text(PIC.FixedFreq_Hz,FIG.DFT.ymax*1.06,'FF','FontSize',FIG.fontSize,'Units','data','VerticalAlignment','bottom','HorizontalAlignment','center')
end

% Plot Features (e.g., Tone Freq, Formants and Troughs)
for i=1:length(PIC.FeatureFreqs_Hz)
   plot(PIC.FeatureFreqs_Hz(i)*ones(1,2),[0 1e4],':','Color',FIG.FeatureColors{mod(i,2)+1})
end

xlim([0 FIG.SynchR.xmax])
ylim([0 FIG.DFT.ymax])
set(gca, 'FontSize', FIG.fontSize);
xlabel('Frequency (Hz)');
ylabel('Synchronized Rate (sp/sec)')

% ADD Synch Coef Labels to right Axis
set(FIG.handles.DFT_right,'YAxisLocation','right','YTick',[0:.2:1],'FontSize',FIG.fontSize,'TickDir','out','TickLength',[.03 0.025])
subplot(FIG.handles.DFT_right)
ht=title('Synch','FontSize',FIG.fontSize);
set(ht,'units','norm','position',[7    1.0427 0])

PIC.DFT.Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
%%%%%%%%%%% Compute Synchs and Phases for each feature
FeatureINDs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureSynchs=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeaturePhases=NaN+zeros(size(PIC.FeatureFreqs_Hz));
FeatureRaySig=zeros(size(PIC.FeatureFreqs_Hz));
for i=1:length(PIC.FeatureFreqs_Hz)
   FEATind=find(round(FFTfreqs*1000)==round(PIC.FeatureFreqs_Hz(i)*1000));
   if ~isempty(FEATind)
      FeatureINDs(i)=FEATind;
      FeatureSynchs(i)=abs(SynchRateDFT(FeatureINDs(i)))/SynchRateDFT(1);
      FeaturePhases(i)=angle(SynchRateDFT(FeatureINDs(i)));
      RayleighStat=2*NumDrivenSpikes*FeatureSynchs(i)^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
      RayleighCRIT=chi2inv(1-PIC.DFT.Rayleigh_P,2);
      if RayleighStat>RayleighCRIT
         FeatureRaySig(i)=1;
      end
   end
end

PIC.DFT.FeatureINDs=FeatureINDs;
PIC.DFT.FeatureSynchs=FeatureSynchs;
PIC.DFT.FeaturePhases=FeaturePhases;
PIC.DFT.FeatureRaySig=FeatureRaySig;
SigText={' ','*'};

%%% List: Numspikes,PIC.FundamentalFreq_Hz,Synch,Phase,Signif
Features={'F0','T0','F1','T1','F2','T2','F3','T3'};
if sum(strcmp(PIC.TEMPLATE,{'EHrBF','EHrFF'}))
   FeatIND=find(strcmp(Features,x.Stimuli.Condition.Feature(1:2)));
   SynchText=sprintf('Nsps=%d, %s =%.2f Hz, Synch=%.2f[%s], Ph= %.2f rad', ...
      NumDrivenSpikes,Features{FeatIND},PIC.FeatureFreqs_Hz(FeatIND),FeatureSynchs(FeatIND),SigText{FeatureRaySig(FeatIND)+1},FeaturePhases(FeatIND));
elseif sum(strcmp(PIC.TEMPLATE,{'TrBF','TrFF'}))
   SynchText=sprintf('Nsps=%d, Tone =%.2f Hz, Synch=%.2f[%s], Ph= %.2f rad', ...
      NumDrivenSpikes,PIC.FundamentalFreq_Hz,FeatureSynchs(1),SigText{FeatureRaySig(1)+1},FeaturePhases(1));
   FeatIND=1;
end
subplot(FIG.handles.psth);
FIG.psth.htext=text(0.99,0.96,SynchText,'FontSize',FIG.fontSize,'HorizontalAlign','right', ...
   'Units','norm','VerticalAlign','top');

subplot(FIG.handles.figBox)
clear FullSynchText
% if sum(strcmp(PIC.TEMPLATE,{'EHrBF','EHrFF'}))
   FullSynchText{1}=sprintf('Feature   Synch   Phase (rad)');
   for i=1:length(PIC.FeatureFreqs_Hz)
      FullSynchText{i+1} = sprintf('    %s%s        %.2f%s       %.2f',Features{i},SigText{(i==FeatIND)+1},FeatureSynchs(i),SigText{FeatureRaySig(i)+1},FeaturePhases(i));
   end
% elseif strcmp(PIC.TEMPLATE,'TrBF')
%    FullSynchText = 'test';
% end
text(0.01, .03,strvcat(FullSynchText), 'Units', 'normalized', 'FontSize', FIG.fontSize,...
   'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color','k');
% text(0.03, synchrYcorner+synchrHeight-Yshift, FullSynchText, 'Units', 'normalized', 'FontSize', FIG.fontSize,...
%    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top','Color','k');

%%%%%%%%%%%%%

return;



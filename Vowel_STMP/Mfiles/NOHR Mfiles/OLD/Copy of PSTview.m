function PSTview(picNums, excludeLines, figNum)
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
% 7) Make basic plot for tone & vowel, to verify basic results
% 8) PLAN: plots to show main points
% 9) PLAN: data to collect, priorities!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8/20/04 TODO
% *- Find Period for stimulus, (tone: f; vowel: F0)
% *- Compute PerHistogram from Period
% *- show raster, PSTH, PerH, FFT, RLFs


%%%%%%%%%%%%%% TODO 8/20/04 %%%%%%%%%%%%%%%%%%%%%%%%
% PERhist is GOOD
% 1) Need to plot Synch Rate plot from PerHist
% 2) Verify we get basically the same info as from PSThist
% 3) Calculate Synch and Phase for each feature?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









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
   [PICS{PICind}.x,errorMSG]=loadPic(PICS{PICind}.num);
   if ~isempty(errorMSG)
      error(errorMSG);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%% CHECK for ERRORS in each picture
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
concatPICs;

if(~exist('figNum', 'var'))
   FIG.num = 100;   % default, if not specified.
else
   FIG.num = figNum;
end
% FIG.psth.binWidth_sec = 100e-6; % matches Miller et al 1999a
FIG.psth.binWidth_sec = 20e-6; % matches Wong et al 1998
FIG.fontSize = 8; % used for axis labelling.


FIG.psth.StartTime_sec=20/1000;  % Skip onset transient
FIG.psth.EndTime_sec=PIC.x.Hardware.Trigger.StmOn/1000;

layoutFigure;
do_raster;
do_rate
do_PST_histo;
do_PerHist

return;


%%################################################################################################
function concatPICs()
global PICS PIC

%%%%%%%%%%%%% LATER: Verify all PICs are the same
PIC=PICS{1};  % For now, just assume so
for PICind=1:length(PICS)
   PIC.nums(PICind)=PICS{PICind}.num;
end
PIC=rmfield(PIC,'num');

%%%%%%%%%%%% Concatenate ALL SPIKES
% verify all spikes in PICS{1} are from fully presented lines
LASTgoodSpike=max(find(PIC.x.spikes{1}(:,1)==PIC.x.Stimuli.fully_presented_lines));
if isempty(LASTgoodSpike)  % if last spike is an earlier line than all presented
   LASTgoodSpike=length(PIC.x.spikes{1}(:,1));
end
PIC.x.spikes{1}=PIC.x.spikes{1}(1:LASTgoodSpike,:);
NumLines=PIC.x.spikes{1}(LASTgoodSpike,1);
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
figure_prop_val =  { 'auto'            ,'inches', [8.7083    4.0208    5.8333    6.1667]};
FIG.handles.main = figure(FIG.num); clf;
set(gcf,figure_prop_name,figure_prop_val);

NameText=sprintf('''pst'' analysis for pictures: p%04d, ', PIC.nums(1));
for PICind=2:length(PIC.nums)
   NameText=[NameText sprintf('p%04d ', PIC.nums(PICind))];
end
NameText=NameText(1:end-1);
set(gcf, 'Name', NameText);

rasterXcorner=0.35; rasterYcorner=0.55; rasterWidth=0.6; rasterHeight=0.4;
rateXcorner=0.07; rateYcorner=rasterYcorner; rateWidth=0.20; rateHeight=rasterHeight;
psthXcorner=rasterXcorner; psthWidth=rasterWidth; psthHeight=0.2; psthYcorner=rasterYcorner-psthHeight-0.04; 
perhXcorner=rateXcorner; perhYcorner=psthYcorner; perhWidth=rateWidth; perhHeight=psthHeight;
synchrXcorner=rasterXcorner; synchrWidth=rasterWidth; synchrHeight=0.2; synchrYcorner=psthYcorner-synchrHeight-0.04; ; 
FIG.handles.rate   = subplot('Position',[rateXcorner rateYcorner rateWidth rateHeight]);
FIG.handles.raster = subplot('Position',[rasterXcorner rasterYcorner rasterWidth rasterHeight]);
FIG.handles.psth   = subplot('Position',[psthXcorner psthYcorner psthWidth psthHeight]);
FIG.handles.perh   = subplot('Position',[perhXcorner perhYcorner perhWidth perhHeight]);
FIG.handles.SynchR = subplot('Position',[synchrXcorner synchrYcorner synchrWidth synchrHeight]);

textString = sprintf('picture %03d recorded on %s\n%s', ...
               x.General.picture_number,...
               x.General.date, ...
               x.Stimuli.description{1});
textString = strrep(textString, '_', '\_');
FIG.handles.figBox = axes('Position',[0 0 1 1],'Visible','off');
text(0.5, 1, textString, 'Units', 'normalized', 'FontSize', FIG.fontSize-1,...
      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
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
% YLabel('raster');
FIG.raster.hx=xlabel('time (sec)');
set(FIG.raster.hx,'units','norm','pos',[0.4926   -0.01         0])
set(gca,'YTick',0:10:FIG.raster.ymax)
% set(gca, 'YTickLabel', '');
% set(gca, 'XTickLabel', '');
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
% text(.95,.95,'Rate','Horiz','right','units','norm')
% set(gca, 'XTickLabel', '');
% set(gca, 'YTickLabel', '');
% set(gca, 'TickDir', 'out');
return;


%%################################################################################################
function do_PST_histo()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.psth);

lastBin_sec = (x.Hardware.Trigger.StmOn + x.Hardware.Trigger.StmOff) / 1000;
pst_X = [0:FIG.psth.binWidth_sec:lastBin_sec];
pst_Y=hist(x.spikes{1}(:,2), pst_X);
%%% Convert PST to spikes per second
pst_Y=pst_Y/PIC.x.Stimuli.fully_presented_lines/FIG.psth.binWidth_sec; % Convert to sp/sec

%%%%%%%%%%%%%%%%%%% Window for Calculating Synch Rate
drivenSpikeIndices = find((pst_X >= FIG.psth.StartTime_sec)&(pst_X <= FIG.psth.EndTime_sec));
% NumDrivenSpikes=sum(pst_Y(drivenSpikeIndices));

plot(pst_X, pst_Y, 'k');
hold on
plot(pst_X(drivenSpikeIndices), pst_Y(drivenSpikeIndices), 'b');
xlim([0 FIG.raster.xmax]);
set(gca, 'FontSize', FIG.fontSize);
% YLabel('PSTH');
% set(gca, 'XTickLabel', '');
% set(gca, 'YTickLabel', '');
set(gca, 'TickDir', 'out');
FIG.psth.hy=ylabel('(sp/s)');
set(FIG.psth.hy,'units','norm','pos',[-0.108    0.4874         0])
FIG.psth.hx=xlabel('time (sec)');
set(FIG.psth.hx,'units','norm','pos',[0.4926   -0.01         0])



%%%%%%%%%% MOVE THIS TO PERHISTOGRAM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find DFT - to get Synch and phase at FundamentalFreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SynchRateCFFT=fft(pst_Y(drivenSpikeIndices))/length(drivenSpikeIndices);
FFTfreqs=(0:length(SynchRateCFFT)-1)*(1/FIG.psth.binWidth_sec)/length(SynchRateCFFT);

subplot(FIG.handles.SynchR);
plot(FFTfreqs,abs(SynchRateCFFT))
hold on
% Plot BF (Target Frequency)
plot(x.Stimuli.Condition.BaseFrequency_kHz*1000*ones(1,2),[0 1e4],'k--')
if isfield(x.Stimuli.Condition,'Feature')
   % Plot Formants
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(2)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'r:')
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(4)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'r:')
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(6)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'r:')
   % Plot Troughs
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(1)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'g:')
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(3)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'g:')
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(5)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'g:')
   plot(x.Stimuli.BASELINE.FeatFreqs_Hz(7)*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz*ones(1,2),[0 1e4],'g:')
else
   % PLot tone frequency   
   plot(x.Stimuli.main.freqs*ones(1,2),[0 1e4],'r:')
end
xlim([0 5000])
ylim([0 max(abs(SynchRateCFFT))])
set(gca, 'FontSize', FIG.fontSize);
xlabel('Frequency (Hz)');
ylabel('Synchronized Rate (sp/sec)')




% Synch=abs(DFT(2))/DFT(1);
% Phase=angle(DFT(2));
% RayleighStat=2*NumDrivenSpikes*Synch^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
% Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
% RayleighCRIT=chi2inv(1-Rayleigh_P,2);
% if RayleighStat>RayleighCRIT
%    SignifText='*';
% else
%    SignifText='x';
% end
%%%%%%%%%%%%%



return;


%%################################################################################################
function do_PerHist()
global FIG PIC

x = PIC.x;
subplot(FIG.handles.perh);

% Find Fundamental Frequency for Stimulus
if isfield(x.Stimuli,'main')
   % Tones
   FundamentalFreq=x.Stimuli.main.freqs;
else
   % Vowels
   FundamentalFreq=x.Stimuli.BASELINE.F0_Hz*x.Stimuli.updateRate_Hz/x.Stimuli.origstimUpdateRate_Hz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% DFT usage for getting Synch/Phase %%%%%%%%%%%
% K=16 used for all frequencies: estimates 1st 7 harmonics, which is plenty for us.
% Johnson (1980) used 30-200 bins/cycle for the PerHist;
% Anderson et al (1971) used 10 bins/cycle for PerHist.
% With K bins/cycle for all frequencies, Synch/Phase from PerHist DFT is identical to from PSTHist.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=256;  %% # bins/cycle
BINwidth=1/FundamentalFreq/K;
M=floor((FIG.psth.EndTime_sec-FIG.psth.StartTime_sec)*FundamentalFreq);  % Integer number of cycles to include in spike window
EndTime=FIG.psth.StartTime_sec+M/FundamentalFreq; % Limit to integer number of cycles

% Find driven spikes to use
spikeTimes = x.spikes{1};
drivenSpikeIndices = find( (spikeTimes(:,2) >= FIG.psth.StartTime_sec)&(spikeTimes(:,2) <= EndTime) );
drivenSpikes=spikeTimes(drivenSpikeIndices,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Period Histogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drivenSpikes_BINS=rem(floor(drivenSpikes/BINwidth),K)+1; %Convert times to PerHist bins (1:K)
[PerHist,x]=hist(drivenSpikes_BINS,(1:K));  % Make Histogram from BINS (1:K)
NumDrivenSpikes=sum(PerHist);



%%%%%%%%%%%%%% TODO 8/20/04 %%%%%%%%%%%%%%%%%%%%%%%%
% PERhist is GOOD
% 1) Need to plot Synch Rate plot from PerHist
% 2) Verify we get basically the same info as from PSThist
% 3) Calculate Synch and Phase for each feature?







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find DFT - to get Synch and phase at FundamentalFreq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFT=fft(PerHist);
Synch=abs(DFT(2))/DFT(1);
Phase=angle(DFT(2));
RayleighStat=2*NumDrivenSpikes*Synch^2;  % Rayleigh criterion for significance of Synch/Phase coefficients
Rayleigh_P=0.001;  %Confidence in Synch/Phase estimates: P<Rayleigh_P;
RayleighCRIT=chi2inv(1-Rayleigh_P,2);
if RayleighStat>RayleighCRIT
   SignifText='*';
else
   SignifText='x';
end
%%%%%%%%%%%%%

bar((1:K)-.5,PerHist,1)
xlim([0 K])
set(gca, 'FontSize', FIG.fontSize);
text(.05,.92,sprintf('F0 = %.f (Hz)',FundamentalFreq),'units','norm', 'FontSize', FIG.fontSize)
YLabel('PerHist');
set(gca, 'XTickLabel', '');
set(gca, 'YTickLabel', '');
set(gca, 'TickDir', 'out');


%%% List: Numspikes,FundamentalFreq,Synch,Phase,Signif
% [y,TextBin]=min(PerHist);
subplot(FIG.handles.psth);
FIG.psth.htext=text(0.99,0.96,sprintf('Nspikes=%d, Freq=%.2f Hz, Synch=%.2f%s, Phase=%.2f rad', ...
   NumDrivenSpikes,FundamentalFreq,Synch,SignifText,Phase),'FontSize',FIG.fontSize,'HorizontalAlign','right', ...
   'Units','norm','VerticalAlign','top');

return;


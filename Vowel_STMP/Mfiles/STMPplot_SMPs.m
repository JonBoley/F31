function STMPplot_RPSs(ExpDate,UnitName,STIMtype,DATAtype)
% FROM: STMPplot_DFTs(ExpDate,UnitName,STIMtype,DATAtype)
% M. Heinz Jun 11, 2007
% Plots Rate, Synch, and Phase data to match STMP POPULATION PERhist and
% DFT plots.  Assumes STMP format for data storage, but this can be either:
% DATAtype=
%     - 1unitNstim (DFTs.Info.STMPshift=0)
%     - Nchans1stim (DFTs.Info.STMPshift=1)
%
% Setup to be general for STMP data, but some specifics for different STIM
% types.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
global FeaturesText

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=1;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03),
%%%% Specify ExpDate if not provided
if ~exist('ExpDate','var')
   if TESTunitNUM==1
      ExpDate='111804'; UnitName='1.28'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';
   elseif TESTunitNUM==2
      ExpDate='041805'; UnitName='2.04';
   elseif TESTunitNUM==3
      ExpDate='111606'; UnitName='2.09';
   elseif TESTunitNUM==4
      ExpDate='041307'; UnitName='1.03';
   end
end

%%%% Find the full Experiment Name
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(STMP_ExpList)
   if ~isempty(strfind(STMP_ExpList{i},ExpDateText))
      ExpName=STMP_ExpList{i};
      break;
   end
end
if ~exist('ExpName','var')
   disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
   disp(strvcat(STMP_ExpList))
   beep
   error('STOPPED');
end

%%%% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Setup related directories
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)
disp(sprintf('... STMP: Plotting ''%s'' RPSs  for:  Experiment: ''%s''; Unit: %d.%02d',DATAtype,ExpName,TrackNum,UnitNum))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify DFTs exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DFTs_filename=sprintf('DFTs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
DFTs_dataname=sprintf('DFTs_%s',DATAtype);
eval(['load ''' fullfile(STMPanal_dir,DFTs_filename) ''''])
eval(['DFTs = ' DFTs_dataname ';']);
eval(['clear ' DFTs_dataname]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this data, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(DFTs.DFTs.DFTs{1,1,1},1);
NumP=size(DFTs.DFTs.DFTs{1,1,1},2);
NumFEATURES=size(DFTs.DFTs.DFTs,1);
NumPOL=size(DFTs.DFTs.DFTs,2);
NumHARMS=size(DFTs.DFTs.DFTs,3);

BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave

XLIMITS_rate=[0 300];
XLIMITS_synch=[0 1];
XLIMITS_phase=[-1 1];
XLIMITS_Nspikes=[0 3000];

FeatureColors={'r','g'};  % 'Formants','Troughs'
PARAMcolors={'b','r','g','c','k','y'};
FIG.LEGfontSize=12;
FIG.ANNOTfontSize=12;
FIG.LABELfontSize=10;
FIG.TITLEfontSize=10;
FIG.FEATfontSize=14;
FIG.VERTLABfontSize=10;
FIG.HORLABfontSize=14;

%%%%%%%%%%%%%%%%%%%%%%%%% TODO RPS
% *1) scale axes Ywidth as PERhist
% *2) Track down phase/synch - why is phase NaN, but synch shown
% *3) Why so much diff in PHASE re: ARO??? Neural Delay: YES!!!!!
% *4) ADD F0 plot at end
% *5) ADD NUMspikes plot at end
% *6) clean up code at END

NUMcols=3;  % R,S,P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FEATURE plot  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
set(gcf,'DefaultTextInterpreter','none')
set(gcf,'units','norm','pos',[0.55   0.07    0.45    0.85])
set(gcf,'Name',sprintf('RPSs(Feat) for Unit: %s',UnitName))
set(gcf,'NumberTitle','off')

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         FIGinfo.title=sprintf('%s: [FEAT] Rate, Synch, Phase Plots [FormsAtHarms: %s, PolarityInverted: %s]\nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
            STIMtype,DFTs.ConditionInfo.FormsAtHarmonics{HarmIND},DFTs.ConditionInfo.polarities{PolIND}, ...
            ExpName,TrackNum,UnitNum,DFTs.Info.BF_kHz,DFTs.Info.Threshold_dBSPL,DFTs.Info.SR_sps,DFTs.Info.Q10);
         FIGinfo.CONDlabel_text=sprintf('Feature: %s',DFTs.ConditionInfo.features{FeatIND});
        
         CHdata.CHvals=DFTs.ChannelInfo.channel_values{FeatIND,PolIND,HarmIND};
         %% FIND BF channel
         % BFind=find(abs(log2(CHdata.CHvals/DFTs.Info.BF_kHz))<BFoctCRIT);
         BFind=find(DFTs.ChannelInfo.Octave_Shifts{FeatIND,PolIND,HarmIND}==0);

         CHdata.Rdata=DFTs.Feature_RSP.rate_sps{FeatIND,PolIND,HarmIND}; % Rate
         CHdata.Sdata=DFTs.Feature_RSP.synch{FeatIND,PolIND,HarmIND}; % Synch
         CHdata.Pdata=DFTs.Feature_RSP.phase_cyc{FeatIND,PolIND,HarmIND}; % Phase
         CHdata.Adata=DFTs.ALSR.ALSRs_sps{FeatIND,PolIND,HarmIND}(BFind,:); % ALSR

         %% Find appropriate YLIMITS to MATCH PERhists: based on channel_values and chGAIN
         chGAIN=2; % # of channels covered by max PERhist - USED TO SCALE CHdataY - NEEDS TO MATCH PERhist
         NumCH=length(CHdata.CHvals);         highCH=max(CHdata.CHvals);         lowCH=min(CHdata.CHvals);
         ch_logCHwidth=log10(highCH/lowCH)/(NumCH-1);  % log10 channel width
         ch_YMAX=10^((NumCH-1+chGAIN)*ch_logCHwidth)*lowCH;    % sets an extra (GAIN-1) log channel widths (to allow for chGAIN>1)
         YLIMITS=[lowCH ch_YMAX];  % Used for all plots
         YTICKS=sort(round(CHdata.CHvals*100)/100);

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Setup Horizontal labels/lines for feature frequencies
         if sum(strcmp(STIMtype,{'EHrBFi','EHvNrBFi'}))
            ii=0;
            for FeatIND2=1:find(strcmp(FeaturesText,'T3'))
               ii=ii+1;
               FIGinfo.HorizontalLabels.Text{ii}=FeaturesText{FeatIND2};
               FIGinfo.HorizontalLabels.Yvals{ii}=DFTs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}(FeatIND2)/1000;
               if strcmp(FeaturesText{FeatIND2},'F0')
                  FIGinfo.HorizontalLabels.color{ii}='k';
               else
                  FIGinfo.HorizontalLabels.color{ii}=FeatureColors{-rem(FeatIND2,2)+2};
               end
               FIGinfo.HorizontalLabels.linestyle{ii}=':';
               FIGinfo.HorizontalLabels.linewidth{ii}=1;
            end
         end
         if sum(strcmp(STIMtype,{'TrBFi'}))
            ii=1;
            FeatIND2=find(strcmp(FeaturesText,'TN'));
            FIGinfo.HorizontalLabels.Text{ii}=FeaturesText{FeatIND2};
            FIGinfo.HorizontalLabels.Yvals{ii}=DFTs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}/1000;
            FIGinfo.HorizontalLabels.color{ii}='k';
            FIGinfo.HorizontalLabels.linestyle{ii}=':';
            FIGinfo.HorizontalLabels.linewidth{ii}=1;
         end

         %%%% Rate Plot
         PLOTnum=(FeatIND-1)*NUMcols+1;
         eval(['h' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Rdata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
            semilogy(CHdata.Adata(Pind),DFTs.Info.BF_kHz,'o','MarkerSize',6,'Color',PARAMcolors{Pind})
         end
         xlabel(sprintf('Rate (sp/sec)  [O: ALSR]'),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_rate)
         ylim(YLIMITS)  % Same Ylimits for all plots
         ylabel(DFTs.ChannelInfo.channel_parameter,'FontSize',FIG.LABELfontSize+2)
         text(0.05,.99,FIGinfo.CONDlabel_text,'Units','norm','Vert','top','FontSize',FIG.FEATfontSize)
         set(gca,'XDir','reverse')
         set(gca,'YTick',YTICKS,'YTickLabel',YTICKS)
         set(gca,'FontSize',FIG.ANNOTfontSize)
         
         %%%% Synch Plot
         PLOTnum=PLOTnum+1;
         eval(['h' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Sdata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
         end
         xlabel(sprintf('Synch Coef to %s',DFTs.ConditionInfo.features{FeatIND}),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_synch)
         ylim(YLIMITS)  % Same Ylimits for all plots
         if (FeatIND==1)
            title(FIGinfo.title,'FontSize',FIG.TITLEfontSize)
         end
         set(gca,'XDir','reverse')
         set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
         set(gca,'YTick',YTICKS,'YTickLabel','')
         set(gca,'FontSize',FIG.ANNOTfontSize)

         %%%% Phase Plot
         PLOTnum=PLOTnum+1;
         eval(['h' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Pdata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
         end
         xlabel(sprintf('Phase to %s (cyc)',DFTs.ConditionInfo.features{FeatIND}),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_phase)
         ylim(YLIMITS)  % Same Ylimits for all plots
         set(gca,'XDir','reverse')
         set(gca,'XTick',[-1 -1/2 0 1/2 1],'XTickLabel',[-1 -1/2 0 1/2 1])
         set(gca,'YTick',YTICKS,'YTickLabel','')
         set(gca,'FontSize',FIG.ANNOTfontSize)

         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Plot Horizontal labels/lines for feature frequencies
         for PLOToffset=0:-1:-2
            eval(['set(gcf,''CurrentAxes'',h' num2str(PLOTnum+PLOToffset) ')'])
            if PLOToffset == 0
               XLIMITS_temp=XLIMITS_phase;
            elseif PLOToffset == -1
               XLIMITS_temp=XLIMITS_synch;
            elseif PLOToffset == -2
               XLIMITS_temp=XLIMITS_rate;
            end
            for ii=1:length(FIGinfo.HorizontalLabels.Text)
               if (FIGinfo.HorizontalLabels.Yvals{ii}>YLIMITS(1))&(FIGinfo.HorizontalLabels.Yvals{ii}<YLIMITS(2))
                  if PLOToffset==0
                     text(XLIMITS_temp(1),FIGinfo.HorizontalLabels.Yvals{ii},FIGinfo.HorizontalLabels.Text{ii}, ...
                        'Units','data','Vert','Middle','Horiz','left','FontSize',FIG.HORLABfontSize,'Color',FIGinfo.HorizontalLabels.color{ii})
                  end
                  plot(XLIMITS_temp,FIGinfo.HorizontalLabels.Yvals{ii}*ones(1,2),'Color',FIGinfo.HorizontalLabels.color{ii}, ...
                     'LineStyle',FIGinfo.HorizontalLabels.linestyle{ii},'LineWidth',FIGinfo.HorizontalLabels.linewidth{ii})
               end
            end
            hold off
         end
      end
   end
end

%% Clean up axes
Xcorner=0.1;
Xwidth=.27;
Xshift=0.025;
Ycorner=0.05;
Yshift=0.05;
Ywidth=.94*(1-NumFEATURES*(Yshift+.01))/NumFEATURES;   %.26 for 3; .42 for 2
TICKlength=0.02;

for FeatIND=1:NumFEATURES
   % Rate
   eval(['HAND=h' num2str((FeatIND-1)*NUMcols+1) ';'])
   set(HAND,'Position',[Xcorner Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
   % Synch
   eval(['HAND=h' num2str((FeatIND-1)*NUMcols+2) ';'])
   set(HAND,'Position',[Xcorner+(Xwidth+Xshift) Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
   % Phase
   eval(['HAND=h' num2str((FeatIND-1)*NUMcols+3) ';'])
   set(HAND,'Position',[Xcorner+2*(Xwidth+Xshift) Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% F0 plot (with NUMspikes) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: for now, just cut and paste - could be cleaned up later since
% so much similar code!!!!  (Not so similar - leave as is)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf
set(gcf,'DefaultTextInterpreter','none')
set(gcf,'units','norm','pos',[0.46   0.07    0.45    0.85])
set(gcf,'Name',sprintf('RPSs(F0) for Unit: %s',UnitName))
set(gcf,'NumberTitle','off')

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         FIGinfo.title=sprintf('%s: [F0] Rate, Synch, Phase Plots [FormsAtHarms: %s, PolarityInverted: %s]\nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
            STIMtype,DFTs.ConditionInfo.FormsAtHarmonics{HarmIND},DFTs.ConditionInfo.polarities{PolIND}, ...
            ExpName,TrackNum,UnitNum,DFTs.Info.BF_kHz,DFTs.Info.Threshold_dBSPL,DFTs.Info.SR_sps,DFTs.Info.Q10);
         FIGinfo.CONDlabel_text=sprintf('Feature: %s',DFTs.ConditionInfo.features{FeatIND});

         CHdata.CHvals=DFTs.ChannelInfo.channel_values{FeatIND,PolIND,HarmIND};
         %% FIND BF channel
         BFind=find(abs(log2(CHdata.CHvals/DFTs.Info.BF_kHz))<BFoctCRIT);

         CHdata.Rdata=DFTs.F0_RSP.rate_sps{FeatIND,PolIND,HarmIND}; % Rate
         CHdata.Sdata=DFTs.F0_RSP.synch{FeatIND,PolIND,HarmIND}; % Synch
         CHdata.Pdata=DFTs.F0_RSP.phase_cyc{FeatIND,PolIND,HarmIND}; % Phase
         CHdata.Ndata=DFTs.DFTs.PH_orig_NumDrivenSpikes{FeatIND,PolIND,HarmIND}; % Number of Spikes

         %%%% NUMspikes Plot
         PLOTnum=(FeatIND-1)*NUMcols+1;
         eval(['hB' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Ndata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
         end
         xlabel(sprintf('Number of Driven Spikes'),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_Nspikes)
         ylim(YLIMITS)  % Same Ylimits for all plots
         ylabel(DFTs.ChannelInfo.channel_parameter,'FontSize',FIG.LABELfontSize+2)
         text(0.05,.99,FIGinfo.CONDlabel_text,'Units','norm','Vert','top','FontSize',FIG.FEATfontSize)
         set(gca,'XDir','reverse')
         set(gca,'YTick',YTICKS,'YTickLabel',YTICKS)
         set(gca,'FontSize',FIG.ANNOTfontSize)
         
         %%%% Synch Plot
         PLOTnum=PLOTnum+1;
         eval(['hB' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Sdata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
         end
         xlabel(sprintf('Synch Coef to F0'),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_synch)
         ylim(YLIMITS)  % Same Ylimits for all plots
         if (FeatIND==1)
            title(FIGinfo.title,'FontSize',FIG.TITLEfontSize)
         end
         set(gca,'XDir','reverse')
         set(gca,'XTick',[0 .25 .5 .75 1],'XTickLabel',{'0','.25','.5','.75','1'})
         set(gca,'YTick',YTICKS,'YTickLabel','')
         set(gca,'FontSize',FIG.ANNOTfontSize)

         %%%% Phase Plot
         PLOTnum=PLOTnum+1;
         eval(['hB' num2str(PLOTnum) '=subplot(NumFEATURES,NUMcols,PLOTnum);'])
         for Pind=1:NumP
            semilogy(CHdata.Pdata(:,Pind),CHdata.CHvals,'*-','Color',PARAMcolors{Pind})
            hold on
         end
         xlabel(sprintf('Phase to F0 (cyc)'),'FontSize',FIG.LABELfontSize)
         xlim(XLIMITS_phase)
         ylim(YLIMITS)  % Same Ylimits for all plots
         set(gca,'XDir','reverse')
         set(gca,'XTick',[-1 -1/2 0 1/2 1],'XTickLabel',[-1 -1/2 0 1/2 1])
         set(gca,'YTick',YTICKS,'YTickLabel','')
         set(gca,'FontSize',FIG.ANNOTfontSize)

         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Plot Horizontal labels/lines for feature frequencies
         for PLOToffset=0:-1:-2
            eval(['set(gcf,''CurrentAxes'',hB' num2str(PLOTnum+PLOToffset) ')'])
            if PLOToffset == 0
               XLIMITS_temp=XLIMITS_phase;
            elseif PLOToffset == -1
               XLIMITS_temp=XLIMITS_synch;
            elseif PLOToffset == -2
               XLIMITS_temp=XLIMITS_rate;
            end
            for ii=1:length(FIGinfo.HorizontalLabels.Text)
               if (FIGinfo.HorizontalLabels.Yvals{ii}>YLIMITS(1))&(FIGinfo.HorizontalLabels.Yvals{ii}<YLIMITS(2))
                  if PLOToffset==0
                     text(XLIMITS_temp(1),FIGinfo.HorizontalLabels.Yvals{ii},FIGinfo.HorizontalLabels.Text{ii}, ...
                        'Units','data','Vert','Middle','Horiz','left','FontSize',FIG.HORLABfontSize,'Color',FIGinfo.HorizontalLabels.color{ii})
                  end
                  plot(XLIMITS_temp,FIGinfo.HorizontalLabels.Yvals{ii}*ones(1,2),'Color',FIGinfo.HorizontalLabels.color{ii}, ...
                     'LineStyle',FIGinfo.HorizontalLabels.linestyle{ii},'LineWidth',FIGinfo.HorizontalLabels.linewidth{ii})
               end
            end
            hold off
         end
      end
   end
end

%% Clean up axes
% set(get(h2,'Title'),'String','')
Xcorner=0.1;
Xwidth=.27;
Xshift=0.025;
Ycorner=0.05;
Yshift=0.05;
Ywidth=.94*(1-NumFEATURES*(Yshift+.01))/NumFEATURES;   %.26 for 3; .42 for 2

TICKlength=0.02;

for FeatIND=1:NumFEATURES
   % Rate
   eval(['HAND=hB' num2str((FeatIND-1)*NUMcols+1) ';'])
   set(HAND,'Position',[Xcorner Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
   % Synch
   eval(['HAND=hB' num2str((FeatIND-1)*NUMcols+2) ';'])
   set(HAND,'Position',[Xcorner+(Xwidth+Xshift) Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
   % Phase
   eval(['HAND=hB' num2str((FeatIND-1)*NUMcols+3) ';'])
   set(HAND,'Position',[Xcorner+2*(Xwidth+Xshift) Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
end


return;



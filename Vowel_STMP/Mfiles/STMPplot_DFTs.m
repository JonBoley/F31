function STMPplot_DFTs(ExpDate,UnitName,STIMtype,DATAtype)
% FROM: STMPplot_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% M. Heinz Jun 10, 2007
% Plots Synchronized Rates (DFTs) for a STMP POPULATION.  Assumes STMP format for
% data storage, but this can be either: DATAtype=
%     - 1unitNstim (DFTs.Info.STMPshift=0)
%     - Nchans1stim (DFTs.Info.STMPshift=1)
%
% Setup to be general for STMP data, but some specifics for different STIM
% types.
%
% Calls basic function: neurogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
global FeaturesText

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=2;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03),
%%%% Specify ExpDate if not provided
if ~exist('ExpDate','var')
   if TESTunitNUM==1
      ExpDate='111804'; UnitName='1.28'; STIMtype='EHrBFi'; DATAtype='Nchans1stim';
   elseif TESTunitNUM==2
      ExpDate='041805'; UnitName='2.04'; STIMtype='EHvNrBFi'; DATAtype='Nchans1stim'
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
disp(sprintf('... STMP: Plotting ''%s'' PERhists  for:  Experiment: ''%s''; Unit: %d.%02d',DATAtype,ExpName,TrackNum,UnitNum))

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

FIGinfo.Ylog=1;  FIGinfo.Xlog=0;  % LOG y-axis
FIGinfo.Ylabel_text=DFTs.ChannelInfo.channel_parameter;
FIGinfo.Xlabel_text='Frequency (kHz)';
FIGinfo.param_units = DFTs.ParameterInfo.parameter_units;
FIGinfo.LineStyle='-';
FIGinfo.Marker='x';

FIGinfo.XLIMITS=[0 5];   % 5 kHz XMAX
FIGinfo.chGAIN=2;
PARAMInfo=DFTs.ParameterInfo;

BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave
CHdata.TriFiltWidth=1;
FeatureColors={'r','g'};  % 'Formants','Troughs'

global BASE10_FigureNUMS
figure(2+BASE10_FigureNUMS); clf
set(gcf,'DefaultTextInterpreter','none')
set(gcf,'units','norm','pos',[0.46   0.07    0.45    0.85])
set(gcf,'Name',sprintf('DFTs for Unit: %s',UnitName))
set(gcf,'NumberTitle','off')

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         FIGinfo.title=sprintf('%s: Synchronized Rate Plots [TriFiltWidth=%d] [FormsAtHarms: %s, PolarityInverted: %s]\nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
            STIMtype,CHdata.TriFiltWidth,DFTs.ConditionInfo.FormsAtHarmonics{HarmIND},DFTs.ConditionInfo.polarities{PolIND}, ...
            ExpName,TrackNum,UnitNum,DFTs.Info.BF_kHz,DFTs.Info.Threshold_dBSPL,DFTs.Info.SR_sps,DFTs.Info.Q10);
         FIGinfo.CONDlabel_text=sprintf('Feature: %s',DFTs.ConditionInfo.features{FeatIND});

         CHdata.CHvals=DFTs.ChannelInfo.channel_values{FeatIND,PolIND,HarmIND};
         CHdata.Ydata=DFTs.DFTs.DFTs{FeatIND,PolIND,HarmIND};
         % Take magnitude of DFT to get Synchronized_Rate_sps
         for ChanIND=1:NumCH
            for ParamIND=1:NumP
               CHdata.Ydata{ChanIND,ParamIND}=abs(CHdata.Ydata{ChanIND,ParamIND});
            end
         end
         CHdata.Xdata=cell(size(CHdata.Ydata));  % Frequency
         %% FIND BF channel to make BOLD
         CHdata.BOLDchs=find(abs(log2(CHdata.CHvals/DFTs.Info.BF_kHz))<BFoctCRIT);

         % Setup frequency_kHz vector
         for ChanIND=1:NumCH
            for ParamIND=1:NumP
               CHdata.Xdata{ChanIND,ParamIND}=(1/1000)*(0:length(CHdata.Ydata{ChanIND,ParamIND})-1) * DFTs.DFTs.DFT_freqbinwidth_Hz{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND);
            end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Setup Vertical labels for temporal periods of features
         if sum(strcmp(STIMtype,{'EHrBFi','EHvNrBFi'}))
            ii=0;
            for FeatIND2=1:find(strcmp(FeaturesText,'T3'))
               ii=ii+1;
               FIGinfo.VerticalLabels.Text{ii}=FeaturesText{FeatIND2};
               FIGinfo.VerticalLabels.Xvals{ii}=DFTs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}(FeatIND2)/1000;
               if strcmp(FeaturesText{FeatIND2},'F0')
                  FIGinfo.VerticalLabels.color{ii}='k';
               else
                  FIGinfo.VerticalLabels.color{ii}=FeatureColors{-rem(FeatIND2,2)+2};
               end
               FIGinfo.VerticalLabels.linestyle{ii}=':';
               FIGinfo.VerticalLabels.linewidth{ii}=2.5;
            end
         end

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
                  FIGinfo.HorizontalLabels.linestyle{ii}='--';
                  FIGinfo.HorizontalLabels.linewidth{ii}=1;
            end
         end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Setup Other lines
         ii=1;
         FIGinfo.OtherLines.Xvals{ii}=CHdata.CHvals;
         FIGinfo.OtherLines.Yvals{ii}=CHdata.CHvals;
         FIGinfo.OtherLines.color{ii}='k';
         FIGinfo.OtherLines.linestyle{ii}='-';
         FIGinfo.OtherLines.linewidth{ii}=1;

         eval(['h'  num2str(FeatIND) '=subplot(NumFEATURES,1,FeatIND);'])         
         neurogram(CHdata,FIGinfo,PARAMInfo)

         if FeatIND~=1
            legend off
            title('')
         end
         if FeatIND~=NumFEATURES
            xlabel('')
         end
      end
   end
end

%% Clean up axes
% set(get(h2,'Title'),'String','')
Xcorner=0.1;
Xwidth=.85;
Ycorner=0.05;
Yshift=0.05;
Ywidth=.94*(1-NumFEATURES*(Yshift+.01))/NumFEATURES;   %.26 for 3; .42 for 2

TICKlength=0.02;

for FeatIND=1:NumFEATURES
   eval(['HAND=h' num2str(FeatIND) ';'])
   set(HAND,'Position',[Xcorner Ycorner+(NumFEATURES-FeatIND)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
end

return;

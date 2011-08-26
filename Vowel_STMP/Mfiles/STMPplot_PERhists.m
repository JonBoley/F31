function STMPplot_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% FROM: STMPcompute_PERhists(ExpDate,UnitName,STIMtype,DATAtype)
% M. Heinz Jun 09, 2007
% Plots PERhistograms for a STMP POPULATION.  Assumes STMP format for
% data storage, but this can be either: DATAtype=
%     - 1unitNstim (PHs.Info.STMPshift=0)
%     - Nchans1stim (PHs.Info.STMPshift=1)
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
      ExpDate='041805'; UnitName='2.04'; STIMtype='EHvNrBFi'; DATAtype='Nchans1stim';
   elseif TESTunitNUM==3
      ExpDate='111606'; UnitName='2.09';
   elseif TESTunitNUM==4
      ExpDate='041307'; UnitName='1.03';
   elseif TESTunitNUM==5
      ExpDate='041805'; UnitName='2.04'; STIMtype='TrBFi'; DATAtype='Nchans1stim';  % TONE re BFI, with Tone RLV
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
%%%% Verify PHs exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PHs_filename=sprintf('PHs_%s.%d.%02d.%s.mat',DATAtype,TrackNum,UnitNum,STIMtype);
PHs_dataname=sprintf('PHs_%s',DATAtype);
eval(['load ''' fullfile(STMPanal_dir,PHs_filename) ''''])
eval(['PHs = ' PHs_dataname ';']);
eval(['clear ' PHs_dataname]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this data, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(PHs.PERhists.PERhists{1,1,1},1);
NumP=size(PHs.PERhists.PERhists{1,1,1},2);
NumFEATURES=size(PHs.PERhists.PERhists,1);
NumPOL=size(PHs.PERhists.PERhists,2);
NumHARMS=size(PHs.PERhists.PERhists,3);

%% Find min F0 across all features - to set PERhist XLIMITS
F0min=Inf;
for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         for ChanIND=1:NumCH
            for ParamIND=1:NumP
               F0min=min([F0min PHs.ChannelInfo.channel_F0_Hz{FeatIND,PolIND,HarmIND}(ChanIND)]);
            end  % end Param
         end  % end Freq
      end % End FormsatHarms
   end % End Invert Polarity
end % End: FeatIND
if strcmp(STIMtype,'TrBFi')
   NUM_F0_cycles=5;
else
   NUM_F0_cycles=1;
end
PERhist_XMAX=NUM_F0_cycles*1/F0min*1000;  % take the PERhist_Xlims on vowels (in ms)

FIGinfo.Ylog=1;  FIGinfo.Xlog=0;  % LOG y-axis
FIGinfo.Ylabel_text=PHs.ChannelInfo.channel_parameter;
FIGinfo.Xlabel_text='Time (ms)';
FIGinfo.param_units = PHs.ParameterInfo.parameter_units;

FIGinfo.XLIMITS=[0 PERhist_XMAX];
FIGinfo.chGAIN=2;
PARAMInfo=PHs.ParameterInfo;

BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave
CHdata.TriFiltWidth=7;
FeatureColors={'r','g'};  % 'Formants','Troughs'

global BASE10_FigureNUMS

if isempty(BASE10_FigureNUMS)
	BASE10_FigureNUMS=2000;
end

figure(1+BASE10_FigureNUMS); clf
set(gcf,'DefaultTextInterpreter','none')
set(gcf,'units','norm','pos',[0.005    0.07    0.45    0.85])
set(gcf,'Name',sprintf('PERhists for Unit: %s',UnitName))
set(gcf,'NumberTitle','off')

for FeatIND=1:NumFEATURES
   for PolIND=1:NumPOL
      for HarmIND=1:NumHARMS
         FIGinfo.title=sprintf('%s: Period Histograms [TriFiltWidth=%d] [FormsAtHarms: %s, PolarityInverted: %s]\nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
            STIMtype,CHdata.TriFiltWidth,PHs.ConditionInfo.FormsAtHarmonics{HarmIND},PHs.ConditionInfo.polarities{PolIND}, ...
            ExpName,TrackNum,UnitNum,PHs.Info.BF_kHz,PHs.Info.Threshold_dBSPL,PHs.Info.SR_sps,PHs.Info.Q10);
         FIGinfo.CONDlabel_text=sprintf('Feature: %s',PHs.ConditionInfo.features{FeatIND});

         CHdata.CHvals=PHs.ChannelInfo.channel_values{FeatIND,PolIND,HarmIND};
         CHdata.Ydata=PHs.PERhists.PERhists{FeatIND,PolIND,HarmIND};
         % Show 'NUM_F0_cycles' of stim
         if NUM_F0_cycles>1
            for ChanIND=1:NumCH
               for ParamIND=1:NumP
                  CHdata.Ydata{ChanIND,ParamIND}=repmat(CHdata.Ydata{ChanIND,ParamIND},1,NUM_F0_cycles);
               end
            end
         end
         
         CHdata.Xdata=cell(size(CHdata.Ydata));
         %% FIND BF channel to make BOLD
         CHdata.BOLDchs=find(abs(log2(CHdata.CHvals/PHs.Info.BF_kHz))<BFoctCRIT);
         
         % Setup time_mesc vector
         for ChanIND=1:NumCH
            for ParamIND=1:NumP
               CHdata.Xdata{ChanIND,ParamIND}=1000*(0.5:length(CHdata.Ydata{ChanIND,ParamIND})) * PHs.PERhists.PH_binwidth_sec{FeatIND,PolIND,HarmIND}(ChanIND,ParamIND);
            end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Setup Vertical labels for temporal periods of features
         if sum(strcmp(STIMtype,{'EHrBFi','EHvNrBFi'}))
            ii=0;
            for FeatIND2=1:find(strcmp(FeaturesText,'T3'))
               if 1/PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}(FeatIND2) > 0.5/1000
                  ii=ii+1;
                  FIGinfo.VerticalLabels.Text{ii}=['1/' FeaturesText{FeatIND2}];  %% 1/
                  FIGinfo.VerticalLabels.Xvals{ii}=1/PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}(FeatIND2)*1000;
                  if strcmp(FeaturesText{FeatIND2},'F0')
                     FIGinfo.VerticalLabels.color{ii}='k';
                  else
                     FIGinfo.VerticalLabels.color{ii}=FeatureColors{-rem(FeatIND2,2)+2};
                  end
                  FIGinfo.VerticalLabels.linestyle{ii}='';
                  FIGinfo.VerticalLabels.linewidth{ii}=1;
               end
            end
         end
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% Setup Horizontal labels/lines for feature frequencies
         if sum(strcmp(STIMtype,{'EHrBFi','EHvNrBFi'}))
            ii=0;
            for FeatIND2=1:find(strcmp(FeaturesText,'T3'))
                  ii=ii+1;
                  FIGinfo.HorizontalLabels.Text{ii}=FeaturesText{FeatIND2};
                  FIGinfo.HorizontalLabels.Yvals{ii}=PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}(FeatIND2)/1000;
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
            FIGinfo.HorizontalLabels.Yvals{ii}=PHs.ChannelInfo.channel_FeatureFreqs_Hz{FeatIND,PolIND,HarmIND}{1}/1000;
            FIGinfo.HorizontalLabels.color{ii}='k';
            FIGinfo.HorizontalLabels.linestyle{ii}=':';
            FIGinfo.HorizontalLabels.linewidth{ii}=1;
         end

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

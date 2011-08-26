function STMPplot_RLFs(ExpDate,UnitName)
% FROM: STMPsave_STs_1unitNstim_EHrBFi(ExpDate,UnitName,STIMtype)
% M. Heinz Jun 11, 2007
% Plots all RLFs for 1 unit for STMP stimuli
%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% STIMtype Possibilities
%
% NOTE - all RLVs are shown (TB, EHrlv, EHINrlv, SACrlvQ)
% labeling depends on STMP data STIMtypes available
%
% *1) EHrlv alone (no levels labeled) [111804, 1.33]
% *2) EHrlv and EHINrlv (Vowel level labeled from EHINrlv) [041805, 2.05]
% *3) EHrlv and EHrBFi (show EH levels in EHrBFI from STs_1unitNstim)
% *4) EHrlv, EHINrlv, and EHvNrBFi (Vowel level labeled from EHINrlv; show Noise Attens in EHvNrBFI from STs_1unitNstim)
% *5) TrBFI and ToneRLV[TB] (show TONE levels in TrBFI from STs_1unitNstim)
% 6) SACrlv
%
% TODO LATER:  
% *1) Add - ToneRLV w/&wo/ TrBFi  {ISH: 041804, 2.04 - has
% both]***************
% 2) Add - SACrlvQ w/&wo/ SACrlv
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global STMP_dir STMP_ExpList
global FeaturesText FormsAtHarmonicsText InvertPolarityText

% close previous RLF plots
OPENfigs=get(0,'Children');
close(intersect([100 101],OPENfigs))

YLIMITS_ALL=[0 300];  % SAME YLIMITS FOR ALL RLF PLOTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TESTunitNUM=3;  % 1: ARO (111804, 1.28), 2: ISH (041805, 2.04), 3: Purdue1 (111606, 2.09); 4: Purdue2 (041306, 1.03),
%%%% Specify ExpDate if not provided
if ~exist('ExpDate','var')
   if TESTunitNUM==1
      ExpDate='111804'; UnitName='1.28';  % EHrlv w/EHrBFi [ARO2005]
   elseif TESTunitNUM==2
      ExpDate='111804'; UnitName='1.33'; % EHrlv alone
   elseif TESTunitNUM==3
      ExpDate='041805'; UnitName='2.04'; % EHrlv, EHINrlv, w/EHvNrBFi [ISH2006]
                                         % ALSO TB with TrBFi
                                         % ALSO SACrlv
   elseif TESTunitNUM==4
      ExpDate='041805'; UnitName='2.05'; % EHrlv, EHINrlv alone
   elseif TESTunitNUM==5
      ExpDate='111606'; UnitName='2.09'; % 
   elseif TESTunitNUM==6
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
RAWdata_dir=fullfile(STMP_dir,'ExpData',ExpName);
eval(['cd ''' RAWdata_dir ''''])
UNITinfo_dir=fullfile(STMP_dir,'ExpData',ExpName,'UNITinfo');   % For general unit info (BF, SR, bad lines, ...)
if ~exist('UNITinfo','dir')
   mkdir('UNITinfo');
end
STMPanal_dir=fullfile(STMP_dir,'ExpData',ExpName,'STMPanalyses');   % For STMP analyses (Spike Trains, PSTs, PerHist, DFTs, SAC/SCCs, ...)
if ~exist('UNITinfo','dir')
   mkdir('UNITinfo');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify unitINFO exists, if not create
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unitINFO_filename=sprintf('unitINFO.%d.%02d.mat',TrackNum,UnitNum);
eval(['ddd=dir(''' fullfile(UNITinfo_dir,unitINFO_filename) ''');'])
if isempty(ddd)
   eval(['cd ''' RAWdata_dir ''''])
   STMPsave_UnitINFO(TrackNum,UnitNum,1);  % skip verify
end
eval(['load ''' fullfile(UNITinfo_dir,unitINFO_filename) ''''])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Decide STIMtype based on what STMP data is available
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('STIMtype','var')
   TrBFi_PicList=findPics('TrBFi',[TrackNum UnitNum]);
   EHrBFi_PicList=findPics('EHrBFi',[TrackNum UnitNum]);
   EHvNrBFi_PicList=findPics('EHvNrBFi',[TrackNum UnitNum]);
   SACrlv_PicList=findPics('SACrlv',[TrackNum UnitNum]);
   STIMtypes_TEXT={'TrBFi','EHrBFi','EHvNrBFi','SACrlv'};
   STMPdataYES=[~isempty(TrBFi_PicList) ~isempty(EHrBFi_PicList) ~isempty(EHvNrBFi_PicList) ~isempty(SACrlv_PicList)];
   STIMtype=STIMtypes_TEXT(find(STMPdataYES));
end

CalibPic=unitINFO.Info.CalibPIC_ToUse;
TextFontSize=12;
AnnotFonotSize=16;
LineWidth=3;
TITLEFontSize=10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show EHrlfs for this unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EHrlf_piclist=findPics('EHrlv',[TrackNum UnitNum]);
if ~isempty(EHrlf_piclist)
   disp(sprintf('... STMP: Plotting "EHrlfs" for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))
   plot_EHrlfs(EHrlf_piclist,CalibPic)

   EHrlf_FIG=gcf;
   set(EHrlf_FIG,'Name',sprintf('EHrlfs for Unit: %s',UnitName))
   set(EHrlf_FIG,'NumberTitle','off')
   ylim(YLIMITS_ALL);
   YLIMITS_EHrlf=ylim;
   XLIMITS_EHrlf=xlim;
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f', ...
      ExpDate,UnitName,unitINFO.Info.BF_kHz,unitINFO.Info.Threshold_dBSPL,unitINFO.Info.SR_sps,unitINFO.Info.Q10), ...
      'units','norm','FontSize',TITLEFontSize)

   
   % Label L_EH - vowel level from EHINrlv if it exists, otherwise from EHrBFI
   % it is exists
   VOWELlevels_dBSPL=[];
   EHINrlf_piclist=findPics('EHINrlv',[TrackNum UnitNum]);
   if ~isempty(EHINrlf_piclist)
      x=loadPic(EHINrlf_piclist(1));
      VOWELlevels_dBSPL=x.Stimuli.Condition.Level_dBSPL;    
   elseif sum(strcmp(STIMtype,'EHrBFi'))
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% Load STs_1unitNstim if it exists
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      STs_1unitNstim_filename=sprintf('STs_1unitNstim.%d.%02d.%s.mat',TrackNum,UnitNum,'EHrBFi');
      eval(['ddd=dir(''' fullfile(STMPanal_dir,STs_1unitNstim_filename) ''');'])
      if ~isempty(ddd)
         eval(['load ''' fullfile(STMPanal_dir,STs_1unitNstim_filename) ''''])
         VOWELlevels_dBSPL=STs_1unitNstim.ParameterInfo.param_values;
      else
         beep
         disp(sprintf('File: %s HAS NOT BEEN CREATED YET - unable to label vowel levels',STs_1unitNstim_filename))
      end
   end
   for ii=1:length(VOWELlevels_dBSPL)
      plot(VOWELlevels_dBSPL(ii)*ones(2,1),[YLIMITS_EHrlf],'k:','LineWidth',LineWidth)
      text(VOWELlevels_dBSPL(ii)+0.01*diff(XLIMITS_EHrlf),50-(ii-1)*15,sprintf('%.f',VOWELlevels_dBSPL(ii)),'HorizontalAlignment','left','FontSize',TextFontSize,'units','data')
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show EHINrlfs for this unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EHrlf_piclist=findPics('EHINrlv',[TrackNum UnitNum]);
if ~isempty(EHrlf_piclist)
   disp(sprintf('... STMP: Plotting "EHINrlfs" for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))
   dBAtt_2_SNR = plot_EHINrlfs(EHINrlf_piclist,CalibPic);

   EHINrlf_FIG=gcf;
   set(EHINrlf_FIG,'Name',sprintf('EHINrlfs for Unit: %s',UnitName))
   set(EHINrlf_FIG,'NumberTitle','off')
   ylim(YLIMITS_ALL);
   YLIMITS_EHINrlf=ylim;
   XLIMITS_EHINrlf=xlim;
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f', ...
      ExpDate,UnitName,unitINFO.Info.BF_kHz,unitINFO.Info.Threshold_dBSPL,unitINFO.Info.SR_sps,unitINFO.Info.Q10), ...
      'units','norm','FontSize',TITLEFontSize)

   % Label A_N - noise Atten for rest of data
   EHvNreBFi_piclist=findPics('EHvNrBFi',[TrackNum UnitNum]);
   if ~isempty(EHvNreBFi_piclist)
      x=loadPic(EHvNreBFi_piclist(1));
      dBattns=x.Stimuli.Condition.NoiseAttens_dB;
      dBattns=dBattns(find(dBattns~=120));  % Take out 120 dB - In Quiet
      SNRs_dB=dBattns+dBAtt_2_SNR;
      for ii=1:length(SNRs_dB)
         plot(SNRs_dB(ii)*ones(2,1),[0 1000],'k:','LineWidth',LineWidth)
         text(SNRs_dB(ii)-0.01*diff(XLIMITS_EHINrlf),50,sprintf('%.f',SNRs_dB(ii)),'HorizontalAlignment','left','units','data','FontSize',TextFontSize)
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show TONErlfs(TB) for this unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TONErlf_piclist=findPics('TB',[TrackNum UnitNum]);
if ~isempty(TONErlf_piclist)
   disp(sprintf('... STMP: Plotting "TONErlfs" for:  Experiment: ''%s''; Unit: %d.%02d',ExpName,TrackNum,UnitNum))
   plot_EHrlfs(TONErlf_piclist,CalibPic,102)  % Use FIGnum=102 to avoid 100: EHrlf, 101: EHINrlf

   TONErlf_FIG=gcf;
   set(TONErlf_FIG,'Name',sprintf('TONErlfs for Unit: %s',UnitName))
   set(TONErlf_FIG,'NumberTitle','off')
   ylim(YLIMITS_ALL);
   YLIMITS_TONErlf=ylim;
   XLIMITS_TONErlf=xlim;
   title(sprintf('     Exp%s, Unit %s: BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f', ...
      ExpDate,UnitName,unitINFO.Info.BF_kHz,unitINFO.Info.Threshold_dBSPL,unitINFO.Info.SR_sps,unitINFO.Info.Q10), ...
      'units','norm','FontSize',TITLEFontSize)

   % Label Tone levels for STMP data
   TreBFi_piclist=findPics('TrBFi',[TrackNum UnitNum]);
   if ~isempty(TreBFi_piclist)
      x=loadPic(TreBFi_piclist(1));
      TONElevels_dBSPL=x.Stimuli.Condition.Levels_dBSPL;
      for ii=1:length(TONElevels_dBSPL)
         plot(TONElevels_dBSPL(ii)*ones(2,1),[0 1000],'k:','LineWidth',LineWidth)
         text(TONElevels_dBSPL(ii)+0.01*diff(XLIMITS_TONErlf),50,sprintf('%.f',TONElevels_dBSPL(ii)),'HorizontalAlignment','left','units','data','FontSize',TextFontSize)
      end
   end
   set(TONErlf_FIG,'units','norm','pos',[0.2    0.5324    0.4000    0.39])
end

%%% LATER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Show SACrlfs for this unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

function STMPplot_SCCfunctions(ExpDate,UnitName,FeatureIndices,AttenIndices)
% STMPplot_SCCfunctions(ExpDate,UnitName)
% Plots Cross-Correlation Functions for a STMP POPULATION.
% Modified specifically for vowel STMP in noise (EHIN) data
% 
% J. Boley (April 12, 2012)
% FROM: STMPplot_SCCfunctions (M. Heinz May 28, 2009)
% 
% Calls basic function: neurogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ExpList ROOT_dir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('ExpDate','var')
    ExpDate='062311'; UnitName='1.01'; FeatureIndices=[3]; AttenIndices=[1];
end

%%%% Find the full Experiment Name
ExpDateText=strcat('20',ExpDate(end-1:end),'_',ExpDate(1:2),'_',ExpDate(3:4));
for i=1:length(ExpList)
    if ~isempty(strfind(ExpList{i},ExpDateText))
        ExpName=ExpList{i};
        break;
    end
end
if ~exist('ExpName','var')
    disp(sprintf('***ERROR***:  Experiment: %s not found\n   Experiment List:',ExpDate))
    disp(strvcat(ExpList))
    beep
    error('STOPPED');
end

%%%% Parse out the Track and Unit Number
TrackNum=str2num(UnitName(1:strfind(UnitName,'.')-1));
UnitNum=str2num(UnitName(strfind(UnitName,'.')+1:end));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Setup related directories
SACSCCanal_dir=fullfile(ROOT_dir,'ExpData',ExpName,'UNITSdata');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load SCCs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCCs_filename=sprintf('UnitLook_EHIN.%d.%02d.mat',TrackNum,UnitNum);
SCC_datanames = {'unit','NSCCs','NSACs','NSCC_delays_usec','NSAC_BFs_kHz',...
    'NSCC_BFs_kHz','NSCC_peaks','NSCC_CDs_usec','NSCC_BFinds','Nattens_dB'};

if ~exist(fullfile(SACSCCanal_dir,SCCs_filename),'file')
    error('File does not exist [%s]',fullfile(SACSCCanal_dir,SCCs_filename));
end
load(fullfile(SACSCCanal_dir,SCCs_filename),SCC_datanames{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this data, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCH=size(NSAC_BFs_kHz{1,1},2); % number of freq channels
NumFEATURES=length(FeatureIndices);
NumP=length(AttenIndices)*2; % Noise levels + hand-picked peaks

FIGinfo.Ylog=1;  FIGinfo.Xlog=0;  % LOG y-axis
FIGinfo.Ylabel_text='BF (kHz)';
FIGinfo.Xlabel_text='Delay (usec)';
FIGinfo.param_units = 'usec';
FIGinfo.LineStyle='-';
FIGinfo.Marker='none';

% FIGinfo.XLIMITS=[-1 1]*SCCs.paramsOUT{1,1,1}{1,1}.MAXdelay_sec*1e6;   %
FIGinfo.XLIMITS=[-1 1]*6000;   %
FIGinfo.chGAIN=.8;
PARAMInfo.param_values = reshape(repmat(Nattens_dB,2,1),1,2*length(Nattens_dB));

BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave
CHdata.TriFiltWidth=1;

figure(651); clf
set(gcf,'DefaultTextInterpreter','none')
if ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'units','norm','pos',[0.46   0.07    0.45    0.85])
end
set(gcf,'Name',sprintf('SCCs for Unit: %s',UnitName))
set(gcf,'NumberTitle','off')

FIGinfo.title=sprintf('SCC Plots [TriFiltWidth=%d] \nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
    CHdata.TriFiltWidth, ...
    ExpName,TrackNum,UnitNum,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10);

for FeatNUM=1:NumFEATURES 
    FeatIND=FeatureIndices(FeatNUM);
    switch FeatIND
        case 1
            [HarmIND,PolIND]=find(~cellfun(@isempty,unit.EHvLTASS_reBF_simFF.F1));
            FIGinfo.CONDlabel_text=sprintf('Feature: %s','F1');
        case 3
            [HarmIND,PolIND]=find(~cellfun(@isempty,unit.EHvLTASS_reBF_simFF.F2));
            FIGinfo.CONDlabel_text=sprintf('Feature: %s','F2');
    end


    [SCCpos,SCCs_belowBF,SCCs_aboveBF,BFind,BF]=...
        getSCCindsWithBF(unit.Info.BF_kHz,NSAC_BFs_kHz{1,1},NSCC_BFs_kHz{1,1});

    % CHvals = BFs
    for i=1:size(SCCpos,1)+1
        if i==1
            CHdata.CHvals(i) = NSAC_BFs_kHz{FeatIND,1}{BFind};
        else
            CHdata.CHvals(i) = NSCC_BFs_kHz{FeatIND,1}{SCCpos(i-1,1)}(mod(SCCpos(i-1,2),2)+1);
        end
    end
    % Find center BF (to make BOLD)... should be ch.1
    CHdata.BFref_ind=find(abs(log2(CHdata.CHvals/unit.Info.BF_kHz))<BFoctCRIT);
    CHdata.BOLDchs=CHdata.BFref_ind;

    for ChanIND=1:NumCH
        for ParamIND=1:NumP
            if mod(ParamIND,2) % odd params = full SCC sequence
                if CHdata.CHvals(ChanIND)<CHdata.CHvals(CHdata.BFref_ind)
                    % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
                    % but we want other way around, so need to flip 1st set here
                    CHdata.Ydata{ChanIND,ParamIND}= ...
                        fliplr(NSCCs{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)});
                    CHdata.Xdata{ChanIND,ParamIND}= ...
                        NSCC_delays_usec{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                elseif ChanIND==CHdata.BFref_ind % should be ch.1
                    CHdata.Ydata{ChanIND,ParamIND}= ...
                        NSACs{FeatIND,round(ParamIND/2)}{BFind};
                    CHdata.Xdata{ChanIND,ParamIND}= ...
                        NSCC_delays_usec{FeatIND,round(ParamIND/2)}{BFind};
                elseif CHdata.CHvals(ChanIND)>CHdata.CHvals(CHdata.BFref_ind)
                    CHdata.Ydata{ChanIND,ParamIND}= ...
                        NSCCs{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                    CHdata.Xdata{ChanIND,ParamIND}= ...
                        NSCC_delays_usec{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                end
            else  % even params = hand-picked CDs
                if CHdata.CHvals(ChanIND)<CHdata.CHvals(CHdata.BFref_ind)
                    % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
                    % but we want other way around, so need to flip 1st set here
                    CHdata.Ydata{ChanIND,ParamIND}= ...
                        NSCC_peaks{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                    CHdata.Xdata{ChanIND,ParamIND}= ...
                        -NSCC_CDs_usec{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                elseif ChanIND==CHdata.BFref_ind % should be ch.1
                    CHdata.Ydata{ChanIND,ParamIND}= NaN;
                    CHdata.Xdata{ChanIND,ParamIND}= NaN;
                elseif CHdata.CHvals(ChanIND)>CHdata.CHvals(CHdata.BFref_ind)
                    CHdata.Ydata{ChanIND,ParamIND}= ...
                        NSCC_peaks{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                    CHdata.Xdata{ChanIND,ParamIND}= ...
                        NSCC_CDs_usec{FeatIND,round(ParamIND/2)}{SCCpos(ChanIND-1,1)};
                end
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup Vertical labels for temporal periods of features
    ii=1;
    FIGinfo.VerticalLabels.Text{ii}='';
    FIGinfo.VerticalLabels.Xvals{ii}=0;
    FIGinfo.VerticalLabels.color{ii}='k';
    FIGinfo.VerticalLabels.linestyle{ii}='-';
    FIGinfo.VerticalLabels.linewidth{ii}=1;

    % 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 	%% Setup Horizontal labels/lines for feature frequencies
    % 	for ii=1:NumCH
    % 		FIGinfo.HorizontalLabels.Text{ii}='';
    % 		FIGinfo.HorizontalLabels.Yvals{ii}=SCCs.ChannelInfo.channel_values{FeatIND,HarmIND}(ii);
    % 		FIGinfo.HorizontalLabels.color{ii}='k';
    % 		FIGinfo.HorizontalLabels.linestyle{ii}='-';
    % 		FIGinfo.HorizontalLabels.linewidth{ii}=1;
    % 	end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup Other lines
    ii=1;
    FIGinfo.OtherLines.Xvals{ii}=CHdata.CHvals;
    FIGinfo.OtherLines.Yvals{ii}=CHdata.CHvals;
    FIGinfo.OtherLines.Yvals{ii}(end)=1;
    FIGinfo.OtherLines.Yvals{ii}(1)=2;
    FIGinfo.OtherLines.color{ii}='k';
    FIGinfo.OtherLines.linestyle{ii}='-';
    FIGinfo.OtherLines.linewidth{ii}=2;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %% Setup Other markers
    % ii=1;
    % FIGinfo.OtherMarkers.Xvals{ii}=0;
    % FIGinfo.OtherMarkers.Yvals{ii}=1;  % SACSCCmetrics.CDscc_usec/1000,SACSCCfunctions.SCC_AB_avg(find(SACSCCfunctions.delays_usec==SACSCCmetrics.CDscc_usec))
    % NormFact=(10^(chGAIN*ch_logCHwidth)-1)*CHdata.CHvals(CHind)/Y_MAXval;
    % semilogy(CHdata.Xdata{CHind,Pind}, ...
    % 	trifilt(CHdata.Ydata{CHind,Pind},CHdata.TriFiltWidth)*NormFact+CHdata.CHvals(CHind), ...
    % 	'LineWidth',LINEwidth,'LineStyle',FIGinfo.LineStyle,'Marker',FIGinfo.Marker,'Color',PARAMcolors{Pind})
    % FIGinfo.OtherMarkers.color{ii}='r';
    % FIGinfo.OtherMarkers.marker{ii}='x';
    % FIGinfo.OtherMarkers.markerSize{ii}=14;
    % FIGinfo.OtherMarkers.linewidth{ii}=2;



    eval(['h'  num2str(FeatNUM) '=subplot(NumFEATURES,1,FeatNUM);'])
    neurogram(CHdata,FIGinfo,PARAMInfo)

    legend off
end % for FeatNUM=1:NumFEATURES


%% Clean up axes
% set(get(h2,'Title'),'String','')
Xcorner=0.1;
Xwidth=.85;
Ycorner=0.15;
Yshift=0.2;
Ywidth=.7;
Ywidth=.94*(1-NumFEATURES*(Yshift+.01))/NumFEATURES;   %.26 for 3; .42 for 2

YLIMtemp=ylim;
% ylim([1 YLIMtemp(2)])

TICKlength=0.02;
% for FeatNUM=1:NumFEATURES
% 	eval(['HAND=h' num2str(FeatNUM) ';'])
% 	set(HAND,'Position',[Xcorner Ycorner+(NumFEATURES-FeatNUM)*(Ywidth+Yshift) Xwidth Ywidth],'TickLength',[TICKlength 0.025])
% % 	set(HAND,'Position',[Xcorner Ycorner Xwidth Ywidth],'TickLength',[TICKlength 0.025])
% end

set(gca,'LineWidth',2)
set(gca,'Box','on')

cd(ROOT_dir)
print('-dpng',fullfile(ROOT_dir,'ExpData',sprintf('SCCreBF_%s_%s.png',ExpDate,UnitName)));

return;

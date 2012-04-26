function varargout = calcCDreBF_manual(varargin)
% CALCCDREBF_MANUAL M-file for calcCDreBF_manual.fig
%      CALCCDREBF_MANUAL, by itself, creates a new CALCCDREBF_MANUAL or raises the existing
%      singleton*.
%
%      H = CALCCDREBF_MANUAL returns the handle to a new CALCCDREBF_MANUAL or the handle to
%      the existing singleton*.
%
%      CALCCDREBF_MANUAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALCCDREBF_MANUAL.M with the given input arguments.
%
%      CALCCDREBF_MANUAL('Property','Value',...) creates a new CALCCDREBF_MANUAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before calcCDreBF_manual_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to calcCDreBF_manual_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help calcCDreBF_manual

% Last Modified by GUIDE v2.5 26-Apr-2012 08:34:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @calcCDreBF_manual_OpeningFcn, ...
                   'gui_OutputFcn',  @calcCDreBF_manual_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before calcCDreBF_manual is made visible.
function calcCDreBF_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to calcCDreBF_manual (see VARARGIN)

% Choose default command line output for calcCDreBF_manual
handles.output = hObject;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Verify parameters and experiment, unit are valid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(varargin,2)<2, ExpDate='062311'; UnitName='1.01'; end
if size(varargin,2)<4, 
    FeatureIndices=[3]; AttenIndices=[1]; 
else
    ExpDate=varargin{1}; 
    UnitName=varargin{2}; 
    FeatureIndices=varargin{3}; 
    AttenIndices=varargin{4};
end

handles.data.ExpDate=ExpDate;
handles.data.UnitName=UnitName;
handles.data.FeatureIndices=FeatureIndices;
handles.data.AttenIndices=AttenIndices;
guidata(hObject, handles);

updateData(hObject, handles);

% UIWAIT makes calcCDreBF_manual wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = calcCDreBF_manual_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on selection change in CFlistbox.
function CFlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to CFlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns CFlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from CFlistbox

set(handles.savedText,'Visible','off');

handles = updateData(hObject, handles);

contents = get(hObject,'String');
handles.data.selectedCFkHz = str2num(contents{get(hObject,'Value')});

set(handles.detailsText,'String',sprintf('CF = %1.3f kHz\nPlease Select CD',...
    handles.data.selectedCFkHz));

% RE-COLOR THE SELECTED CF ?
% get the SCC data for this curve
CFindexTEMP = find(round(1e4*handles.data.CHdata.CHvals)/1e4==handles.data.selectedCFkHz);
handles.data.selectedCFind = CFindexTEMP;
corrTEMP = handles.data.CHdata.Ydata{CFindexTEMP,1};
delaysTEMP_usec = handles.data.CHdata.Xdata{CFindexTEMP,1};

% get user input
set(handles.instructionsText,'Visible','on');
% set(handles.CFlistbox,'Enable','off');
[x,y,button] = ginput(1);
% set(handles.CFlistbox,'Enable','on');
set(handles.instructionsText,'Visible','off');
if ~isempty(x) && button<2
    CD=x;
    %%% Need to select peak for the appropriate curve
    PEAK=corrTEMP(find(delaysTEMP_usec>=CD,1));
    %%%
     if ((x<handles.data.FIGinfo.XLIMITS(1)) || (x>handles.data.FIGinfo.XLIMITS(end)))
        CD=NaN;
        PEAK=NaN;
    end

    % update axes (after updating data)
    handles.data.CHdata.Xdata{CFindexTEMP,2}=CD;
    handles.data.CHdata.Ydata{CFindexTEMP,2}=PEAK;
    cla(handles.axes1,'reset');
    neurogram(handles.data.CHdata,handles.data.FIGinfo,handles.data.PARAMInfo);
    legend off

    % update detailsText
    set(handles.detailsText,'String',sprintf('CF = %1.3f kHz\nCD = %1.0f us',...
        handles.data.selectedCFkHz,CD));

    % SAVE CD & PEAK LOCALLY
    handles.data.selectedCD = CD;
    handles.data.selectedPeak = PEAK;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function CFlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CFlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in savebutton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE CD & PEAK TO FILE (ALONG WITH CONFIDENCE)
load(handles.data.SCCs_filename,'NSCC_peaks','NSCC_CDs_usec');

ChanIND = handles.data.selectedCFind;
FeatIND = handles.data.FeatureIndices;
ParamIND = handles.data.AttenIndices;

if handles.data.CHdata.CHvals(ChanIND)<handles.data.CHdata.CHvals(handles.data.CHdata.BFref_ind)
    % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
    % but we want other way around, so need to flip 1st set here
    NSCC_peaks{FeatIND,round(ParamIND/2)}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedPeak;
    NSCC_CDs_usec{FeatIND,round(ParamIND/2)}{handles.data.SCCpos(ChanIND-1,1)}=...
        -handles.data.selectedCD;
elseif handles.data.CHdata.CHvals(ChanIND)>handles.data.CHdata.CHvals(handles.data.CHdata.BFref_ind)
    NSCC_peaks{FeatIND,round(ParamIND/2)}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedPeak;
    NSCC_CDs_usec{FeatIND,round(ParamIND/2)}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedCD;
end

save(handles.data.SCCs_filename,'NSCC_peaks','NSCC_CDs_usec','-append');

set(handles.savedText,'Visible','on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on key press with focus on CFlistbox and none of its controls.
function CFlistbox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to CFlistbox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in doneButton.
function doneButton_Callback(hObject, eventdata, handles)
% hObject    handle to doneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume;
delete(handles.figure1);


function handles = updateData(hObject, handles)
% update the data from the saved file
global ExpList ROOT_dir % setup_Vowel_STMP.m should have been called already

ExpDate = handles.data.ExpDate;
UnitName = handles.data.UnitName;
FeatureIndices = handles.data.FeatureIndices;
AttenIndices = handles.data.AttenIndices;

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
handles.data.SCCs_filename = fullfile(SACSCCanal_dir,SCCs_filename);

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

FIGinfo.title=sprintf('SCC Plots [TriFiltWidth=%d] \nExp: ''%s''; Unit: %d.%02d  BF=%.2f kHz, Thr=%.f dB SPL, SR=%.1f sps, Q10=%.1f\n', ...
    CHdata.TriFiltWidth, ...
    ExpName,TrackNum,UnitNum,unit.Info.BF_kHz,unit.Info.Threshold_dBSPL,unit.Info.SR_sps,unit.Info.Q10);
FIGinfo.title = regexprep(FIGinfo.title, '[_]', '\\_'); % fix underscores

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
    handles.data.SCCpos=SCCpos;

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

    cla(handles.axes1,'reset');
    neurogram(CHdata,FIGinfo,PARAMInfo);
    legend off
end % for FeatNUM=1:NumFEATURES


%%% Clean up axes
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

% set(gca,'LineWidth',2)
% set(gca,'Box','on')

%%%%%

% Populate CFlistbox with CFs
set(handles.CFlistbox,'String',fliplr(cellfun(@num2str,num2cell(sort(CHdata.CHvals)),'UniformOutput',false)));

% Update handles structure
handles.data.CHdata=CHdata;
handles.data.FIGinfo=FIGinfo;
handles.data.PARAMInfo=PARAMInfo;
handles.data.FeatureIndices = FeatureIndices;
handles.data.AttenIndices = AttenIndices;
guidata(hObject, handles);

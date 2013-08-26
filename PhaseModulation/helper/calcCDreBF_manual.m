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

handles.SACSCCfunctions = varargin{1};
handles.SACSCCmetrics = varargin{2};
handles.BF_kHz = varargin{3};
handles.numFibers = numel(handles.SACSCCfunctions);

for ii=1:handles.numFibers
    handles.rho(ii)=...
        handles.SACSCCmetrics{ii}.SCCpeak_AB/...
        sqrt(handles.SACSCCmetrics{ii}.SACpeak_A*...
        handles.SACSCCmetrics{ii}.SACpeak_B);
    handles.CD(ii)=...
        handles.SACSCCfunctions{ii}.delays_usec(...
        handles.SACSCCfunctions{ii}.SCC_AB_avg==...
        max(handles.SACSCCfunctions{ii}.SCC_AB_avg));
end

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
varargout{1} = handles.rho;
varargout{2} = handles.CD;

delete(handles.figure1);


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
BFoctCRIT = 1/128;
CFindexTEMP = find(abs(log2(handles.data.CHdata.CHvals/handles.data.selectedCFkHz))<BFoctCRIT);
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
    neurogram2(handles.data.CHdata,handles.data.FIGinfo,handles.data.PARAMInfo);
    legend off
    
    % update detailsText
    set(handles.detailsText,'String',sprintf('CF = %1.3f kHz\nCD = %1.0f us',...
        handles.data.selectedCFkHz,CD));
    
    % SAVE CD & PEAK LOCALLY
    handles.data.selectedCD = CD;
    handles.data.selectedPeak = PEAK;
    handles.CD(CFindexTEMP) = CD;
    handles.rho(CFindexTEMP) = PEAK;
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
ChanIND = handles.data.selectedCFind;
FeatIND = handles.data.FeatureIndices;
AttIND = handles.data.AttenIndices;

if handles.data.CHdata.CHvals(ChanIND)<handles.data.CHdata.CHvals(handles.data.CHdata.BFref_ind)
    % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
    % but we want other way around, so need to flip 1st set here
    NSCC_peaks{FeatIND,AttIND}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedPeak;
    NSCC_CDs_usec{FeatIND,AttIND}{handles.data.SCCpos(ChanIND-1,1)}=...
        -handles.data.selectedCD;
elseif handles.data.CHdata.CHvals(ChanIND)>handles.data.CHdata.CHvals(handles.data.CHdata.BFref_ind)
    NSCC_peaks{FeatIND,AttIND}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedPeak;
    NSCC_CDs_usec{FeatIND,AttIND}{handles.data.SCCpos(ChanIND-1,1)}=...
        handles.data.selectedCD;
end

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


function handles = updateData(hObject, handles)
% update the data from the saved file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get general parameters for this data, e.g., all BFs, levels, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIGinfo.CONDlabel_text = '';
FIGinfo.Ylog=1;  FIGinfo.Xlog=0;  % LOG y-axis
FIGinfo.Ylabel_text='BF (kHz)';
FIGinfo.Xlabel_text='Delay (usec)';
FIGinfo.param_units = 'usec';
FIGinfo.LineStyle='-';
FIGinfo.Marker='none';

% FIGinfo.XLIMITS=[-1 1]*SCCs.paramsOUT{1,1,1}{1,1}.MAXdelay_sec*1e6;   %
FIGinfo.XLIMITS=[-1 1]*2000;   %
FIGinfo.chGAIN=.8;
PARAMInfo.param_values = [NaN NaN];

BFoctCRIT=1/128;  % Chooses as BF channel is within 1/128 octave
CHdata.TriFiltWidth=1;

FIGinfo.title='';

% CHvals = BFs
CHdata.CHvals = handles.BF_kHz;
% Find center BF (to make BOLD)... should be ch.1
CHdata.BFref_ind=1;
CHdata.BOLDchs=CHdata.BFref_ind;

for ChanIND=1:handles.numFibers
    for ParamIND=1:2
        AttIND=1;
        if mod(ParamIND,2) % odd params = full SCC sequence
            CHdata.Xdata{ChanIND,ParamIND}= ...
                handles.SACSCCfunctions{ChanIND}.delays_usec;
            if CHdata.CHvals(ChanIND)<CHdata.CHvals(CHdata.BFref_ind)
                % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
                % but we want other way around, so need to flip 1st set here
                CHdata.Ydata{ChanIND,ParamIND}= ...
                    handles.SACSCCfunctions{ChanIND}.SCC_AB_avg;
%                     fliplr(handles.SACSCCfunctions{ChanIND}.SCC_AB_avg);
            elseif ChanIND==CHdata.BFref_ind % should be ch.1
                CHdata.Ydata{ChanIND,ParamIND}= ...
                    handles.SACSCCfunctions{ChanIND}.SAC_A_avg;
            elseif CHdata.CHvals(ChanIND)>CHdata.CHvals(CHdata.BFref_ind)
                CHdata.Ydata{ChanIND,ParamIND}= ...
                    handles.SACSCCfunctions{ChanIND}.SCC_AB_avg;
            end
        else  % even params = hand-picked CDs
            if CHdata.CHvals(ChanIND)<CHdata.CHvals(CHdata.BFref_ind)
                % SCCs: CF2 re CF1, so if CF1>CF2, CD>0
                % but we want other way around, so need to flip 1st set here
                CHdata.Ydata{ChanIND,ParamIND}= ...
                    CHdata.Ydata{ChanIND,ParamIND-1}(...
                    min(find(handles.SACSCCfunctions{ChanIND}.delays_usec>=handles.CD(ChanIND))));
                CHdata.Xdata{ChanIND,ParamIND}= ...
                    handles.CD(ChanIND);
            elseif ChanIND==CHdata.BFref_ind % should be ch.1
                CHdata.Ydata{ChanIND,ParamIND}= NaN;
                CHdata.Xdata{ChanIND,ParamIND}= NaN;
            elseif CHdata.CHvals(ChanIND)>CHdata.CHvals(CHdata.BFref_ind)
                CHdata.Ydata{ChanIND,ParamIND}= ...
                    CHdata.Ydata{ChanIND,ParamIND-1}(...
                    min(find(handles.SACSCCfunctions{ChanIND}.delays_usec>=handles.CD(ChanIND))));
                CHdata.Xdata{ChanIND,ParamIND}= ...
                    handles.CD(ChanIND);
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
neurogram2(CHdata,FIGinfo,PARAMInfo);
legend off

NumFEATURES=1;
%%% Clean up axes
% set(get(h2,'Title'),'String','')
Xcorner=0.1;
Xwidth=.85;
Ycorner=0.15;
Yshift=0.2;
Ywidth=.7;
Ywidth=.94*(1-NumFEATURES*(Yshift+.01))/NumFEATURES;   %.26 for 3; .42 for 2

YLIMtemp=ylim;

TICKlength=0.02;

%%%%%

% Populate CFlistbox with CFs
set(handles.CFlistbox,'String',...
    fliplr(cellfun(@num2str,num2cell(sort(CHdata.CHvals)),'UniformOutput',false)));

% Update handles structure
handles.data.CHdata=CHdata;
handles.data.FIGinfo=FIGinfo;
handles.data.PARAMInfo=PARAMInfo;
guidata(hObject, handles);

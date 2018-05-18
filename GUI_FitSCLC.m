function varargout = GUI_FitSCLC(varargin)
% GUI_FITSCLC MATLAB code for GUI_FitSCLC.fig
%      GUI_FITSCLC, by itself, creates a new GUI_FITSCLC or raises the existing
%      singleton*.
%
%      H = GUI_FITSCLC returns the handle to a new GUI_FITSCLC or the handle to
%      the existing singleton*.
%
%      GUI_FITSCLC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_FITSCLC.M with the given input arguments.
%
%      GUI_FITSCLC('Property','Value',...) creates a new GUI_FITSCLC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_FitSCLC_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_FitSCLC_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FitSCLC

% Last Modified by GUIDE v2.5 15-May-2018 13:10:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @GUI_FitSCLC_OpeningFcn, ...
    'gui_OutputFcn',  @GUI_FitSCLC_OutputFcn, ...
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


% --- Executes just before GUI_FitSCLC is made visible.
function GUI_FitSCLC_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FitSCLC (see VARARGIN)

% Choose default command line output for GUI_FitSCLC
handles.output = hObject;

%default parameter values
handles.In.project = {'default'}; %name of the project = output file ('default')
handles.In.mA_cm2 = 1; %[1] 0=user data in mA/cm2;1=A/m2 (SI)
%ranging
handles.In.ranging = 3; %[1] toggle V/I/auto-range for least-squares fitting (1=voltage;2=current;3=auto)
handles.In.I_range = [1e1 1e3]; %[A/m2] Current density range [min max]
handles.In.V_range = [1 10]; %[V] Voltage range [min max]
% handles.Fit.info = 1; %[1] level of graphical info during auto range, 0=none;1=problems;2=all (1)
handles.In.SmoothPar = [13;... [1] nr of points in Savitzky-Golay moving average filter; must be odd (13)
    7;... %order of polynomial model in Savitzky-Golay, <nPts (7)
    0.1]; %span: fraction of data in rloess smooth (1 is all, 0 is none) (0.1)
handles.In.AutoRangePar = [10;... [1] minimum nr of points for meaningful SCLC fit (10)
    0.2;...  [V] minimum voltage to accept (0.2)
    1.7;2;...[1] minimum and preferred lowest slope of jV ([1.7 2])
    0];  %[1/V] minimum acceptable curvature of power law slope (0)
handles.In.OutputPar = [0.1 10 10]; %min V, max V, #points/decade of export bias
%fitting
handles.In.FitModel = 1; %1=Murgatroyd/Gil, 2=drift-diffusion
handles.In.muModel = 1; %1=GDM, 2=Arrhenius or 1=eGDM, 2=cGDM, 3=ET-GDM
handles.In.PositiveGamma = 1; %0=free gamma; 1=gamma must be positive, check line below
    handles.In.MinGamma = 0; %minimum gamma value (0 or -Inf), depends on the line above
handles.In.Constrained = 0; %0=free mu, gamma; 1=constrained mu,gamma
handles.In.generalPar = [3.6 0 0;... [1] relative dielectric constant & margin
    0.1 0 0.1;...                    [eV] phi1
    0.1 0 0.1];                     %[eV] phi2
handles.In.GillPar = [2.8e-5 0 -3 3;... [eV(V/m)^-0.5] field enhancement prefactor B
    600 0 1 2000];                     %[K] field enhancement offset
handles.In.GDMpar = [1e-7 0 -3 3;... [m^2/Vs] mobility at zero field at T=inf
    0.075 0 0.025 0.200];           %[eV] Gaussian disorder sigma
handles.In.ArrheniusPar = [1e-7 0 -3 3;... [m^2/Vs] mobility at zero field at T=inf
    0.250 0 0.01 0.500];                   %[eV] activation energy
handles.In.eGDMpar = [0.075 0 0.025 0.120;... [eV] Gaussian disorder sigma
    1e11 0 -2 2;...                           [1/s] attempt to hop frequency
    1.8e-9 0 1.0e-9 3.0e-9];                 %[m] inter-site distance
handles.In.ET_GDMpar = [0.2e-9 0 0.1e-9 2.0e-9]; %[m] localization length
handles.In.SpeedMode = 1; %[1] 0=all exp bias points, 1=reduced bias mesh+interpolation
handles.In.DDpar = [1.0e-9;... [m] dx - mesh size (0.1...1e-9)
    5;...    [1] In.MinIt -  min nr of iterations in self consistent loop (5)
    500;...  [1] In.MaxIt - max nr. of iterations in self-consistent loop (200)
    10e-9;... [m] In.dxCont - thickness of contact region (5e-9)
    10;...   [1] In.Qexcess - factor for exceeding SCLC charge concentration (10)
    0.10;...  [1] In.mix - weight of newest solution in mixing with old (0.1...0.2)
    1e-5;... [1] In.MaxChange - convergence criterion: max relative change in potential (1e-5)
    0];     %[1] In.UseImPot - toggles use of image potential; false is more stable (false)
%define structures for various user data
handles.DATA = [];
handles.Fit = [];
handles.curves = [];
handles.status = [];
handles.status.DataLoaded = false; %no data have yet been loaded
handles.status.RangeSet = false; %no points to fit have been identified
handles.status.Fitted = false; %no fit has been performed
handles.status.Saved = false; %data & fit have not been saved
handles.status.LastFit = ''; %name of last fit

%load in GUI
DisplayPar(handles)

%set (in)active elements
handles = SetActive(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_FitSCLC wait for user response (see UIRESUME)
% uiwait(handles.the_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_FitSCLC_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ProjectName_Callback(hObject, eventdata, handles)
% hObject    handle to ProjectName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ProjectName as text
%        str2double(get(hObject,'String')) returns contents of ProjectName as a double
handles.In.project = get(hObject,'String');
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function ProjectName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProjectName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in LoadedFiles.
function LoadedFiles_Callback(hObject, eventdata, handles)
% hObject    handle to LoadedFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LoadedFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LoadedFiles



% --- Executes during object creation, after setting all properties.
function LoadedFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadedFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadNewData.
function LoadNewData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadNewData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check if non-saved fitted data exist
if handles.status.Fitted && ~handles.status.Saved 
    choice = questdlg('Continue without saving?', ...
        'Warning','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end

%remove old figures
cla(handles.jVplot)
cla(handles.slopePlot)
cla(handles.gammaPlot)
cla(handles.mu0Plot)

%load jV data
[FileName,PathName]=uigetfile('*.txt','MultiSelect','on');

if isnumeric(FileName) %no file selected
    return
end

if ischar(FileName) %make sure FileName consists of cell arrays
    DATA.NrTemp = 1;
    FileName = mat2cell(FileName,1);
else
    DATA.NrTemp = size(FileName,2);
end

%order according to T (low to high, just because it looks nicer)
DATA.T = zeros(DATA.NrTemp,1);
DATA.L = zeros(DATA.NrTemp,1);
for k=1:DATA.NrTemp
    hlp = FileName{1,k}; %turn cell array into character array
    DATA.T(k) = str2double(hlp(end-7:end-5)); %[K] get T from file name
    DATA.L(k) = str2double(hlp(end-13:end-11)); %[nm] get L from file name
end
[DATA.T,I] = sort(DATA.T);
DATA.L = 1e-9*DATA.L(I); %[m]
FileName(1,:) = FileName(1,I);

%show in listbox
set(handles.LoadedFiles,'String',FileName);

DATA.MaxPts = 1e4; %[1] maximum number of points in jV curve
DATA.jVT = NaN(DATA.MaxPts+1,2,DATA.NrTemp); %preallocate empty array for j(V,T)
DATA.slope = NaN(DATA.MaxPts,DATA.NrTemp); %power law slope
DATA.slopeS = NaN(DATA.MaxPts,DATA.NrTemp); %smoothened slope
DATA.ind = false(DATA.MaxPts,DATA.NrTemp); %indices of points to be used in fitting
% DATA.T = zeros(DATA.NrTemp,1);

%read data and store in DATA structure
for k=1:DATA.NrTemp
    dummy = dlmread([PathName FileName{1,k}]); %,'\t' is for tab delimiter. Needed?
    dummy = abs(dummy); %program doesn't like negative bias & current
    nData = size(dummy,1);
    DATA.jVT(1:nData,:,k) = dummy;
    DATA.slope(1:nData-1,k) = ... raw power law slope
        diff(log10(dummy(:,2)),1,1)./diff(log10(dummy(:,1)),1,1);
    DATA.slopeS(1:nData-1,k) = DATA.slope(1:nData-1,k); %i.e. raw slope
%     hlp = FileName{1,k}; %turn cell array into character array
%     DATA.T(k) = str2double(hlp(end-7:end-5)); %[K] get T from file name, see header
end
if handles.In.mA_cm2==1 %user data in mA/cm2
    DATA.jVT(:,2,:) = 10*DATA.jVT(:,2,:); %[A/m2] 
end

%store data in handles
handles.DATA = DATA;
handles.curves.jVfit = zeros(DATA.NrTemp,1); %handles for fitted jV curves
handles.curves.slopeExp = zeros(DATA.NrTemp,1); %handles for raw experimental slope
handles.curves.slopeExpS = zeros(DATA.NrTemp,1); %handles for smoothened exp. slope
handles.curves.slopeFit = zeros(DATA.NrTemp,1); %handles for fitted jV slopes
handles.In.OutputParN = ceil(handles.In.OutputPar(3)*...
    log10(handles.In.OutputPar(2)/handles.In.OutputPar(1))+1); %nr of output points
handles.In.FileName = FileName;
handles.In.PathName = PathName;
handles.Fit.jV = NaN(handles.In.OutputParN,handles.DATA.NrTemp+1); %array with V,j(T1),j(T2) etc.
handles.Fit.slope = NaN(handles.In.OutputParN-1,handles.DATA.NrTemp+1); %same, for slope

%plot the raw data
handles = Make_jVplot(handles);
handles = Make_slopePlot(handles);

handles.status.DataLoaded = true; %data have been loaded
handles.status.RangeSet = false; %no points to fit have been identified
handles.status.Fitted = false; %no fit has been performed
handles.status.Saved = false; %data & fit have not been saved
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in FindRange.
function FindRange_Callback(hObject, eventdata, handles)
% hObject    handle to FindRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check if non-saved fitted data exist
if handles.status.Fitted && ~handles.status.Saved 
    choice = questdlg('Continue without saving?', ...
        'Warning','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end

handles.Fit.jV(:) = NaN;
handles.Fit.slope(:) = NaN;
for k=1:handles.DATA.NrTemp %loop over T
    nData = sum(~isnan(handles.DATA.jVT(:,2,k))); %nr of points in jV @ T
    data = handles.DATA.jVT(1:nData,:,k); %all existing points in jV @ T
    switch handles.In.ranging
        case 1 %voltage range
            ind = data(:,1)>=handles.In.V_range(1) & data(:,1)<=handles.In.V_range(2);
        case 2 %V's corresponding to current range
            ind = data(:,2)>=handles.In.I_range(1) & data(:,2)<=handles.In.I_range(2);
        otherwise %auto range
            %some data processing
            slope = handles.DATA.slope(1:nData-1,k); %power law slope @ T
            ind_voltage = data(:,1)>=handles.In.AutoRangePar(2); %indices of points above threshold voltage
            ind_slope = [false; slope>=handles.In.AutoRangePar(4)]; %indices of points with desired slope
            slopeS = handles.DATA.slopeS(1:nData-1,k);
            curvature = diff(slopeS,1,1)./diff(log10(data(2:end,1)),1,1); %curvature
            ind_curv = [false; curvature>=handles.In.AutoRangePar(5); false]; %indices of points with (roughly) increasing slope
            
            %find suitable range
            ind = ind_voltage&ind_slope&ind_curv; %Murgatroyd requires an increasing slope>2 for jV
            if sum(ind,1)<handles.In.AutoRangePar(1) %problems...
                if (sum(ind_slope,1)>=handles.In.AutoRangePar(1))&&...
                        (sum(ind_curv,1)>=handles.In.AutoRangePar(1)) %...due combination of factors
                    %no solution!
                    disp(['Warning: questionable data in dataset',num2str(k)])
                end
                if sum(ind_slope,1)<handles.In.AutoRangePar(1) %...due to too low slope
                    ind_slope = [false; slope>=handles.In.AutoRangePar(3)]; %solution: try lower target slope
                    if sum(ind_slope,1)<handles.In.AutoRangePar(1) %still problems...
                        %no solution!
                        disp(['Warning: data very far from SCLC in dataset',num2str(k)])
                    else %problem sort of solved...
                        disp(['Warning: SCLC not reached in dataset ',num2str(k)])
                    end
                end
                if sum(ind_curv,1)<handles.In.AutoRangePar(1) %...due to decreasing slope
                    ind_curv(:) = true; %solution: ignore the curvature
                    disp(['Warning: wrong curvature in dataset',num2str(k)])
                end
                ind = ind_voltage&ind_slope&ind_curv; %update ind with relaxed constraints
            end
    end
    handles.DATA.jVT(end,1,k) = sum(ind,1); %number of meaningful points
    handles.DATA.ind(1:nData,k) = ind; %corresponding indices
end

%update jV and slope plots
handles = Make_jVplot(handles);
handles = Make_slopePlot(handles);

if max(handles.DATA.jVT(end,1,:))>=1
    handles.status.RangeSet = true; %at least 1 point to fit has been identified
else
    handles.status.RangeSet = false; %no points to fit have been identified
end
handles.status.Fitted = false; %no fit has been performed
handles.status.Saved = false; %data & fit have not been saved
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure




% --- Executes on selection change in RangeMethod.
function RangeMethod_Callback(hObject, eventdata, handles)
% hObject    handle to RangeMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns RangeMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from RangeMethod
handles.In.ranging = get(handles.RangeMethod,'Value');
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function RangeMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RangeMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IrangeMin_Callback(hObject, eventdata, handles)
% hObject    handle to IrangeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IrangeMin as text
%        str2double(get(hObject,'String')) returns contents of IrangeMin as a double
handles.In.I_range(1) = str2double(get(hObject,'String'));
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function IrangeMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IrangeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function IrangeMax_Callback(hObject, eventdata, handles)
% hObject    handle to IrangeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IrangeMax as text
%        str2double(get(hObject,'String')) returns contents of IrangeMax as a double
handles.In.I_range(2) = str2double(get(hObject,'String'));
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function IrangeMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IrangeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VrangeMin_Callback(hObject, eventdata, handles)
% hObject    handle to VrangeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VrangeMin as text
%        str2double(get(hObject,'String')) returns contents of VrangeMin as a double
handles.In.V_range(1) = str2double(get(hObject,'String'));
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function VrangeMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VrangeMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VrangeMax_Callback(hObject, eventdata, handles)
% hObject    handle to VrangeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VrangeMax as text
%        str2double(get(hObject,'String')) returns contents of VrangeMax as a double
handles.In.V_range(2) = str2double(get(hObject,'String'));
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function VrangeMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VrangeMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when entered data in editable cell(s) in AutoRangeTable.
function AutoRangeTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to AutoRangeTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.AutoRangePar = get(handles.AutoRangeTable,'Data');
guidata(hObject, handles); %Update handles structure



function handles = Make_jVplot(handles)

cla(handles.jVplot)
hold(handles.jVplot,'on')
for k=1:handles.DATA.NrTemp
    sz = 5*ones(handles.DATA.MaxPts,1);  %[pts] size of unused data points (5)
    sz(handles.DATA.ind(:,k)) = 25; %[pts] size of data points used in fit (25)
    if handles.DATA.NrTemp==1
        c = [1 0 0];
    else
        c = [(k-1)/(handles.DATA.NrTemp-1),0,1-(k-1)/(handles.DATA.NrTemp-1)]; %color as RGB triplet
    end
    scatter(handles.DATA.jVT(1:end-1,1,k),handles.DATA.jVT(1:end-1,2,k),...
        sz,c,'Parent',handles.jVplot); %experiment
    handles.curves.jVfit(k) = ... fit
        line(handles.Fit.jV(:,1),handles.Fit.jV(:,k+1),'LineWidth',2,'Color',c,'Parent',handles.jVplot);
end
set(handles.jVplot,'Xscale','log');
set(handles.jVplot,'Yscale','log');
grid on
xlabel(handles.jVplot,'Voltage [V]')
ylabel(handles.jVplot,'Current Density [A/m^{2}]')
title(handles.jVplot,'jV curves + fits')



function handles = Make_slopePlot(handles)

cla(handles.slopePlot)
hold(handles.slopePlot,'on')
for k=1:handles.DATA.NrTemp
    sz = 5*ones(handles.DATA.MaxPts,1);  %[pts] size of unused data points (5)
    sz(handles.DATA.ind(:,k)) = 25; %[pts] size of data points used in fit (25)
    if handles.DATA.NrTemp==1
        c = [1 0 0];
    else
        c = [(k-1)/(handles.DATA.NrTemp-1),0,1-(k-1)/(handles.DATA.NrTemp-1)]; %color as RGB triplet
    end
    handles.curves.slopeExp(k) = ... raw experimental slope
        scatter(handles.DATA.jVT(1:end-1,1,k),handles.DATA.slope(:,k),...
        sz,c,'Parent',handles.slopePlot);
    handles.curves.slopeExpS(k) = ... smoothened experimental slope
        line(handles.DATA.jVT(1:end-1,1,k),handles.DATA.slopeS(:,k),...
        'LineWidth',0.5,'Color','black','LineStyle','--','Parent',handles.slopePlot);
    handles.curves.slopeFit(k) = ... fit
        line(handles.Fit.slope(:,1),handles.Fit.slope(:,k+1),'LineWidth',2,'Color',c,'Parent',handles.slopePlot);
end
set(handles.slopePlot,'Xscale','log');
set(handles.slopePlot,'Yscale','lin');
set(handles.slopePlot,'YLim',[0,5]); %fixed y-range, comment-out for auto-range
grid on
xlabel(handles.slopePlot,'Voltage [V]')
ylabel(handles.slopePlot,'d(log(j))/d(log(V)) [-]')
title(handles.slopePlot,'power law slope + fits')



function Make_gammaPlot(handles)

cla(handles.gammaPlot)
hold(handles.gammaPlot,'on')
x = handles.DATA.T.^-1;
scatter(x,handles.Fit.gamma,'Parent',handles.gammaPlot)
gammaFit = Gill(handles.In.GillPar(1,2),handles.In.GillPar(2,2),...
    handles.DATA.T,handles.In.MinGamma);
line(x,gammaFit,'Parent',handles.gammaPlot)
set(handles.gammaPlot,'Xscale','lin');
set(handles.gammaPlot,'Yscale','lin');
grid(handles.gammaPlot,'on')
xlabel(handles.gammaPlot,'1/T [K^{-1}]')
ylabel(handles.gammaPlot,'gamma [(V/m)^{-1/2}]')



function Make_mu0Plot(handles)

cla(handles.mu0Plot)
hold(handles.mu0Plot,'on')
x = handles.DATA.T.^-handles.In.C(2);
scatter(x,handles.Fit.mu0,'Parent',handles.mu0Plot)
if handles.In.muModel==1 %mu0 according to GDM
    muFit = mu0_T(handles.In.GDMpar(1,2),handles.In.C,...
        handles.In.GDMpar(2,2),handles.DATA.T);
else %Arrhenius
    muFit = mu0_T(handles.In.ArrheniusPar(1,2),handles.In.C,...
        handles.In.ArrheniusPar(2,2),handles.DATA.T);
end
line(x,muFit,'Parent',handles.mu0Plot)
set(handles.mu0Plot,'Xscale','lin');
set(handles.mu0Plot,'Yscale','log');
grid(handles.mu0Plot,'on')
xlabel(handles.mu0Plot,['1/T^' num2str(handles.In.C(2)) ' [K^{' num2str(-handles.In.C(2)) '}]'])
ylabel(handles.mu0Plot,'mu0 [m^{2}/ Vs]')



% --- Executes when entered data in editable cell(s) in SmoothTable.
function SmoothTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to SmoothTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.SmoothPar = get(handles.SmoothTable,'Data');
%some protection against wrong input: nPts must be odd and larger than order
handles.In.SmoothPar(1) = round(max(handles.In.SmoothPar(1),handles.In.SmoothPar(2)+1));
if mod(handles.In.SmoothPar(1),2)==0 
    handles.In.SmoothPar(1) = handles.In.SmoothPar(1)+1;
end
set(handles.SmoothTable,'Data',handles.In.SmoothPar);
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in Smooth.
function Smooth_Callback(hObject, eventdata, handles)
% hObject    handle to Smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%check if non-saved fitted data exist
if handles.status.Fitted && ~handles.status.Saved 
    choice = questdlg('Continue without saving?', ...
        'Warning','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end
%data processing
handles.DATA.ind(:) = false;
handles.Fit.jV(:) = NaN;
handles.Fit.slope(:) = NaN;
for k=1:handles.DATA.NrTemp
    ind =~isnan(handles.DATA.jVT(:,2,k));
    data = handles.DATA.jVT(ind,:,k);
    dataS = smooth(data(:,1),data(:,2),handles.In.SmoothPar(1),...
        'sgolay',handles.In.SmoothPar(2)); %noise reduction
    slope = diff(log10(dataS),1,1)./diff(log10(data(:,1)),1,1); %power law slope
    handles.DATA.slope(1:sum(ind)-1,k) = slope;        
    handles.DATA.slopeS(1:sum(ind)-1,k) = ...
        smooth(data(2:end,1),slope,handles.In.SmoothPar(3),'rloess'); %smooth slope
end
handles = Make_jVplot(handles); %update plots
handles = Make_slopePlot(handles);

handles.status.RangeSet = false; %no points to fit have been identified
handles.status.Fitted = false; %no fit has been performed
handles.status.Saved = false; %data & fit have not been saved
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes on selection change in FitModel.
function FitModel_Callback(hObject, eventdata, handles)
% hObject    handle to FitModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FitModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FitModel
handles.In.FitModel = get(hObject,'Value');
if handles.In.FitModel==1 %Murgatroyd/Gill
    handles.muModel.String = {'GDM';'Arrhenius'};
else %drift-diffusion
    handles.muModel.String = {'eGDM';'cGDM';'ET-GDM'};
end
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function FitModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FitModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in mumodel.
function muModel_Callback(hObject, eventdata, handles)
% hObject    handle to mumodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mumodel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mumodel
handles.In.muModel = get(hObject,'Value');
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes during object creation, after setting all properties.
function muModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mumodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ConstrainedFit.
function ConstrainedFit_Callback(hObject, eventdata, handles)
% hObject    handle to ConstrainedFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ConstrainedFit
handles.In.Constrained = get(hObject,'Value');
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in ShowGuess.
function ShowGuess_Callback(hObject, eventdata, handles)
% hObject    handle to ShowGuess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Basically this is a short version of Fit_Callback() - see there
guess.Vbi = handles.In.generalPar(2,1)-handles.In.generalPar(3,1); %[V] built-in voltage
if handles.In.muModel==1
    handles.In.C = [2/3;2]; %makes mu0(T) follow GDM, see function mu0_T
else %Arrhenius
    handles.In.C = [1;1];
end
handles.In.OutputParN = ceil(handles.In.OutputPar(3)*...
    log10(handles.In.OutputPar(2)/handles.In.OutputPar(1))+1); %nr of output points
guess.jV = NaN(handles.In.OutputParN,handles.DATA.NrTemp+1); %array with V,j(T1),j(T2) etc.
guess.jV(:,1) = logspace(log10(handles.In.OutputPar(1)),log10(handles.In.OutputPar(2)),...
    handles.In.OutputParN)'; %[V] bias range for exported jV curves
guess.slope = NaN(handles.In.OutputParN-1,handles.DATA.NrTemp+1); %same, for slope
guess.slope(:,1) = guess.jV(2:end,1);

if handles.In.FitModel==2
    In.phi0 = [handles.In.generalPar(2,1),handles.In.generalPar(3,1)];
    In.sigma = handles.In.eGDMpar(1,1);
    In.nu0 = handles.In.eGDMpar(2,1);
    In.aNN = handles.In.eGDMpar(3,1);
    if handles.In.muModel==3 %ET-GDM
        In.alpha = handles.In.ET_GDMpar(1,1); %[m]
        In.L = min(handles.DATA.L); %worst-case, L is overruled below
        In = MakeLookup(In,handles.DATA.T,handles.In.OutputPar(2)); %requires lookup table
    end
    In.muModel = handles.In.muModel;
    In.DDpar = handles.In.DDpar;
    In.bias = guess.jV(:,1);
end

for k=1:handles.DATA.NrTemp
    if handles.In.FitModel==1 %Murgatroyd/Gill
        %get mu0 and gamma vs. T from fitted parameterss
        guess.gamma(k) = Gill(handles.In.GillPar(1,1),handles.In.GillPar(2,1),...
            handles.DATA.T(k),handles.In.MinGamma);
        if handles.In.muModel==1 %mu0 according to GDM
            guess.mu0(k) = mu0_T(handles.In.GDMpar(1,1),handles.In.C,...
                handles.In.GDMpar(2,1),handles.DATA.T(k));
        else %Arrhenius
            guess.mu0(k) = mu0_T(handles.In.ArrheniusPar(1,1),handles.In.C,...
                handles.In.ArrheniusPar(2,1),handles.DATA.T(k));
        end
        guess.jV(:,k+1) = SCLC(guess.mu0(k),guess.gamma(k),...
            handles.DATA.L(k),handles.In.generalPar(1,1),...
            guess.jV(:,1)-guess.Vbi); %[A/m2] j(V)
    else %drift-diffusion/eGDM
        In.T = handles.DATA.T(k); %load required parameters in a temporary In. structure
        In.epsR = handles.In.generalPar(1,1);
        In.L = handles.DATA.L(k);
        [~,Out] = ODDD(In); %[In,Out], run ODDD, output In is not needed
        guess.jV(:,k+1) = Out.j; %extract j from Out structure
    end
    guess.slope(:,k+1) = ... power law slope
        diff(log10(guess.jV(:,k+1)),1,1)./diff(log10(guess.jV(:,1)),1,1);
    set(handles.curves.jVfit(k),'XData',guess.jV(:,1),'YData',guess.jV(:,k+1)) %update jV plot
    set(handles.curves.slopeFit(k),'XData',guess.slope(:,1),'YData',guess.slope(:,k+1)) %and slope
end



% --- Executes on button press in FitIt.
function FitIt_Callback(hObject, eventdata, handles)
% hObject    handle to FitIt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Init
handles.In.epsR = handles.In.generalPar(1,1); %[1] relative dielectric constant
handles.In.C = [0;0]; %just so it exists
handles.Fit.mu0 = zeros(handles.DATA.NrTemp,1); %[m2/Vs] zero field&density mobility at each T
if handles.In.FitModel==1 %Murgatroyd/Gill
    handles.In.Vbi = handles.In.generalPar(2,1)-handles.In.generalPar(3,1); %[V] built-in voltage
    handles.In.dVbi = handles.In.generalPar(2,3)+handles.In.generalPar(3,3); %[V] margin in Vbi
    if handles.In.muModel==1
        handles.In.C = [2/3;2]; %makes mu0(T) follow GDM, see function mu0_T
    else
        handles.In.C = [1;1]; %makes mu0(T) follow Arrhenius behavior
    end
    handles.Fit.gamma = zeros(handles.DATA.NrTemp,1); %[(V/m)^-0.5] field enhancement factor
end
handles.In.OutputParN = ceil(handles.In.OutputPar(3)*...
    log10(handles.In.OutputPar(2)/handles.In.OutputPar(1))+1); %nr of output points
handles.Fit.jV = NaN(handles.In.OutputParN,handles.DATA.NrTemp+1); %array with V,j(T1),j(T2) etc.
handles.Fit.jV(:,1) = logspace(log10(handles.In.OutputPar(1)),log10(handles.In.OutputPar(2)),...
    handles.In.OutputParN)'; %[V] bias range for exported jV curves
handles.Fit.slope = NaN(handles.In.OutputParN-1,handles.DATA.NrTemp+1); %same, for slope
handles.Fit.slope(:,1) = handles.Fit.jV(2:end,1);

%Set fitit parameter guess & limits, first assuming GDM
if handles.In.FitModel==1 %Murgatroyd/Gill
    Fit.par = [...
        [handles.In.Vbi-handles.In.dVbi handles.In.Vbi handles.In.Vbi+handles.In.dVbi];... 1) built-in voltage
        [log10(handles.In.GDMpar(1,1))+handles.In.GDMpar(1,3)... 2) mobility prefactor (only constrained fit)
        log10(handles.In.GDMpar(1,1)) ...
        log10(handles.In.GDMpar(1,1))+handles.In.GDMpar(1,4)];...
        [handles.In.GDMpar(2,3) handles.In.GDMpar(2,1) handles.In.GDMpar(2,4)];... 3) disorder (only constrained fit)
        [log10(handles.In.GillPar(1,1))+handles.In.GillPar(1,3)... 4) field enhancement prefactor (only constrained fit)
        log10(handles.In.GillPar(1,1))...
        log10(handles.In.GillPar(1,1))+handles.In.GillPar(1,4)];...
        [handles.In.GillPar(2,3) handles.In.GillPar(2,1) handles.In.GillPar(2,4)]]; %5) field enhancement offset (only constrained)
    handles.status.LastFit = 'GDM'; %name of last fit
    if handles.In.muModel==2 %update if Arrhenius
        Fit.par(2,:) = [...
            log10(handles.In.ArrheniusPar(1,1))+handles.In.ArrheniusPar(1,3)... 2) mobility prefactor
            log10(handles.In.ArrheniusPar(1,1)) ...
            log10(handles.In.ArrheniusPar(1,1))+handles.In.ArrheniusPar(1,4)];
        Fit.par(3,:) = [... 3) activation energy
            handles.In.ArrheniusPar(2,3) handles.In.ArrheniusPar(2,1) handles.In.ArrheniusPar(2,4)];
        handles.status.LastFit = 'Arrhenius'; %name of last fit
    end
else %drift-diffusion/eGDM
    Fit.par = [...
        handles.In.generalPar(2,1)-handles.In.generalPar(2,3)... 1) left electron injection barrier In.phi0(1)
        handles.In.generalPar(2,1)...
        handles.In.generalPar(2,1)+handles.In.generalPar(2,3);...
        handles.In.generalPar(3,1)-handles.In.generalPar(3,3)... 2) right electron injection barrier In.phi0(2)
        handles.In.generalPar(3,1)...
        handles.In.generalPar(3,1)+handles.In.generalPar(3,3);...
        handles.In.eGDMpar(1,3) handles.In.eGDMpar(1,1) handles.In.eGDMpar(1,4);... 3) disorder In.sigma
        log10(handles.In.eGDMpar(2,1))+handles.In.eGDMpar(2,3)... 4) attempt frequency In.nu0
        log10(handles.In.eGDMpar(2,1))...
        log10(handles.In.eGDMpar(2,1))+handles.In.eGDMpar(2,4);...
        1e9*[handles.In.eGDMpar(3,3) handles.In.eGDMpar(3,1) handles.In.eGDMpar(3,4)]]; %5) nearest neighbor distance in nm In.aNN
    if handles.In.muModel==3 %extend if effective temperature version of GDM
        Fit.par(6,:) = 1e9*...
            [handles.In.ET_GDMpar(1,3) handles.In.ET_GDMpar(1,1) handles.In.ET_GDMpar(1,4)]; %6) localization length in nm In.alpha
    end
    handles.status.LastFit = 'eGDM'; %name of last fit
end

%Set fmincon inputs
if handles.In.FitModel==2 || handles.In.Constrained %DD & constrained: all parameters are globally fitted
    lb = Fit.par(:,1)'; %lower bound
    p0 = Fit.par(:,2)'; %initial guess
    ub = Fit.par(:,3)';	%upper bound
else %free Murgatroyd fit: only keep Vbi in global fitting routine
    lb = Fit.par(1,1)'; %lower bound
    p0 = Fit.par(1,2)'; %initial guess
    ub = Fit.par(1,3)'; %upper bound
end
A = [];
b = [];
Aeq = [];
beq = [];

%% FitIt!
options = optimoptions('fmincon','Display','none','PlotFcns',@optimplotfval);
[Fit.p,Fit.Error] = ...
    fmincon(@(p0) FitSCLCerror(p0,handles.DATA,handles.In,Fit.par,handles.In.C),p0,...
    A,b,Aeq,beq,lb,ub,[],options);
handles.status.Fitted = true; %a fit has been performed
handles.status.Saved = false; %data & fit have not been saved

%% post-processing
%store FitIt parameters
if handles.In.FitModel==1 %Murgatroyd/Gill
    handles.In.generalPar(2,2) = handles.In.generalPar(2,1); %phi1, undetermined
    handles.In.generalPar(3,2) = handles.In.generalPar(2,1)-Fit.p(1); %Vbi
    if handles.In.Constrained
        handles.In.GillPar(1,2) = 10^Fit.p(4); %B
        handles.In.GillPar(2,2) = Fit.p(5); %T0
        if handles.In.muModel==1 %GDM
            handles.In.GDMpar(1,2) = 10^Fit.p(2); %muStar
            handles.In.GDMpar(2,2) = Fit.p(3); %sigma
        else %Arrhenius
            handles.In.ArrheniusPar(1,2) = 10^Fit.p(2); %muStar
            handles.In.ArrheniusPar(2,2) = Fit.p(3); %E_act
        end
    end
else %drift-diffusion + e/c/ET-GDM
    %store fitit parameters in GUI table
    handles.In.generalPar(2,2) = Fit.p(1); %phi1
    handles.In.generalPar(3,2) = Fit.p(2); %phi2
    handles.In.eGDMpar(1,2) = Fit.p(3); %sigma
    handles.In.eGDMpar(2,2) = 10^Fit.p(4); %nu0
    handles.In.eGDMpar(3,2) = 1e-9*Fit.p(5); %aNN
    if handles.In.muModel==3 %ET-GDM
        handles.In.ET_GDMpar(1,2) = 1e-9*Fit.p(6); %alpha
    end
    %store fitit parameters in In. to calculate corresponding jV
    In.epsR  = handles.In.generalPar(1,1);
    In.phi0  = [Fit.p(1),Fit.p(2)];
    In.sigma = Fit.p(3);
    In.nu0   = 10^Fit.p(4);
    In.aNN   = 1e-9*Fit.p(5);
    if handles.In.muModel==3 %ET-GDM
        In.alpha = 1e-9*Fit.p(6);
        In.L = min(handles.DATA.L); %worst-case, L is overruled below
        In = MakeLookup(In,handles.DATA.T,handles.In.OutputPar(2)); %requires lookup table
    end
    In.muModel = handles.In.muModel;
    In.DDpar = handles.In.DDpar;
end

%calculate jV from FitIt parameters
for k=1:handles.DATA.NrTemp
    if handles.In.FitModel==1 %Murgatroyd/Gill
        if handles.In.Constrained
            %get mu0 and gamma vs. T directly from fitted parameters
            handles.Fit.gamma(k) = Gill(handles.In.GillPar(1,2),handles.In.GillPar(2,2),...
                handles.DATA.T(k),handles.In.MinGamma);
            if handles.In.muModel==1 %mu0 according to GDM
                handles.Fit.mu0(k) = mu0_T(handles.In.GDMpar(1,2),handles.In.C,...
                    handles.In.GDMpar(2,2),handles.DATA.T(k));
            else %Arrhenius
                handles.Fit.mu0(k) = mu0_T(handles.In.ArrheniusPar(1,2),handles.In.C,...
                    handles.In.ArrheniusPar(2,2),handles.DATA.T(k));
            end
        else
            %mu0 and gamma need to be calculated by re-fitting at final Vbi
            j0 = handles.DATA.jVT(handles.DATA.ind(:,k),2,k); %[A/m2] current to fit
            bias = handles.DATA.jVT(handles.DATA.ind(:,k),1,k); %[V] corresponding bias
            [~,handles.Fit.mu0(k),handles.Fit.gamma(k)] = ... [j,mu,gamma], j not needed here
                FreeFit(handles.DATA.L(k),handles.In.generalPar(1,1),...
                bias-Fit.p(1),handles.DATA.T(k),handles.In.MinGamma,j0,Fit.par,handles.In.C); %
        end
        handles.Fit.jV(:,k+1) = SCLC(handles.Fit.mu0(k),handles.Fit.gamma(k),...
            handles.DATA.L(k),handles.In.generalPar(1,1),...
            handles.Fit.jV(:,1)-Fit.p(1)); %[A/m2] j(V)
    else %drift-diffusion & e/c/ET-GDM
        In.T = handles.DATA.T(k); %more storage in In. structure
        In.L = handles.DATA.L(k);
        In.bias = handles.Fit.jV(:,1);
        [In,Out] = ODDD(In);
        handles.Fit.jV(:,k+1) = Out.j; %[A/m2] j(V)
        handles.Fit.mu0(k) = In.mu0 ;%zero-field & density mu
    end
    handles.Fit.slope(:,k+1) = ... power law slope
        diff(log10(handles.Fit.jV(:,k+1)),1,1)./diff(log10(handles.Fit.jV(:,1)),1,1);
    set(handles.curves.jVfit(k),'XData',handles.Fit.jV(:,1),'YData',handles.Fit.jV(:,k+1)) %update jV plot
    set(handles.curves.slopeFit(k),'XData',handles.Fit.slope(:,1),'YData',handles.Fit.slope(:,k+1)) %and slope
end

%further analysis in case of Murgatroyd
if handles.In.FitModel==1 && ~handles.In.Constrained && handles.DATA.NrTemp>1
    %analyze gamma(T) in terms of Gill expression eq.2
    fun = @(p) 1e8*sum( (handles.Fit.gamma - ...
        Gill(p(1),p(2),handles.DATA.T,handles.In.MinGamma) ).^2 );
    p = fminsearch(fun,[handles.In.GillPar(1,1) handles.In.GillPar(2,1)]);
    handles.In.GillPar(1,2) = p(1);
    handles.In.GillPar(2,2) = p(2);
    
    %analyze mu0(T)
    x = handles.DATA.T.^-handles.In.C(2);
    y = log(handles.Fit.mu0);
    p = polyfit(x,y,1);
    if handles.In.muModel==1 %GDM: mu0(T)=muStar*exp[-(C*sigma/kT)^2)
        handles.In.GDMpar(1,2) = exp(p(2));
        handles.In.GDMpar(2,2) =(1.380650324e-23/1.602176487e-19)*(-p(1))^(1/handles.In.C(2))/handles.In.C(1);
    else
        handles.In.ArrheniusPar(1,2) = exp(p(2));
        handles.In.ArrheniusPar(2,2) =(1.380650324e-23/1.602176487e-19)*(-p(1))^(1/handles.In.C(2))/handles.In.C(1);
    end
end

if handles.In.FitModel==1 && handles.DATA.NrTemp>1
    Make_gammaPlot(handles)
    Make_mu0Plot(handles)
end

DisplayPar(handles) %load fitted parameters into GUI
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure
h = msgbox('Fit Completed');

% --- Executes when entered data in editable cell(s) in generalFitTable.
function generalFitTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to generalFitTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.generalPar = get(hObject,'Data');
if handles.In.generalPar(1,3)~=0 % protection to show that epsR cannot be fitted
    handles.In.generalPar(1,3) = 0;
    set(hObject,'Data',handles.In.generalPar);
end
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in GillTable.
function GillTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to GillTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.GillPar = get(hObject,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in GDMtable.
function GDMtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to GDMtable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.GDMpar = get(handles.GDMtable,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in ArrheniusTable.
function ArrheniusTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ArrheniusTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.ArrheniusPar = get(hObject,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in eGDMtable.
function eGDMtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to eGDMtable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.eGDMpar = get(hObject,'Data');
guidata(hObject, handles); %Update handles structure

function DisplayPar(handles)
%displays all current parameters in the GUI
set(handles.ProjectName,'String',handles.In.project);
set(handles.mA_cm2,'Value',handles.In.mA_cm2);
set(handles.SmoothTable,'Data',handles.In.SmoothPar);
set(handles.RangeMethod,'Value',handles.In.ranging);
set(handles.IrangeMin,'String',num2str(handles.In.I_range(1)));
set(handles.IrangeMax,'String',num2str(handles.In.I_range(2)));
set(handles.VrangeMin,'String',num2str(handles.In.V_range(1)));
set(handles.VrangeMax,'String',num2str(handles.In.V_range(2)));
set(handles.AutoRangeTable,'Data',handles.In.AutoRangePar);
set(handles.OutputTable,'Data',handles.In.OutputPar);

set(handles.FitModel,'Value',handles.In.FitModel);
set(handles.muModel,'Value',handles.In.muModel);
set(handles.PositiveGamma,'Value',handles.In.PositiveGamma);
set(handles.ConstrainedFit,'Value',handles.In.Constrained);
set(handles.generalFitTable,'Data',handles.In.generalPar);
set(handles.GillTable,'Data',handles.In.GillPar);
set(handles.GDMtable,'Data',handles.In.GDMpar);
set(handles.ArrheniusTable,'Data',handles.In.ArrheniusPar);
set(handles.eGDMtable,'Data',handles.In.eGDMpar);
set(handles.ET_GDMtable,'Data',handles.In.ET_GDMpar);
set(handles.SpeedMode,'Value',handles.In.SpeedMode);
set(handles.DDtable,'Data',handles.In.DDpar);


% --- Executes when entered data in editable cell(s) in OutputTable.
function OutputTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to OutputTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.OutputPar = get(hObject,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in DDtable.
function DDtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to DDtable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.DDpar = get(hObject,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%save jV and slope for all T
dlmwrite([handles.In.PathName,char(handles.In.project),' jV fit','.txt'],...
    [NaN,handles.DATA.T']); %header with T's
dlmwrite([handles.In.PathName,char(handles.In.project),' jV fit','.txt'],...
    handles.Fit.jV,'-append'); %txt-file with fitted jV
dlmwrite([handles.In.PathName,char(handles.In.project),' slope fit','.txt'],...
    [NaN,handles.DATA.T']); %header with T's
dlmwrite([handles.In.PathName,char(handles.In.project),' slope fit','.txt'],...
    handles.Fit.slope,'-append'); %slope of fitted jV
DATA = handles.DATA; In = handles.In; Fit = handles.Fit;
save([handles.In.PathName,char(handles.In.project),' all'],...
    'DATA','In','Fit'); %mat-file

%save parameters of last completed fit
switch handles.status.LastFit
    case 'GDM'
        T = table([get(handles.generalFitTable,'RowName');...
            get(handles.GillTable,'RowName');...
            get(handles.GDMtable,'RowName')],...
            [handles.In.generalPar(:,1);handles.In.GillPar(:,1);handles.In.GDMpar(:,1)],...
            [handles.In.generalPar(:,2);handles.In.GillPar(:,2);handles.In.GDMpar(:,2)],...
        'VariableNames',{'parameter','guess','fit'});
    case 'Arrhenius'
        T = table([get(handles.generalFitTable,'RowName');...
            get(handles.GillTable,'RowName');...
            get(handles.ArrheniusTable,'RowName')],...
            [handles.In.generalPar(:,1);handles.In.GillPar(:,1);handles.In.ArrheniusPar(:,1)],...
            [handles.In.generalPar(:,2);handles.In.GillPar(:,2);handles.In.ArrheniusPar(:,2)],...
        'VariableNames',{'parameter','guess','fit'});
    case 'eGDM'
        if handles.In.muModel<3 %eGDM or cGDM
            T = table([get(handles.generalFitTable,'RowName');...
                get(handles.eGDMtable,'RowName')],...
                [handles.In.generalPar(:,1);handles.In.eGDMpar(:,1)],...
                [handles.In.generalPar(:,2);handles.In.eGDMpar(:,2)],...
                'VariableNames',{'parameter','guess','fit'});
        else %ET-GDM
            T = table([get(handles.generalFitTable,'RowName');...
                get(handles.eGDMtable,'RowName');...
                get(handles.ET_GDMtable,'RowName')],...
                [handles.In.generalPar(:,1);handles.In.eGDMpar(:,1);handles.In.ET_GDMpar(:,1)],...
                [handles.In.generalPar(:,2);handles.In.eGDMpar(:,2);handles.In.ET_GDMpar(:,2)],...
                'VariableNames',{'parameter','guess','fit'});
        end
end
writetable(T,[handles.In.PathName,char(handles.In.project),' par','.txt'])

%save mu0 and gamma (if applicable) vs. T
if ~strcmp(handles.status.LastFit,'eGDM') %i.e. for Murgatroyd/Gill
    T = table(handles.DATA.T,handles.Fit.mu0,handles.Fit.gamma,...
        'VariableNames',{'T','mu0','gamma'});
else
    T = table(handles.DATA.T,handles.Fit.mu0,...
        'VariableNames',{'T','mu0'});
end
writetable(T,[handles.In.PathName,char(handles.In.project),' mu_gamma','.txt'])

handles.status.Saved = true; %data & fit have been saved
guidata(hObject, handles); %Update handles structure

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%check if non-saved fitted data exist
if handles.status.Fitted && ~handles.status.Saved 
    choice = questdlg('Continue without saving?', ...
        'Warning','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end
handles.DATA.ind(:) = false;
handles.DATA.slope(:) = NaN;
handles.DATA.slopeS(:) = NaN;
handles.Fit.jV(:) = NaN;
handles.Fit.slope(:) = NaN;
for k=1:handles.DATA.NrTemp
    nData = sum(~isnan(handles.DATA.jVT(:,2,k))); %nr of points in jV @ T
    handles.DATA.slope(1:nData-1,k) = ... raw power law slope
        diff(log10(handles.DATA.jVT(1:nData,2,k)),1,1)./...
        diff(log10(handles.DATA.jVT(1:nData,1,k)),1,1);
    handles.DATA.slopeS(1:nData-1,k) = handles.DATA.slope(1:nData-1,k); %i.e. raw slope
%     set(handles.curves.slopeExpS(k),'YData',handles.DATA.slopeS(:,k)) %update plot
end
handles = Make_jVplot(handles);
handles = Make_slopePlot(handles);

handles.status.RangeSet = false; %no points to fit have been identified
handles.status.Fitted = false; %no fit has been performed
handles.status.Saved = false; %data & fit have not been saved
handles = SetActive(handles); %set (in)active elements
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in SpeedMode.
function SpeedMode_Callback(hObject, eventdata, handles)
% hObject    handle to SpeedMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SpeedMode
handles.In.SpeedMode = get(hObject,'Value');
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in mA_cm2.
function mA_cm2_Callback(hObject, eventdata, handles)
% hObject    handle to mA_cm2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mA_cm2
handles.In.mA_cm2 = get(hObject,'Value');
guidata(hObject, handles); %Update handles structure


% --- Executes on button press in PositiveGamma.
function PositiveGamma_Callback(hObject, eventdata, handles)
% hObject    handle to PositiveGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PositiveGamma
handles.In.PositiveGamma = get(hObject,'Value');
if handles.In.PositiveGamma
    handles.In.MinGamma = 0; %allow only gamma>0
else
    handles.In.MinGamma =-Inf; %allow all gamma's
end
guidata(hObject, handles); %Update handles structure


% --- Executes when entered data in editable cell(s) in ET_GDMtable.
function ET_GDMtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ET_GDMtable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
handles.In.ET_GDMpar = get(handles.ET_GDMtable,'Data');
guidata(hObject, handles); %Update handles structure


% --- Executes when user attempts to close the_GUI.
function the_GUI_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to the_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if handles.status.Fitted && ~handles.status.Saved 
    choice = questdlg('Continue without saving?', ...
        'Warning','Yes','No','No');
    if strcmp(choice,'No')
        return
    end
end
delete(hObject);

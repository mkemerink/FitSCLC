function handles = SetActive(handles)
%{
Activates/de-activates the fields in the GUI according to the current
status and the options chosen.

Linköping University
Martijn Kemerink, November 20, 2018
%}

%% first: deactivate everything
set(handles.ProjectName,'Enable','off');
set(handles.mA_cm2,'Enable','off');
set(handles.LoadNewData,'Enable','off');
set(handles.SmoothTable,'Enable','off');
set(handles.Reset,'Enable','off');
set(handles.Smooth,'Enable','off');
set(handles.RangeMethod,'Enable','off');
set(handles.IrangeMin,'Enable','off');
set(handles.IrangeMax,'Enable','off');
set(handles.VrangeMin,'Enable','off');
set(handles.VrangeMax,'Enable','off');
set(handles.AutoRangeTable,'Enable','off');
set(handles.FindRange,'Enable','off');
set(handles.OutputTable,'Enable','off');
set(handles.RequestText,'Enable','off');
set(handles.FitModel,'Enable','off');
set(handles.muModel,'Enable','off');
set(handles.PositiveGamma,'Visible','off');
set(handles.ConstrainedFit,'Visible','off');
set(handles.generalFitTable,'Enable','off');
set(handles.GillTable,'Enable','off');
set(handles.GDMtable,'Enable','off');
set(handles.ArrheniusTable,'Enable','off');
set(handles.eGDMtable,'Enable','off');
set(handles.ET_GDMtable,'Enable','off');
set(handles.ShowGuess,'Enable','off');
set(handles.FitIt,'Enable','off');
set(handles.Save,'Enable','off');
set(handles.SpeedMode,'Visible','off');
set(handles.DDtable,'Enable','off');

%% second: enable what is relevant and/or possible
%always
set(handles.ProjectName,'Enable','on');
set(handles.mA_cm2,'Enable','on');
set(handles.LoadNewData,'Enable','on');
set(handles.OutputTable,'Enable','on');
%when data have been loaded
if handles.status.DataLoaded
    set(handles.SmoothTable,'Enable','on');
    set(handles.Reset,'Enable','on');
    set(handles.Smooth,'Enable','on');
    set(handles.RangeMethod,'Enable','on');
    switch handles.In.ranging
        case 1 %voltage range
            set(handles.VrangeMin,'Enable','on');
            set(handles.VrangeMax,'Enable','on');
        case 2 %current range
            set(handles.IrangeMin,'Enable','on');
            set(handles.IrangeMax,'Enable','on');
        otherwise %auto range
            set(handles.AutoRangeTable,'Enable','on');
    end
    set(handles.FindRange,'Enable','on');
end
%when points to fit have been identified
if handles.status.RangeSet
    set(handles.FitModel,'Enable','on');
    set(handles.generalFitTable,'Enable','on');
    set(handles.ShowGuess,'Enable','on');
    set(handles.FitIt,'Enable','on');
    if handles.In.FitModel==1 %Murgatroyd/Gill
        set(handles.muModel,'Enable','on');
        set(handles.PositiveGamma,'Visible','on');
        set(handles.ConstrainedFit,'Visible','on');
        set(handles.GillTable,'Enable','on');
        if handles.In.muModel==1 %GDM
            set(handles.GDMtable,'Enable','on');
        else %Arrhenius
            set(handles.ArrheniusTable,'Enable','on');
        end
    else %drift-diffusion/eGDM
        set(handles.muModel,'Enable','on');
        set(handles.eGDMtable,'Enable','on');
        if (handles.In.muModel==3)||(handles.In.muModel==4) %ET-GDM (lattice) or (random)
            set(handles.ET_GDMtable,'Enable','on');
        end
        set(handles.SpeedMode,'Visible','on');
        set(handles.DDtable,'Enable','on');
    end
end
%when a fit has been performed
if handles.status.Fitted
    set(handles.Save,'Enable','on');
    set(handles.RequestText,'Enable','on');
end

end
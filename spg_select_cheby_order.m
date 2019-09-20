%SPG_SELECT_CHEBY_ORDER Provides an interactive means to select an 
% appropriate Chebyshev polynomial order for estimating a set of 
% spectral kernels.
%
% The function is the code for the GUI spg_select_cheby_order.fig 
% SPG_SELECT_CHEBY_ORDER, by itself, creates a new 
% SPG_SELECT_CHEBY_ORDER or raises the existing singleton.
%
% Example:
%
% Selecting required Chbyshev polynomila order for approximating 
% the Meyer-like system of spectral kernels with 4 wavelet scales, 
% defined over the spectral range [0,2]:
%
% lmax = 2;
% g = spg_filter_design(lmax,5,'designtype','meyer');
% chebyOrder = spg_select_cheby_order(g,lmax); % close the GUI for to have output returned. 
%
% ------------------------------------------------------------
% Copyright (C) 2016, Hamid Behjat.
% This file is part of SPG (signal processing on graphs) package - v1.00
%
% Download: miplab.epfl.ch/software/
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = spg_select_cheby_order(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spg_select_cheby_order_OpeningFcn, ...
                   'gui_OutputFcn',  @spg_select_cheby_order_OutputFcn, ...
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

function spg_select_cheby_order_OpeningFcn(hObject, ~, handles, varargin)

maxOrder =200;
set(handles.s_chebyOrder,'Max',maxOrder);
set(handles.s_chebyOrder, 'Min',10);
set(handles.s_chebyOrder,'SliderStep',[1/maxOrder, 0.1000]);
set(handles.s_chebyOrder, 'Value',103);
set(handles.t_chebyOrder,'String',103);
handles.firstChebyChange = 1;
handles.chebyOrder = 103;

nChunks =10;
set(handles.s_zoom,'Max',nChunks);
set(handles.s_zoom, 'Min',1);
set(handles.s_zoom,'SliderStep',[1/nChunks, 0.1000]);
set(handles.s_zoom, 'Value',1);

if  ~isempty(varargin)
    handles.g = varargin{1};
    handles.arange(1) = 0;
    handles.arange(2) = varargin{2};
    
    if nargin==6
        if isstruct(varargin{3}) && isfield(varargin{3},'E')
            handles.G = varargin{3};
        end
    end
    
    % spectral kernels
    spg_view_design(handles.g,handles.arange,...
        'plotLineWidth',2,'guiHandle',handles.axes1);
    
    spg_view_design(handles.g,handles.arange,...
        'plotLineWidth',2,'guiHandle',handles.axes2);
    
    dummy = [0, handles.arange(2)];
    
    set(handles.axes1,'XLim',dummy,'YLim',[0,1.2],...
        'YTick',[0 1],'Box','off');
    
    set(handles.axes2,'XLim',[...
        dummy(2)/10*(get(handles.s_zoom, 'Value')-1),...
        dummy(2)/10*get(handles.s_zoom, 'Value')],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
    
    % Cheby. approx. of kernels
    spg_view_design(handles.g,handles.arange,...
        'plotLineWidth',2,'guiHandle',handles.axes3,...
        'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
    
    spg_view_design(handles.g,handles.arange,...
        'plotLineWidth',2,'guiHandle',handles.axes4,...
        'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
    
    set(handles.axes3,'XLim',dummy,'YLim',[0,1.2],...
        'YTick',[0 1],'Box','off');
    
    set(handles.axes4,'XLim',[...
        dummy(2)/10*(get(handles.s_zoom, 'Value')-1),...
        dummy(2)/10*get(handles.s_zoom, 'Value')],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
else
    error('Missing inputs; the function handles to the spectral kernels & the maximum spectral range need to specified as input to this function.')
end

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = spg_select_cheby_order_OutputFcn(~, ~, handles) 
varargout{1} = handles.chebyOrder;
delete(handles.figure1);

function s_chebyOrder_Callback(hObject, ~, handles)
if handles.firstChebyChange
    message = [...
        sprintf('Select an order, which makes the ''dashed'' line approximately straight. \n\n'),...
        sprintf('This ensures that the atoms that are constructed using the Chebyshev polynomial approximaton of the set of kernels form a ''tight'' frame.\n\n'),...
        sprintf('Close the GUI when you are satisfied with the order.')];
    uiwait(msgbox(message));
    pause(1)
    handles.firstChebyChange = 0;
end

handles.chebyOrder = round(get(hObject,'Value'));
set(handles.t_chebyOrder,'String',round(get(hObject,'Value')));
handles = updateFrame(handles,'chebyChange');
guidata(hObject, handles);

function s_chebyOrder_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function t_chebyOrder_Callback(hObject, ~, handles)
if handles.firstChebyChange
    message = [...
        sprintf('Select an order, which makes the ''dashed'' line approximately straight. \n\n'),...
        sprintf('This ensures that the atoms that are constructed using the Chebyshev polynomial approximaton of the set of kernels form a ''tight'' frame.')];
    uiwait(msgbox(message));
    pause(1)
    handles.firstChebyChange = 0;
end

handles.chebyOrder = str2double(get(hObject,'String'));
set(handles.s_chebyOrder,'Value',handles.chebyOrder);
handles = updateFrame(handles,'chebyChange');
guidata(hObject, handles);

function t_chebyOrder_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_zoom_Callback(hObject, ~, handles)
handles.zoom = round(get(hObject,'Value'));
handles = updateFrame(handles,'zoomChange');
guidata(hObject, handles);

function s_zoom_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function handles = updateFrame(handles,option)

if strcmp(option,'chebyChange')

    fprintf('Updating Chebyshev approx. plots for order: %d ...',...
        handles.chebyOrder)
    
    if ~isfield(handles,'G')
        
        spg_view_design(handles.g,handles.arange,...
            'plotLineWidth',2,'guiHandle',handles.axes3,...
            'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
        
        spg_view_design(handles.g,handles.arange,...
            'plotLineWidth',2,'guiHandle',handles.axes4,...
            'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
    else
        N_eigsToInterp = 1000;
        spg_view_design(handles.g,handles.arange,'G',handles.G,...
            'eigsToInterp',0:1/N_eigsToInterp:handles.arange(2),...
            'plotLineWidth',2,'guiHandle',handles.axes3,...
            'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
        
        spg_view_design(handles.g,handles.arange,'G',handles.G,...
            'eigsToInterp',0:1/N_eigsToInterp:handles.arange(2),...
            'plotLineWidth',2,'guiHandle',handles.axes4,...
            'chebyOrder',handles.chebyOrder,'onlyChebyApprox','yes');
    end
    
    dummy = [0, handles.arange(2)];
    
    set(handles.axes3,'XLim',dummy,'YLim',[0,1.2],...
        'YTick',[0 1],'Box','off');
    
    set(handles.axes4,'XLim',[...
        dummy(2)/10*(get(handles.s_zoom, 'Value')-1),...
        dummy(2)/10*get(handles.s_zoom, 'Value')],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
    fprintf(' done.\n\n')

elseif strcmp(option,'zoomChange')
    
    dummy = [0, handles.arange(2)];
    set(handles.axes2,'XLim',[...
        dummy(2)/10*(get(handles.s_zoom, 'Value')-1),...
        dummy(2)/10*get(handles.s_zoom, 'Value')],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
    
    set(handles.axes4,'XLim',[...
        dummy(2)/10*(get(handles.s_zoom, 'Value')-1),...
        dummy(2)/10*get(handles.s_zoom, 'Value')],...
        'YLim',[0,1.2],'YTick',[],'Box','off');
end

function figure1_CloseRequestFcn(hObject, ~, ~)
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
      delete(hObject);
end

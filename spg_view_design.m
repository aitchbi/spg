%SPG_VIEW_DESIGN Illustrate a system of spectral kernels.
% 
% This function is an extended version of the function:
%
% sgwt_filter_design 
%
% from the the SGWT toolbox.
%
% To see examples of usage and the required inputs, please refer to 
% either of the following demo functions:
%
% spg_demo_ (minnesota | alameda |cerebellum)
% spg_demo_uniformMeyerType
% spg_demo_ (your_data_frame | your_data_decompose)
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

function [f,x,maxY] = spg_view_design(g,arange,varargin)

control_params={'Graph',0,'eigsToInterp',[],'subSampleWarping',0,...
    'warping','none','plotLineWidth',1,'lambda',[],'guiHandle',[],...
    'chebyOrder',[],'onlyChebyApprox','no'};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if isstruct(Graph)
    theEigs = Graph.E;
    diffOfEigs = theEigs(2:end)-theEigs(1:end-1);
    indNonEqual = [1; find(logical(diffOfEigs)) + 1];
    
    if strcmp(warping,'none')
        if isempty(eigsToInterp)
            f = Graph.E(indNonEqual);
            x = f;
        else 
            f = eigsToInterp;
            x = f;
        end
    else
        if isempty(eigsToInterp)
            f = Graph.E(indNonEqual);
            x = warping(indNonEqual)*Graph.lmax;
        else
            if subSampleWarping
                dummy1 = Graph.E(indNonEqual);
                interp_x = dummy1(1:subSampleWarping:end-1);
                interp_x(end+1) = dummy1(end);
                
                dummy2 = warping(indNonEqual)*Graph.lmax;
                interp_y = dummy2(1:subSampleWarping:end-1);
                interp_y(end+1) = dummy2(end);
                
                f = eigsToInterp(eigsToInterp <= interp_x(end));
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f)); 
            else
                f = eigsToInterp;
                interp_x = Graph.E(indNonEqual);
                interp_y = warping(indNonEqual)*Graph.lmax;
                x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,f));
            end
        end
    end
    
elseif ~isempty(eigsToInterp) && strcmp(warping,'none')
    f = eigsToInterp;
    x = f;
else
    x=linspace(arange(1),arange(2),1e3);
    f = x;
end

J=numel(g)-1;
co=get(gca,'ColorOrder'); 
co=[co;co;co];
G=0*x;
G_cheby = 0*x;

if ~isempty(chebyOrder)
    for k=1:numel(g)
        c{k}=sgwt_cheby_coeff(g{k},chebyOrder,chebyOrder+1,arange);
    end
end

for n=0:J
    if isempty(guiHandle)
        if isempty(chebyOrder)
            plot(f,g{1+n}(x),'Color',co(1+n,:),'LineWidth',plotLineWidth);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'Color',co(1+n,:),...
                'LineWidth',plotLineWidth);
        else
            plot(f,g{1+n}(x),'Color',co(1+n,:),'LineWidth',plotLineWidth);
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'-.','Color',...
                co(1+n,:),'LineWidth',plotLineWidth);
        end
    else
        if isempty(chebyOrder)
            plot(f,g{1+n}(x),'Color',co(1+n,:),'LineWidth',plotLineWidth,...
                'Parent',guiHandle);
        elseif strcmp(onlyChebyApprox,'yes')
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'Color',co(1+n,:),...
                'LineWidth',plotLineWidth,'Parent',guiHandle);
        else
            plot(f,g{1+n}(x),'Color',co(1+n,:),'LineWidth',plotLineWidth,...
                'Parent',guiHandle);
            plot(f,sgwt_cheby_eval(x,c{1+n},arange),'-.','Color',...
                co(1+n,:),'LineWidth',plotLineWidth,'Parent',guiHandle);
        end
    end
    
    if n==0
        if isempty(guiHandle)
            hold on
        else
            hold(guiHandle,'on')
        end
    end
    
    if isempty(chebyOrder)
        G=G+g{1+n}(x).^2;
    else
        G=G+g{1+n}(x).^2;
        G_cheby = G_cheby+sgwt_cheby_eval(x,c{1+n},arange).^2;
    end
end

if isempty(guiHandle)
    if isempty(chebyOrder)
        plot(f,G,'k:','LineWidth',plotLineWidth);
    elseif strcmp(onlyChebyApprox,'yes')
        plot(f,G_cheby,'k:','LineWidth',plotLineWidth);
    else
        plot(f,G,'k:','LineWidth',plotLineWidth);
        plot(f,G_cheby,'k-.','LineWidth',plotLineWidth);
    end
else
    if isempty(chebyOrder)
        plot(f,G,'k:','LineWidth',plotLineWidth,'Parent',guiHandle);
    elseif strcmp(onlyChebyApprox,'yes')
        plot(f,G_cheby,'k:','LineWidth',plotLineWidth,'Parent',guiHandle);
    else
        plot(f,G,'k:','LineWidth',plotLineWidth,'Parent',guiHandle);
        plot(f,G_cheby,'k-.','LineWidth',plotLineWidth,'Parent',guiHandle);
    end
end

leglabels{1}='h';
for j=1:J
    leglabels{1+j}=sprintf('g_{%d}',j);
end
leglabels{J+2}='G';

if ~isempty(lambda)
    y=-.5*ones(size(lambda,1),1);
    stem(lambda,y, 'k');
    leglabels{J+3}='Lambda';
end

if isempty(guiHandle)
    hold off
else
    hold(guiHandle,'off')
end

maxY = max(G(:));

function hline(y,varargin)
xl=xlim;
plot(xl,y*[1 1],varargin{:});

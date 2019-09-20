% SPG_DEMO_CEREBELLUM Constructs signal-adapted system of of spectral
% kernels on the Cerebellum graph as presented in:
%
% Behjat, et al., 'Signal-adapted tight frames on graphs', 
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
%
%
% Examples: 
%
% 1. Results as shown in the above reference:   
% >> SPG_DEMO_CEREBELLUM
%
% 2. Changing the number of subbands from the default:
% >> SPG_DEMO_CEREBELLUM('N_subbands',6)
%
% 3. Changing the slice of the cerebellum graph for display:
% >> SPG_DEMO_CEREBELLUM('slNum',12)
%
% 4. Changing the number of subbands and the amount of 
%     smoothing of the warping function:
% >> SPG_DEMO_CEREBELLUM('N_subbands',3,'sFactor',25)
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

function spg_demo_cerebellum(varargin)

control_params={...
    'N_subbands',5,...
    'warpType','smoothed',...
    'slNum',15,...
    'sFactor',21,...
    'N_eigsToInterp',1000};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
fprintf('--------------------------------------------\n')

fprintf('Cerebellum graph. \n')

fprintf('Constructing graph... ')
G = spg_construct_cerebellum_graph;
fprintf(' done.\n')

pause(0.2)

fprintf('Loading graph signals... ')
dummy = load('signal_set_cerebellum.mat');
signalSets{1}.signals = dummy.F2; % F2 cerebellum signal set (cf. section V.C. paper [A])
fprintf(' done.\n')
pause(0.2)

G.L = sgwt_laplacian(G.A,'opt','normalized');

fprintf('Performing eigen-decomposition...')
[G.U,eigsFull] = eig(full(G.L));
fprintf(' done.\n')
G.E    = diag(eigsFull);
G.lmax = max(eigsFull(:));

% construct warpings ------------------------------------------------------
fprintf('Constructing warpings...\n')
[signalSets, G] = spg_construct_warping(signalSets,G,'sFactor',sFactor);
fprintf('Energy-equalizing warping constructed.\n')

noWarp = G.E/(G.E(end));
warpL = [0:1/(G.N -1):1]';

% plot warpings -----------------------------------------------------------
hFig = figure('Name','Spectral Transformations','NumberTitle','off');
set(hFig, 'Position', [715, 475, 320, 260]);
hold(gca,'on')
plot(G.E,noWarp*G.lmax,'--r','LineWidth',2)
plot(G.E,warpL*G.lmax, 'Color',[0.5 0.5 0.5],'LineWidth',2)

switch warpType
    case 'original'
        warpE1 = signalSets{1}.warpE;
    case 'smoothed',
        warpE1 = signalSets{1}.warpE_smooth;
end
plot(G.E,warpE1*G.lmax, 'k','LineWidth',2)

switch warpType
    case 'original'
        legend('noWarp','spectrum-adapted warping',...
            'EE transformation',...
            'Location','southeast')
    case 'smoothed',
        legend('noWarp','spectrum-adapted warping',...
            'smoothed EE transformation',...
            'Location','southeast')
end
title('Spectral Transformations')
set(gca,'XLim',[0 G.lmax],...
    'YLim',[0 G.lmax],'XGrid','on','YGrid','on')
hold(gca,'off')

% plot distribution of eigenvalues 
% -------------------------------------------------------------------------
hFig = figure('Name','Distribution of eigenvalues','NumberTitle','off');
set(hFig, 'Position', [715, 260, 320, 140]);
N_bins =20;
[counts, centers] = hist(G.E,N_bins);
bar(centers,counts,'k')
set(gca,'Box','off')


% plot distribution of ensemble energy 
% -------------------------------------------------------------------------
hFig = figure('Name','Distribution of ensemble energies- F1','NumberTitle','off');
set(hFig, 'Position', [715, 10, 320, 140]);

Es = signalSets{1}.En(2:end);
dummy = zeros(1,N_bins);
for i =1:N_bins
    ind1 = find(G.E(2:end) <= (G.lmax/N_bins)*i);
    ind2 = find(G.E(2:end) > (G.lmax/N_bins)*(i-1));
    ind = find(ismember(ind1,ind2));
    dummy(i) = sum(Es(ind));
end
bar(centers,dummy,'k')
set(gca,'Box','off')

% plots graph -------------------------------------------------------------
G_ss = spg_construct_cerebellum_graph_ss(slNum);
gsp_plot_graph(G_ss)
set(gcf,'Position', [1040, 10, 350, 200],...
    'Name',sprintf('Cerebellum graph, layer %d (zoom to explore)',slNum),...
    'NumberTitle','off')

% construct spectral kernels ----------------------------------------------

% SGWT spectral kernels
g1 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','abspline3');

% Meyer-like spectral kernels
g2 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','meyer');

% Half Cosine Uniform Translates Frame
g3 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','half_cosine_uniform_translates');

% Spectrum-adapted spectral kernels
g4 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','spectrum_adapted','G',G);

% Meyer-type Uniform spectral kernels
g5 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','uniform_meyer_type');

%  Signal-adapted spectral kernels
switch warpType
    case 'original', dummy = signalSets{1}.warpE;
    case 'smoothed', dummy = signalSets{1}.warpE_smooth;
end

g6 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','signal_adapted','warping',...
    dummy,'E',G.E);



% plot spectral kernels ---------------------------------------------------
hFig = figure('Name','Systems of Spectral Kernels Forming Frames',...
    'NumberTitle','off');
set(hFig, 'Position', [10, 10, 700, 105*(6)]);
lineWidth =2;

% Spline-based system of spectral kernels
subplot(6,1,1)
[~,~,yLim] = spg_view_design(g1,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('SGWT system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2*yLim])

% Meyer-like Wavelet system of spectral kernels
subplot(6,1,2)
spg_view_design(g2,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Meyer-like wavelet system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Half Cosine Uniform Translates sytem of spectral kernels
subplot(6,1,3)
spg_view_design(g3,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Half cosine uniform translates system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Spectrum-adapted system of spectral kernels
subplot(6,1,4)
spg_view_design(g4,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Spectrum-adapted system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Meyer-type Uniform system of spectral kernels
subplot(6,1,5)
spg_view_design(g5,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Uniform Meyer-type system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Signal-adapted system of spectral kernels
subplot(6,1,6)
spg_view_design(g6,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
   'plotLineWidth',lineWidth);
title('Signal-adapted system of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

if ~nargin
    fprintf('--------------------------------------------\n')
    message = [...
        sprintf('Cerebellum graph.\n\n'),...
        sprintf('Systems of spectral kernels with 5 subbands are'),...
        sprintf(' displayed as presented in: \n\n'),...
        sprintf('Behjat, et al., ''Signal-adapted tight frames on graphs'','),...
        sprintf('\nIEEE Trans. Signal Process., 2016, doi: 10.1109/TSP.2016.2591513 \n\n'),...
        sprintf('You may specify a different number of subbands by changing'),...
        sprintf(' the ''N_subbands'' parameter when running the demo;'),...
        sprintf(' for instance:\n\n'),...
        sprintf('>> spg_demo_minnesota(''N_subbands'',6) \n\n'),...
        ];
    uiwait(msgbox(message));
    
    message = [...
        sprintf('The cerebellum graph is defined within a volume, i.e.,'),...
        sprintf(' 3D space. Here, we only display one of its layers, with'),...
        sprintf(' its vertices and edges.\n\n'),...
        sprintf('You can view the other layers by changing the ''slNum'' '),...
        sprintf(' parameter when running the demo:\n\n'),...
        sprintf('>> spg_demo_cerebellum(''slNum'',xx) \n\n'),...
        sprintf('where xx can be a value between 1-32.\n\n'),...
        sprintf('Note that in this way of illustrating the graph, edges'),...
        sprintf(' connecting vertices in adjacent layers are not displayed.\n\n')];
    uiwait(msgbox(message));
end

end
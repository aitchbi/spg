% SPG_DEMO_YOUR_DATA_FRAME Construct signal-adapted system of 
% spectral kernels for a given graph and signal set. 
%
% For comparison, the SGWT frame, Meyer-like and spectrum-adapted 
% systems of spectral kernels will also be constructed. 
% Distributions of eigenvalues and ensemble signal energies will 
% also be displayed.
%
% [g,S,G] = SPG_DEMO_YOUR_DATA(A,S,LType)
%
% Inputs:
%
% A     : Graph adjacency matrix
%
% S     : Graph signal set, consisting of M signals, each of 
%          dimension N*1, formed in an N*M matrix. This signal set is used
%          for constructing the energy-equalizing transformation (warping),
%          and thus, to construct the signal-adapted frame.
%
% LType : Laplacian type. 'normalized' | 'combinatorial'
%
%
% Optional input parameters:
%
% 'N_subbands' : number of subbands [default: 7]
%               
%                For example, to set the number of subbands to 5:
%                >> SPG_DEMO_YOUR_DATA(A,S,LType,'N_subbands',5)
%
% 'warpType'   : the type of energy-equalizing transformation to use.
%                'original' | 'smoothed' [default: 'original']
%  
%                Depending on the dataset, the smoothed version can 
%                be beneficial to obtain smoother spectral kernels.
%                using the smoothed version, a moving average will be
%                implemented across the ensemble energies, using a 
%                boxcar window.
% 
%                To use the smoothed version:
%                >> SPG_DEMO_YOUR_DATA(A,S,LType,'warpType','smoothed')
%
% 'sFactor'    : the length of boxcar window used for smoothing.  
%
%                For example, to set the length to 10:
%                >> SPG_DEMO_YOUR_DATA(A,S,LType,'sFactor',10)
%
% Outputs:
% g                 : system of signal-adapted spectral kernels, where 
%                      g{1} : function handle for first subband kernel
%                      g{2} : function handle for second subband kernel
%                       ...
%                      g{N_subbands} : function handle for last subband kernel 
%
% S                 : structure, with fields:
%                      - signals      : input signal set S
%                      - En             : ensemble energy spectral density of S
%                      - C_energy   : cumulative function of En
%                      - warpE        : energy equalizing spectral transformation
%                      - warpE_smooth : smoothed version of warpE, using moving
%                                                    average with a boxcar of lenbth sFactor.
%
% G              : structure containg info about the constructed graph:
%                  - A        : Adjacency matrix
%                  - N        : Number of graph vertices
%                  - lap_type : Laplacian type used to define spectrum
%                  - L        : Laplacian matrix
%                  - U        : Laplacian matrix Eigenvectors, one in each column
%                  - E        : Laplacian matrix eigenvalues
%                  - lmax     : maximum eigenvalue
%                  - repetitiveEigs       : repetitive eigenvalues
%                  - repetitiveEigsCounts : multiplicitity of repetitive eigenvalues
%
% Example:
%
% load('A_minnesota.mat')
% load('signal_sets_minnesota.mat') % consists of 2 signal sets: F1 & F2
% g = SPG_DEMO_YOUR_DATA_FRAME(A,F1,'normalized');
%
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

function [S,G] = spg_demo_your_data_frame(A,S,LType,varargin)

control_params={'N_subbands',7,'warpType','original','sFactor',21,...
    'N_eigsToInterp',1000};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
fprintf('-----------------------------------\n')

G.A = A;
G.N = size(A,2);

if strcmp(LType,'normalized')
    G.lap_type = LType;
    G.L = sgwt_laplacian(G.A,'opt',LType);
elseif strcmp(LType,'combinatorial')
    G.lap_type = LType;
    G.L = sgwt_laplacian(G.A,'opt','raw');
else
    error('unknown Laplacian type.')
end

fprintf('Performing eigen-decomposition...')
[G.U,eigsFull] = eig(full(G.L));
fprintf(' done.\n')
G.E    = diag(eigsFull);
G.lmax = max(eigsFull(:));

if G.N ~= size(S,1)
    error('The signals must have a length equal to the number of graph vertices.')
end

% construct warpings --------------------------
signalSets{1}.signals = S;
fprintf('Constructing warpings...\n')
[signalSets, G] = spg_construct_warping(signalSets,G,'sFactor',sFactor);
fprintf('Energy-equalizing warping constructed.\n')

noWarp = G.E/(G.E(end));
warpL = [0:1/(G.N -1):1]';

% plot warpings -------------------------------
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
            'Energy Equalizing transformation',...
            'Location','southeast')
    case 'smoothed',
        legend('noWarp','spectrum-adapted warping',...
            'Energy Equalizing transformation',...
            'Location','southeast')
end
title('Spectral Transformations')
set(gca,'XLim',[0 G.lmax],...
    'YLim',[0 G.lmax],'XGrid','on','YGrid','on')
hold(gca,'off')

% plot distribution of eigenvalues 
% -------------------------------------------
hFig = figure('Name','Distribution of eigenvalues','NumberTitle','off');
set(hFig, 'Position', [715, 260, 320, 140]);
N_bins = 20;
[counts, centers] = hist(G.E,N_bins);
bar(centers,counts,'k')
set(gca,'Box','off')


% plot distribution of ensemble energy 
% -------------------------------------------
hFig = figure('Name','Distribution of ensemble energies','NumberTitle','off');
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

% construct spectral kernels ----------------
% -------------------------------------------

% SGWT spectral kernels
g1 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','abspline3');

% Meyer-like spectral kernels
g2 = spg_filter_design(G.lmax,N_subbands,...
    'designtype','meyer');

% Half Cosine Uniform Translates sytem of spectral kernels
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
g = spg_filter_design(G.lmax,N_subbands,...
    'designtype','signal_adapted','warping',...
    dummy,'E',G.E);

% plot spectral kernels ---------------------
% -------------------------------------------

hFig = figure('Name','Systems of Spectral Kernels Forming Frames','NumberTitle','off');
set(hFig, 'Position', [10, 10, 700, 105*(6)]);
lineWidth =2;

% SGWT frame
subplot(6,1,1)
[~,~,yLim] = spg_view_design(g1,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('SGWT sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2*yLim])

% Meyer-like Wavelet sytem of spectral kernels
subplot(6,1,2)
spg_view_design(g2,[0,G.lmax],'Graph',G,...
   'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
   'plotLineWidth',lineWidth);
title('Meyer-like wavelet sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Half Cosine Uniform Translates sytem of spectral kernels
subplot(6,1,3)
spg_view_design(g3,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Half cosine uniform translates sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Spectrum-adapted sytem of spectral kernels
subplot(6,1,4)
spg_view_design(g4,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Spectrum-adapted sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Meyer-type Uniform frame
subplot(6,1,5)
spg_view_design(g5,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Uniform Meyer-type sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])

% Signal-adapted sytem of spectral kernels
subplot(6,1,6)
spg_view_design(g,[0,G.lmax],'Graph',G,...
    'eigsToInterp',0:1/N_eigsToInterp:G.lmax,...
    'plotLineWidth',lineWidth);
title('Signal-adapted sytem of spectral kernels')
set(gca,'Box','off','XLim',[0 G.lmax],'YLim',[0 1.2])


S = signalSets;
end
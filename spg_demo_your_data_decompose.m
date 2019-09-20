%SPG_DEMO_YOUR_DATA_DECOMPOSE Construct signal-adapted system 
% of spectral kernels for a given graph and signal set, and decompose 
% the signal set. 
% 
% SPG_DEMO_YOUR_DATA_DECOMPOSE(A,S,LType)
%
% Inputs:
%
% A     : Graph adjacency matrix
%
% S1   : Graph signal set, consisting of M signals, each of 
%          dimension N*1, formed in an N*M matrix. This signal set is used
%          for constructing the energy-equalizing transformation (warping),
%          and thus, to construct the signal-adapted frame.
%
% S2   : The signal set to be decomposed. consisting of X signals, each of 
%          dimension N*1, formed in an N*X matrix.
%
% LType : Laplacian type. 'normalized' | 'combinatorial'
%
%
% Optional input parameters:
%
% 'N_subbands'   : number of subbands [default: 7]
%               
%                  For example, to set the number of subbands to 5:
%                  >> SPG_DEMO_YOUR_DATA_DECOMPOSE(A,S,LType,'N_subbands',5)
%
% 'decompMethod' : Method used for decomposition.
%                 'direct' | 'chebyApprox'  [default: 'direct']
%
%                  To use Chebyshev approximation scheme:
%                  >> SPG_DEMO_YOUR_DATA_DECOMPOSE(A,S,LType,'decompMethod','chebyApprox')
%                  
% 'warpType'     : the type of energy-equalizing transformation to use.
%                  'original' | 'smoothed' [default: 'original']
%  
%                  Depending on the dataset, the smoothed version can 
%                  be beneficial to obtain smoother spectral kernels.
%                  Using the smoothed version, a moving average will be
%                  implemented across the ensemble energies, using a 
%                  boxcar window.
% 
%                  To use the smoothed version:
%                  >> SPG_DEMO_YOUR_DATA_DECOMPOSE(A,S,LType,'warpType','smoothed')
%
% 'sFactor'      : the length of boxcar window used for smoothing.  
%
%                  For example, to set the length to 10:
%                  >> SPG_DEMO_YOUR_DATA_DECOMPOSE(A,S,LType,'sFactor',10)
%
% Outputs:
%
% S2_decomp : a 1*N_subbands cell array. Each cell j contains an N*X
%                  matrix, where each column i represents the j-th subband 
%                  decomposition coeffcients of signal i (i.e. the signal 
%                  located in i-th column of S2). 
% 
% S1             : structure, with fields:
%                  - signals      : input signal set S1
%                  - En             : ensemble energy spectral density of S1
%                  - C_energy   : cumulative function of En
%                  - warpE        : energy equalizing spectral transformation
%                  - warpE_smooth : smoothed version of warpE, using moving
%                                   average with a boxcar of lenbth sFactor.
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
%
% Examples:
%
% 1. Construct signal-adapted system of spectral kernels based on the F1 
% Minnesota signal set (as described in Behjat, et. al, 'Signal-adapted Tight 
% Frames on Graphs', IEEE TSP, 2016) and then decompose the signal set 
% using the direct approach and the Chbyshev polynomila approximation 
% method:
%
% load('A_minnesota.mat')
% load('signal_sets_minnesota.mat') % consists of 2 signal sets: F1 & F2
% S_decomp_direct = SPG_DEMO_YOUR_DATA_DECOMPOSE(A,F1,F1,'normalized','decompMethod','direct');
% S_decomp_cheby = SPG_DEMO_YOUR_DATA_DECOMPOSE(A,F1,F1,'normalized','decompMethod','chebyApprox');
%
% 2. Construct signal-adapted system of spectral kernels based on the F1 
% Minnesota signal set and then decompose the F2 Minnesota signal set:
%
% load('A_minnesota.mat')
% load('signal_sets_minnesota.mat') % consists of 2 signal sets: F1 & F2
% S_decomp_direct = SPG_DEMO_YOUR_DATA_DECOMPOSE(A,F1,F2,'normalized','decompMethod','direct');
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

function [S2_decomp,S1,G] = spg_demo_your_data_decompose(A,S1,S2,LType,varargin)

control_params={'N_subbands',6,'decompMethod','direct','warpType','original','sFactor',21,...
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

[G.U,eigsFull] = eig(full(G.L));
G.E    = diag(eigsFull);
G.lmax = max(eigsFull(:));

if G.N ~= size(S1,1) || G.N ~= size(S2,1) 
    error('The signals must have a length equal to the number of graph vertices.')
end

% construct warpings ------------------------------------------------------
signalSets{1}.signals = S1;
[signalSets, G] = spg_construct_warping(signalSets,G,'sFactor',sFactor);

% select warping type -----------------------------------------------------
switch warpType
    case 'original', dummy = signalSets{1}.warpE;
    case 'smoothed', dummy = signalSets{1}.warpE_smooth;
end

% construct signal-adapted spectral kernels (function handles) ------------
g = spg_filter_design(G.lmax,N_subbands,...
    'designtype','signal_adapted','warping',...
    dummy,'E',G.E);

% decompose signals -------------------------------------------------------
if strcmp(decompMethod,'direct') 
    % Direct spectral graph decomposition
    fprintf('Deco')
    tic
    for iS=1:size(S2,2)
        fhat=G.U'*S2(:,iS);
        for iK=1:numel(g)
            S2_decomp{iK}(:,iS)=G.U*(fhat.*g{iK}(G.E));
        end
    end
    
elseif strcmp(decompMethod,'chebyApprox') 
    % Decomposition using Chebyshev polynomial approximation
    
    % Select a suitable order for Chebyshev polynomial approximation
    chebyOrder = spg_select_cheby_order(g,G.lmax,G);
    
    tic
    % Compute Chebyshev coefficients
    for k=1:numel(g)
        c{k}=sgwt_cheby_coeff(g{k},chebyOrder,chebyOrder+1,[0,G.lmax]);
    end
    
    for iS=1:size(S2,2)
        dummy = sgwt_cheby_op(S2(:,iS),G.L,c,[0,G.lmax]);
        for iK=1:numel(g)
            S2_decomp{iK}(:,iS) = dummy{iK};
        end
    end
end
t = toc; % to compare computaitional costsof direct & Cheby. approx. decompositions  
S1 = signalSets;

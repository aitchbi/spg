% SPG_FILTER_DESIGN Construct system of spectral kernels for a given graph
%                               (and graph signal set).
%
% Baed on the provided inputs, SPG_FILTER_DESIGN constructs:
%
% (1) Uniform Meyer-type system of spectral kernels    
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','signal_adapted','warping',warping,'E',E)
% 
% &
%
% (2) Signal-adapted system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','uniform_meyer_type')
%
% as proposed in:
%
%                H. Behjat, et al., ''Signal-adapted tight frames on graphs'', 
%                IEEE Trans. Signal Process., 2016, 
%                doi: 10.1109/TSP.2016.2591513.
%
% (3) Spline-based system of spectral kernels (SGWT frame)
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','abspline3')
%
% as proposed in:
%
%                D.K. Hammond, et al., ''Wavelets on graphs via spectral graph
%                theory'', Appl. Comput. Harmon. Anal., vol. 30, 
%                pp. 129-150, 2011.   
%
% (4) Meyer-like system of spectral kernels 
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','meyer')
%
% as proposed in:
%
%                N. Leonardi, et al., ''Tight wavelet frames on multislice graphs'', 
%                IEEE Trans. Signal Process., vol. 61(13), pp. 3357-3367, 2013.  
%
% (5) Half-cosine unform translates system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','half_cosine_uniform_translates')
%
% &
%
% (6) Sectrum-adapted system of spectral kernels
% >> g = SPG_FILTER_DESIGN(lmax,N_subbands, 'designtype','spectrum_adapted','G',G)
% 
% as proposed in:
%
%                D. I. Shuman, et al., ''Spectrum-adapted tight graph wavelet
%                and vertex-frequency frames'', IEEE Trans. Signal Process., 
%                vol. 63(16), pp. 4223-4235, 2015.
%
% ------------------------------------------------------------
%
% Inputs:
%
% lmax             : upper bound on spectrum
%
% N_subbands   : number of subbands [default: 7]
%
%
% Optional input parameters :
%
% 'designtype'      : Type of system of spectral kernels to construct 
%                          [default: 'signal_adapted']
%
% 'warping'          : Energy-equalizing transformation used for warping 
%                          the uniform Meyer-type system of spectral kernels
%                          to construct a signal-adapted syetem of spectral 
%                          kernels.
%                        
% 'G'                   : graph structure containg at least the following 
%                         required fields:
%                         - N        : Number of graph vertices
%                         - L        : Graph Laplacian matrix
%                      
% 'E'                   : eigenvalues of the graph Laplacian matrix 
%
% 'subSampleWarping' : [default: 0]
%                                 If nonzero, constructs energy-equalizing transformation  
%                                 based on every N-th eigenvalue, 
%                                 where subSampleWarping = N. 
%
%
% Outputs:
%
% g{1}                 : function handle for first subband spectral kernel
% g{2}                 : function handle for second subband spectral kernel
% ...
% g{N_subbands} : function handle for last subband spectral kernel 
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

function g = spg_filter_design(lmax,N_subbands,varargin)
control_params = {'designtype','signal_adapted','warping',[],'G',[],...
    'E',[],'subSampleWarping',0};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
N_subbands = N_subbands-1; 
g = cell(N_subbands+1,1);

switch designtype
    
    case 'uniform_meyer_type'
        % Uniform Meyer-type system of spectral kernels
        
        g{1} = @(x) spg_kernel_1(x,lmax,N_subbands);
        for j = 1:N_subbands
            g{j+1} = @(x) spg_kernel_2(x,j,lmax,N_subbands);
        end
        
    case 'signal_adapted'
        % Signal-adapted system of spectral kernels
        
        if subSampleWarping 
            % construct warping based on every N-th lambda;  
            % where N = subSampleWarping .
            g{1} = @(x) spg_kernel_1(x,lmax,N_subbands,...
                'warping',warping,'E',E,...
                'subSampleWarping',subSampleWarping);
            
            for j = 1:N_subbands
                g{j+1} = @(x) spg_kernel_2(x,j,lmax,N_subbands,...
                    'warping',warping,'E',E,...
                    'subSampleWarping',subSampleWarping);
            end
        else
            g{1} = @(x) spg_kernel_1(x,lmax,N_subbands,...
                'warping',warping,'E',E);
            for j = 1:N_subbands
                g{j+1} = @(x) spg_kernel_2(x,j,lmax,N_subbands,...
                    'warping',warping,'E',E);
            end
        end
        
    case 'abspline3'
        % Spline-based system of spectral kernels (SGWT frame)
        
        g = sgwt_filter_design(lmax,N_subbands,'designtype','abspline3');
        
    case 'meyer'
        % Meyer-like system of spectral kernels 
        
        gb = @(x) sgwt_meyer(x);
        gb2 = @(x) sgwt_mey_h(x);
        gb3 = @(x) sgwt_meyer_end(x);
        
        t = sgwt_setscales_mey(lmax,N_subbands);
        
        % wavelet kernels
        for j = 1:N_subbands-1
            g{j+1} = @(x) gb(t(j)*x);
        end
        g{N_subbands+1} = @(x) gb3(t(end)*x);
        
        % sacling function kernel
        g{1} = @(x) gb2(t(1)*x);
        
        
    case 'half_cosine_uniform_translates'
        % Half-cosine unform translates system of spectral kernels 
        
        g = gsp_filter_design('uniform_translates',N_subbands+1,lmax);
        
    case 'spectrum_adapted'
        % Spectrum-adapted system of spectral kernels
        
        G = gsp_spectrum_cdf_approx(G);
        wParam.warp_function_type = 'custom';
        wParam.warp_function = G.spectral_warp_fn;
        wParam.upper_bound_translates = 1;
        
        g = gsp_filter_design('warped_translates',N_subbands+1,lmax,wParam);
        
        otherwise
        error('Unknown design type');
end
function [filters,varargout]=gsp_filter_design(filter_type,num_filters,lmax,varargin)
%GSP_FILTER_DESIGN  Design a system of graph spectral filters.
%   Usage:  filters=gsp_filter_design(filter_type,num_filters,lmax);
%           [filters,warp_fn]=gsp_filter_design(filter_type,num_filters,lmax);
%           filters=gsp_filter_design(filter_type,num_filters,lmax,param);
%           [filters,warp_fn]=gsp_filter_design(filter_type,num_filters,lmax,param);
%
%   Input parameters:
%         filter_type           : String containing the type of filter design.
%         num_filters           : Number of filters to cover the interval [0,lmax]
%         lmax                  : Upper bound on the Laplacian spectrum
%   Output parameters:
%         filters               : A cell array with each cell a function handle to the respective filter.
%         varargout             : varargout{1} returns the warping function if a warping function is used.
%   Additional parameters:
%
%   'gsp_filter_design(filter_type,num_filters,lmax)' computes a set of
%   graph spectral filters.
%
%   See also:  
%
%   Demos:  
%
%   Requires: 
% 
%   References:
%     D. K. Hammond, P. Vandergheynst, and R. Gribonval. Wavelets on graphs
%     via spectral graph theory. Appl. Comput. Harmon. Anal., 30(2):129-150,
%     Mar. 2011.
%     
%     N. Leonardi and D. Van De Ville. Tight wavelet frames on multislice
%     graphs. IEEE Trans. Signal Process., 61(13):3357-3367, Jul. 2013.
%     
%
%
%   Url: http://gspbox.sourceforge.net/doc/filter_design/gsp_filter_design.php

% Copyright (C) 2013 David I Shuman, Nathanael Perraudin.
% This file is part of GSPbox version 0.0.1
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

%   AUTHOR : David I Shuman.
%   TESTING: 
%   REFERENCE:
  
if ~isempty(varargin)
    param=varargin{1};
else
    param=0;
end
varargout{1}=0;

switch filter_type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Section 1: Filters that are only adapted to the length of the spectrum
    
    %%%%%%%%%%%%%%%%%
    % The spectral graph wavelet transform (Hammond et al., 2011),
    % implemented by the SGWT Toolbox
    case 'sgwt'
        if param==0
            [filters,~]=sgwt_filter_design(lmax,num_filters-1);
        else
            [filters,~]=sgwt_filter_design(lmax,num_filters-1,param);
        end
        
    %%%%%%%%%%%%%%%%%
    % Variations of the SGWT, also implemented in the SGWT Toolbox   
    % sgwt simple tight frame
    case 'sgwt_simple_tf'
        designtype='simple_tf';
        [filters,~]=sgwt_filter_design(lmax,num_filters-1,'designtype',designtype);

    % sgwt Meyer frame
    case 'sgwt_meyer'
        designtype='meyer';
        [filters,~]=sgwt_filter_design(lmax,num_filters-1,'designtype',designtype);
    
    % sgwt Mexican Hat frame
    case 'sgwt_mexican_hat'
        designtype='mexican_hat';
        [filters,~]=sgwt_filter_design(lmax,num_filters-1,'designtype',designtype);
       
    %%%%%%%%%%%%%%%%%
    case 'parseval'
        if num_filters>2
            error('At most two kernels are allowed');
        end
        
        if ~isfield(param, 'kernel_type')
            kernel_type = 'meyer';
        else
            kernel_type = param.kernel_type;
        end
        
        switch kernel_type
            
            case 'held'
                % Lowpass Kernel
                filters{1}= @(x) gsp_held_kernel(x*(2/lmax));

                % Highpass Kernel
                if num_filters>1
                    filters{2} = @(x) gsp_held_kernel(x*(2/lmax),1);
                end
                
            case 'simoncelli'
                % Lowpass Kernel
                filters{1}= @(x) gsp_simoncelli_kernel(x*(2/lmax));

                % Highpass Kernel
                if num_filters>1
                    filters{2} = @(x) gsp_simoncelli_kernel(x*(2/lmax),1);
                end
                
            case 'papadakis'
                % Lowpass Kernel
                filters{1}= @(x) gsp_papadakis_kernel(x*(2/lmax));

                % Highpass Kernel
                if num_filters>1
                    filters{2} = @(x) gsp_papadakis_kernel(x*(2/lmax),1);
                end
                
            case 'regular'
                if ~isfield(param, 'degree')
                    degree = 3;
                else
                    degree = param.degree;
                end                
                
                % Lowpass Kernel
                filters{1}= @(x) gsp_regular_kernel(x*(2/lmax),degree);

                % Highpass Kernel
                if num_filters>1
                    filters{2} = @(x) gsp_regular_kernel(x*(2/lmax),degree,1);
                end
                                
            otherwise % default is 'meyer'
                if num_filters>2
                    error('At most two Meyer kernels are allowed');
                end

                % Lowpass Meyer Kernel
                filters{1}= @(x) sgwt_kernel_meyer(x*(2/lmax),'sf');

                % Highpass Meyer Kernel
                if num_filters>1
                    filters{2} = @(x) real(sqrt(1-(sgwt_kernel_meyer(x*(2/lmax),'sf')).^2));
                end
        end
        
    %%%%%%%%%%%%%%%%%
    % Uniform_translates of a half cosine window
    case 'uniform_translates'
    filters = gsp_uniform_translates_filter_design(num_filters,lmax);

    %%%%%%%%%%%%%%%%%    
    % Biorthogonal filter banks of Narang and Ortega (still need to
    % implement)
    % case 'biorthogonal'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Section 2: Filter banks that are adapted to the maximum degree of the graph 
% and (optionally) the length of the spectrum

    % sgwt Meyer frame adapted to degree (Leonardi and Van De Ville, 2013)
    case 'sgwt_meyer_degree'
        if ~isfield(param, 'max_degree')
            error('Must specify the maximum degree of the graph');
        else
            max_degree = param.max_degree;
        end
        if ~isfield(param, 'adapt_to_spectrum_length')
            adapt_to_spectrum_length=1;
        else
            adapt_to_spectrum_length = param.adapt_to_spectrum_length;
        end
        designtype='meyer';
        filters=cell(num_filters,1);

        if adapt_to_spectrum_length
            [outer_filters,~]=sgwt_filter_design(lmax,num_filters-1,'designtype',designtype);
             for i=1:num_filters
                filters{i}=@(x) outer_filters{i}(lmax/acos(1-lmax/max_degree)*acos(1-x/max_degree)); 
             end
        else
            [outer_filters,~]=sgwt_filter_design(pi,num_filters-1,'designtype',designtype);
            for i=1:num_filters
                filters{i}=@(x) outer_filters{i}(acos(1-x/max_degree)); 
            end
        end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Section 3: Filter banks that are adapted to the specific spectrum or a subset 
% approximation of the spectrum
    
    %%%%%%%%%%%%%%%%%     
    % Lapped Kernels
    case 'lapped_kernels'
        if ~isfield(param,'approx_spectrum')
            error('Approximate spectrum is not specified')
        end
        if ~isfield(param, 'kernel_type')
            kernel_type = 'held';
        else
            kernel_type = param.kernel_type;
        end
        switch kernel_type
           
            case 'held'
                low_kernel=@(x) gsp_held_kernel(x);
                high_kernel=@(x) gsp_held_kernel(x,1);
                
            case 'simoncelli'
                low_kernel=@(x) gsp_simoncelli_kernel(x);
                high_kernel=@(x) gsp_simoncelli_kernel(x,1);
                
            case 'papadakis'
                low_kernel=@(x) gsp_papadakis_kernel(x);
                high_kernel=@(x) gsp_papadakis_kernel(x,1);
 
            case 'regular'
                if ~isfield(param, 'degree')
                    degree = 3;
                else
                    degree = param.degree;
                end
                low_kernel=@(x) gsp_regular_kernel(x,degree);
                high_kernel=@(x) gsp_regular_kernel(x,degree,1);
            
            otherwise % default is 'meyer'
                low_kernel=@(x) gsp_meyer_kernel(x);
                high_kernel=@(x) gsp_meyer_kernel(x,1);
        end
        midpoints=param.approx_spectrum;
        filters=cell(num_filters,1);
        midpoints(1)=midpoints(1)-1e-8;
        filters{1}=@(x) (x>=midpoints(1)).*(x<midpoints(2)).*low_kernel(2*(x-midpoints(1))/(midpoints(2)-midpoints(1)));
        for k=2:num_filters-1
           filters{k}=@(x) (x>=midpoints(k-1)).*(x<midpoints(k)).*high_kernel(2*(x-midpoints(k-1))/(midpoints(k)-midpoints(k-1)))+(x>=midpoints(k)).*(x<midpoints(k+1)).*low_kernel(2*(x-midpoints(k))/(midpoints(k+1)-midpoints(k)));
        end
        filters{num_filters}=@(x) (x>=midpoints(num_filters-1)).*(x<midpoints(num_filters)).*high_kernel(2*(x-midpoints(num_filters-1))/(midpoints(num_filters)-midpoints(num_filters-1)));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Section 4: Filters that are warped versions of the uniform half cosine 
% translates described above.  
    
    case 'warped_translates'
        if ~isfield(param, 'warp_function_type')
            warp_function_type = 'log';
        else
            warp_function_type = param.warp_function_type;
        end
        
        switch warp_function_type
            
            %%%%%%%%%%%%%%%%%
            % Log warping functinon. An alterntive way to construct 
            % spectral graph wavelets. Comparing to SGWT, the support of 
            % the higher scale wavelet starts at higher spectral values. 
            case 'log'
                filters=cell(num_filters,1);
                warp_function = @(s) log(s+eps);
                
                % Generate (num_filters-1) uniform translates covering [0,log(lmax)] 
                uniform_filters = gsp_uniform_translates_filter_design(num_filters-1,log(lmax));
                
                % Warp the uniform translates by log to generate the "wavelet" filters
                wavelet_filters=cell(num_filters-1,1);
                for i=1:num_filters-1
                    wavelet_filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
                end
                filters(2:num_filters)=wavelet_filters;
                
                % Generate a "scaling" filter that results in a tight frame
                filters{1}= @(x) gsp_log_wavelet_scaling_fn(wavelet_filters,x);
                
            %%%%%%%%%%%%%%%%%
            % Composition of log warping and a custom warping functinon. 
            % An alterntive way to construct spectral graph wavelets. These
            % are adapted to the specific spectrum, not just the length of
            % the spectrum.
            case 'log_plus'
                % Additional required inputs: 
                % 1) param.warp_function

                if ~isfield(param,'warp_function')
                    error('Custom warp function is not specified')
                else
                    warp_function=@(s) log(param.warp_function(s)+eps);
                end
                filters=cell(num_filters,1);
                upper_bound_translates=warp_function(lmax);
                
                % Generate (num_filters-1) uniform translates covering [0,log(lmax)] 
                uniform_filters = gsp_uniform_translates_filter_design(num_filters-1,upper_bound_translates);
                
                % Warp the uniform translates by log to generate the "wavelet" filters
                wavelet_filters=cell(num_filters-1,1);
                for i=1:num_filters-1
                    wavelet_filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
                end
                filters(2:num_filters)=wavelet_filters;
                
                % Generate a "scaling" filter that results in a tight frame
                filters{1}= @(x) gsp_log_wavelet_scaling_fn(wavelet_filters,x);
            
            %%%%%%%%%%%%%%%%%
            % Warping functions based on approximations of the spectrum
            % (i.e., like the filter banks in Section 2, these are
            % spectrum-adapted filter banks).
            case 'pwl'
                % Additional required inputs:
                % 1) param.approx_spectrum
                if ~isfield(param,'approx_spectrum')
                    error('Approximate spectrum is not specified')
                end
                if ~isfield(param,'warp_start_one') % DIS: Get rid of two cases - just include in interpolation points
                    warp_start_one=0;
                else
                    warp_start_one=param.warp_start_one;
                end
                
                % Generate uniform translates covering [0,upper_bound_translates] % DIS: Get rid of two cases - just include in interpolation points
                if warp_start_one
                    warp_function = @(s) gsp_pwl_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s);
                    upper_bound_translates=max(param.approx_spectrum.y);
                else
                    warp_function = @(s) gsp_pwl_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s)-1;
                    upper_bound_translates=max(param.approx_spectrum.y)-1;
                end
                uniform_filters = gsp_uniform_translates_filter_design(num_filters,upper_bound_translates);
                
                % Generate warped filters 
                filters=cell(num_filters,1);
                for i=1:num_filters
                    filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
                end
                
            case 'mono_cubic'
                % Additional required inputs:
                % 1) param.approx_spectrum
                if ~isfield(param,'approx_spectrum')
                    error('Approximate spectrum is not specified')
                end
                
                % Generate uniform translates covering [0,upper_bound_translates]
                warp_function = @(s) gsp_mono_cubic_warp_fn(param.approx_spectrum.x,param.approx_spectrum.y,s);
                upper_bound_translates=max(param.approx_spectrum.y);
                uniform_filters = gsp_uniform_translates_filter_design(num_filters,upper_bound_translates);
                
                % Generate warped filters 
                filters=cell(num_filters,1);
                for i=1:num_filters
                    filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
                end
            
            %%%%%%%%%%%%%%%%%
            % Other user-specified warping function
            case 'custom'
                % Additional required inputs: 
                % 1) param.warp_function and
                % 2) param.upper_bound_translates (optional)

                if ~isfield(param,'warp_function')
                    error('Custom warp function is not specified')
                else
                    warp_function=param.warp_function;
                end

                if ~isfield(param,'upper_bound_translates')
                    param.upper_bound_translates=warp_function(lmax);
                    %error('Upper bound of the uniform translates is not specified')
                end
                filters=cell(num_filters,1);
                
                % Generate uniform translates covering [0,param.upper_bound_translates] 
                uniform_filters = gsp_uniform_translates_filter_design(num_filters,param.upper_bound_translates);
                
                % Generate custom warped filters (may be adapated to the specific spectrum via the custom warp function)
                for i=1:num_filters
                    filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
                end
                
            otherwise
                error('Warp function type not recognized')
        end
        varargout{1}=warp_function;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Section 5: Filter banks that are adapted to a class of graphs (via 
% warping), even though an approximation of the actual realization of the 
% spectrum is not known    

    %%%%%%%%%%%%%%%%%  
    % Combinatorial Laplacian Eigenvalue Warping Function on a Regular
    % Random Graph
    case 'random_regular_lap'
        % Assume N is large enough
        if ~isfield(param,'r')
            error('The degree is not specified')
        else
            r=param.r;
        end
        warp_function = @(z) real(1/2 + r/(2*pi)*asin((z-r)/(2*sqrt(r-1)))-(r-2)/(2*pi)*atan((r-2)*(z-r)./(r*sqrt(4*(r-1)-(z-r).^2)))).*(z >= r-2*sqrt(r-1)).*(z<r+2*sqrt(r-1))+(z>=r+2*sqrt(r-1));
        
        filters=cell(num_filters,1);

        % Generate uniform translates covering [0,param.upper_bound_translates] 
        uniform_filters = gsp_uniform_translates_filter_design(num_filters,1);

        % Generate warped filters 
        for i=1:num_filters
            filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
        end
        varargout{1}=warp_function;
        
    %%%%%%%%%%%%%%%%%  
    % Combinatorial Laplacian Eigenvalue Warping Function on an Erdos-Renyi
    % Random Graph
    case 'erdos_renyi_lap'
        if ~isfield(param,'N')
            error('The number of vertices N is not specified')
        else
            N=param.N;
        end
        if ~isfield(param,'p')
            error('The probability p of an edge is not specified')
        else
            p=param.p;
        end
        warp_function = @(s) gsp_erdos_renyi_warp(s,N,p);
        
        filters=cell(num_filters,1);

        % Generate uniform translates covering [0,param.upper_bound_translates] 
        uniform_filters = gsp_uniform_translates_filter_design(num_filters,1);

        % Generate warped filters 
        for i=1:num_filters
            filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
        end
        varargout{1}=warp_function;
        
    %%%%%%%%%%%%%%%%%  
    % Normalized Laplacian Eigenvalue Warping Function on an Erdos-Renyi
    % Random Graph
    case 'erdos_renyi_norm_lap'
        if ~isfield(param,'N')
            error('The number of vertices N is not specified')
        else
            N=param.N;
        end
        if ~isfield(param,'p')
            error('The probability p of an edge is not specified')
        else
            p=param.p;
        end
        warp_function = @(s) sqrt((p*N)/(1-p))*(1/(2*pi))*(s>1-2*sqrt((1-p)/(p*N))).*(s<1+2*sqrt((1-p)/(p*N))).*(pi/sqrt(p*N/(1-p))+((s-1)/2).*(sqrt(4-(p*N/(1-p))*(1-s).^2))-2/sqrt(p*N/(1-p))*asin(sqrt(p*N/(1-p))*(1-s)/2))+(s>=1+2*sqrt((1-p)/(p*N))); 
        
        filters=cell(num_filters,1);

        % Generate uniform translates covering [0,param.upper_bound_translates] 
        uniform_filters = gsp_uniform_translates_filter_design(num_filters,1);

        % Generate warped filters 
        for i=1:num_filters
            filters{i} = @(x) uniform_filters{i}(warp_function(x)); 
        end
        varargout{1}=warp_function;
        
    otherwise
        error('Unknown filter type');
end

end



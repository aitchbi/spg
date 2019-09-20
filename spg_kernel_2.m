%SPG_KERNEL_2 Second, third, ..., last spectral kernels of the uniform 
% Meyer-type system of spectral kernels,
% as presented in (cf. Proposition 1):
%
% Behjat, et al., 'Signal-adapted tight frames on graphs', 
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
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

function [r] = spg_kernel_2(x,iScale,lmax,N_scales,varargin)

control_params={'warping','none','E',[],'subSampleWarping',0}; 
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

t = 1;
M = 2.73; 
q= 1/(M-1);
a = lmax/(M*N_scales-N_scales+2);

r=zeros(size(x));

Delta = M*a-a;
scaleShift = Delta*(iScale-1);

% Apply the Warping to x (i.e. x = w(x))
if strcmp(warping,'none')
    
else
    theEigs = E;
    diffOfEigs = theEigs(2:end)-theEigs(1:end-1);
    indNonEqual = [1; find(logical(diffOfEigs)) + 1];
    
    if subSampleWarping
        dummy1 = E(indNonEqual);
        interp_x = dummy1(1:subSampleWarping:end-1);
        interp_x(end+1) = dummy1(end);
        
        dummy2 = warping(indNonEqual)*lmax;
        interp_y = dummy2(1:subSampleWarping:end-1);
        interp_y(end+1) = dummy2(end);
        
        x = x(x <= interp_x(end));
        x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,x)); 
    else
        interp_x = E(indNonEqual);
        interp_y = warping(indNonEqual)*lmax;
        x = abs(gsp_mono_cubic_warp_fn(interp_x,interp_y,x)); 
    end
end

r1 = find(x>= a+Delta*(iScale-1) & x<= M*a + Delta*(iScale-1));
r(r1) = sin(pi/2*(meyeraux(q*(t*(x(r1)- Delta*(iScale-1))/a -1))));

if iScale<N_scales
    r2 = find(x> M*a+Delta*(iScale-1) & x <= M*a + Delta*iScale);
    r(r2) = cos(pi/2*(meyeraux(q*(t*(x(r2)- Delta*iScale)/a -1))));
else
    r2 = find(x> M*a+Delta*(iScale-1));
    dummy = ones(size(x));
    r(r2) = dummy(r2);
end
end

function y = meyeraux(x)
p = [-20 70 -84 35 0 0 0 0];
y = polyval(p,x);
end

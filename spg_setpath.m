% SPG_SETPATH Set paths for SPG package. 
%
% NOTE:
% The SGWT toolbox (v. June 4, 2012) and some functions from 
% GSPbox and other resources are included in the SPG package; 
% Check spg/external/ to see the included files. 
% 
% If you already have these external files you may most likely be
% able to safely remove them from this package, but keep the 
% following file:
%
% spg/external/GSPbox/gsp_plot_graph.m
%
% as it is slightly modified from its original release. You may get 
% an error when running demo_minnesota, demo_alameda and 
% demo_cerebellum if another version of this fucntion is called.
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

function spg_setpath()

addpath(genpath(fileparts(mfilename('fullpath'))))

clc
fprintf('----------------------------------\n')
fprintf('SPG package (v1.00 - July 2016)\n')
fprintf('----------------------------------\n')

message = [...
    sprintf('NOTE:\n\n'),...
    sprintf('The SGWT toolbox (v. June 4, 2012) and some functions from'),...
    sprintf(' GSPbox and other resources are included in the SPG package;\n\n'),...
    sprintf('check spg/external/ to see the included files.\n\n'),...
    sprintf('If you already have these external files you may most likely be'),...
    sprintf(' able to safely remove them from this package, but keep the'),...
    sprintf(' following file: \n\n'),...
    sprintf('spg/external/GSPbox/gsp_plot_graph.m \n\n'),...
    sprintf('as this function is slightly modified from its original release.'),...
    sprintf(' You may get an error when running demo_minnesota,')...
    sprintf(' demo_alameda and demo_cerebellum if another'),...
    sprintf(' version of this fucntion is called.')];
%uiwait(msgbox(message));

fprintf([message,sprintf('\n')])
fprintf('----------------------------------\n')

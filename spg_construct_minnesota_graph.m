% SPG_CONSTRUCT_MINNESOTA_GRAPH This function is a minutely modified 
% & renamed version of the following function from GSPbox:
% 
% GSP_CREATE_PATH_GRAPH  Initialize the Minnesota road network
%   Usage:  G=gsp_create_minnesota_graph();
%
%   Input parameters:
%         none
%   Output parameters:
%         G     : Graph structure.
%
%   'gsp_create_minnesota_graph()' initializes a graph structure containing
%   the weighted adjacency matrix (G.W), the number of vertices (G.N), the 
%   plotting coordinates (G.coords), and the plotting coordinate limits 
%   (G.coord_limits) of the Minnesota road network from the MatlabBGL library.
%   We adjust the adjacency matrix so that all edge weights are equal to 1, 
%   and the graph is connected. 
%
%   References:
%     D. Gleich. The MatlabBGL Matlab library.
%     http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/index.html.
%     
%   Url: http://gspbox.sourceforge.net/doc/graphs/specific_graphs/gsp_create_minnesota_graph.php

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

function [G]=spg_construct_minnesota_graph()

Q=load('minnesota.mat');
G.N=size(Q.A,1);
G.coords=Q.xy;
G.coord_limits=[-98,43;-89,50];

% Edit adjacency matrix
A=Q.A;
% clean minnesota graph
A=A-diag(diag(A));
% missing edge needed to connect graph
A(349,355)=1;
A(355,349)=1;
% change a handful of 2 values back to 1
A(86,88)=1;
A(88,86)=1;
A(345,346)=1;
A(346,345)=1;
A(1707,1709)=1;
A(1709,1707)=1;
A(2289,2290)=1;
A(2290,2289)=1;
G.W=sparse(A);

%  Added (H. Behjat) ---------------------------
thresh=1e-8;
G.A=abs(G.W)>thresh; % unweighted adjacency matrix
G.laplacianType = 'normalized';
G.edge_width = 1;
G.edge_color = [0.5 0.5 0.5];
G.edge_style = '-';
G.vertex_size = 15;
G.vertex_color = 'k';
G.vertex_edge_color = 'k';

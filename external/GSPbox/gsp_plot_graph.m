function gsp_plot_graph(G,varargin)
%GSP_PLOT_GRAPH  Plot a graph in 2D or 3D
%   Usage:  gsp_plot_graph(G);
%           gsp_plot_graph(G,param);
%
%   Input parameters:
%         G          : Graph structure or a cell of graph structures.
%         param      : Optional variable containing additional parameters.
%   Output parameters:
%         none
%
%   Additional parameters:
%         param.show_edges     : Set to 0 to only draw the vertices.
%         param.num_clusters   : Number of clusters for a clustered graph.
%         param.clusters       : Cluster identities for a clustered graph.
%         param.cluster_colors : Cluster colors for a clustered graph.
%         param.cp             : Camera position for a 3D graph
%
%   'gsp_plot_graph(G)' plots a graph (or multiple graphs) in 2D or 3D, 
%   using the adjacency matrix (G.A), the plotting coordinates (G.coords), 
%   the coordinate limits (G.coord_limits), the edge width (G.edge_width), 
%   the edge color (G.edge_color), the edge style (G.edge_style), the 
%   vertex size (G_vertex_size), and the vertex color (G_vertex_color).
%
%   See also:  
%
%   Demos:  
%
%   Requires: gplot23D
%
%   Url: http://gspbox.sourceforge.net/doc/plotting/gsp_plot_graph.php

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
  
% Read input parameters
if nargin>1
    param=varargin{1};
else
    param=0;
end

if iscell(G)
    for i=1:length(G)
        gsp_plot_graph(G{i},param);
    end
else
    if ~isfield(G,'coord_limits')
        G.coord_limits=[min(G.coords);max(G.coords)];
    end

    if ~isfield(param,'show_edges')
        if nnz(G.A) > 85000 % HB changed
            show_edges=0;
        else
            show_edges=1;
        end
    else 
        show_edges=param.show_edges;
    end

    if (~isfield(param,'num_clusters') || ~isfield(param,'clusters') || ~isfield(param,'cluster_colors') )
        num_clusters=1;
        clusters=ones(G.N,1);
        cluster_colors=G.vertex_color;
    elseif length(param.cluster_colors)<param.num_clusters 
        error('Not enough cluster colors specified');
    else
        num_clusters=param.num_clusters;
        clusters=param.clusters;
        cluster_colors=param.cluster_colors;
    end
    
    figure;
    hold on;
    
    if show_edges
        gplot23D(G.A,G.coords,G.edge_style,'LineWidth',G.edge_width,'Color',G.edge_color);
    end

    for i=1:num_clusters

        if size(G.coords,2)==2
            scatter(G.coords((clusters==i),1),G.coords((clusters==i),2),G.vertex_size,'MarkerFaceColor',cluster_colors(i),'MarkerEdgeColor',G.vertex_edge_color);
            %xlim([G.coord_limits(1,1) G.coord_limits(2,1)]); % HB commented
            %ylim([G.coord_limits(1,2) G.coord_limits(2,2)]); % HB commented
        elseif size(G.coords,2)==3
            scatter3(G.coords((clusters==i),1),G.coords((clusters==i),2),G.coords((clusters==i),3),G.vertex_size,'MarkerFaceColor',cluster_colors(i),'MarkerEdgeColor',G.vertex_edge_color);
            %xlim([G.coord_limits(1,1) G.coord_limits(2,1)]); % HB commented
            %ylim([G.coord_limits(1,2) G.coord_limits(2,2)]); % HB commented
            %zlim([G.coord_limits(1,3) G.coord_limits(2,3)]); % HB commented
        end
        
    end
    
    % HB commented
    %     if size(G.coords,2)==3
    %         if ~isfield(param,'cp')
    %             cp=[-1.4,-16.9,3.4];
    %         else
    %             cp=param.cp;
    %         end
    %         set(gca,'CameraPosition',cp);
    %     end
    
    axis off;
end

end

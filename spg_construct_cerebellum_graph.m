% SPG_CONSTRUCT_CEREBELLUM_GRAPH Constructs the cerebellum graph 
% as presented in:
%
% Behjat, et al., 'Signal-adapted tight frames on graphs', 
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
%
% Output:
% 'G'       : Graph structure.
%
%
% Example:
%
% >> G = SPG_CONSTRUCT_CEREBELLUM_GRAPH;
% >> gsp_plot_graph(G)
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

function G = spg_construct_cerebellum_graph

dummy = load('cerebellum.mat');
GM = dummy.cerebellum;
mask = zeros(size(GM));
mask(GM>0)=1;
G.A = compute_adjacency(mask,26);
G.W = G.A;
G.N = size(G.A,1);

G.laplacianType = 'normalized';

G.edge_width = 1;
G.edge_color = [0.5 0.5 0.5];
G.edge_style = '-';
G.vertex_size = 15;
G.vertex_color = 'k';
G.vertex_edge_color = 'k';

[x,y,z] = ind2sub(size(mask),find(mask(:)));
G.coords = [x,y,z];
end


function A = compute_adjacency(mask,conn,varargin)
% A minimized version of the function
%
% gwspm_compute_adjacency.m
%
% from the toolbox:
%
% gwSPM (v1.00)
% Hamid Behjat
% http://miplab.epfl.ch/index.php/software/gwspm

dim=size(mask);
N=numel(mask);

%the indices of the non-zero mask elements
indices=find(mask);
[alli, allj, allk]=ind2sub(dim,indices);

switch conn
    case 6
        nN = 3;
    case 10
        nN = 5;
    case 18
        nN = 9;
    case 26
        nN = 13;
end

%the coordinates of the center points
ci=repmat(alli,nN,1);
cj=repmat(allj,nN,1);
ck=repmat(allk,nN,1);

switch conn
    case 6
        ni=[alli  ; alli+1; alli  ];
        nj=[allj+1; allj  ; allj  ];
        nk=[allk  ; allk  ; allk+1];
    case 10
        ni=[alli  ; alli+1; alli+1; alli+1; alli  ];
        nj=[allj+1; allj  ; allj+1; allj-1; allj  ];
        nk=[allk  ; allk  ; allk  ; allk  ; allk+1];
    case 18
        ni=[alli  ;alli+1;alli+1;alli+1;alli  ;alli  ;alli  ;alli+1;alli+1];
        nj=[allj+1;allj  ;allj-1;allj+1;allj  ;allj+1;allj+1;allj  ;allj  ];
        nk=[allk  ;allk  ;allk  ;allk  ;allk+1;allk-1;allk+1;allk-1;allk+1];
    case 26
        ni=[alli;alli+1;alli+1;alli+1;alli;alli;alli;alli+1;alli+1;alli+1;alli+1;alli+1;alli+1];
        nj=[allj+1;allj;allj-1;allj+1;allj;allj+1;allj+1;allj;allj;allj-1;allj+1;allj+1;allj-1];
        nk=[allk;allk;allk;allk;allk+1;allk-1;allk+1;allk-1;allk+1;allk-1;allk-1;allk+1;allk+1];
end

maskZ=cat(2,mask,zeros(dim(1),1,dim(3)));
maskZ=cat(1,maskZ,zeros(1,dim(2)+1,dim(3)));
maskZ=cat(3,maskZ,zeros(dim(1)+1,dim(2)+1,1));

valid=(ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2) & nk>=1 & nk<=dim(3));
ni=ni(valid);
nj=nj(valid);
nk=nk(valid);
ci=ci(valid);
cj=cj(valid);
ck=ck(valid);

tt=sub2ind(size(maskZ),ni,nj,nk);
ee=maskZ(tt);
valid=logical(ee);
ni=ni(valid);
nj=nj(valid);
nk=nk(valid);
ci=ci(valid);
cj=cj(valid);
ck=ck(valid);

cInd=sub2ind(dim,ci,cj,ck);
nInd=sub2ind(dim,ni,nj,nk);

A=sparse([cInd,nInd],[nInd,cInd],ones(1,2*numel(ni)),N,N);

c=find(~sum(A,1));
c=c(~ismember(c,indices));

A(:,c)=[];
A(c,:)=[];

end
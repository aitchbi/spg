% SPG_DEMO_UNIFORMMEYERTYPE Constructs uniform Meyer-type system 
% of spectral kernels as presented in:
%
% Behjat, et al., 'Signal-adapted tight frames on graphs', 
% IEEE Trans. Signal Process., 2016,  doi: 10.1109/TSP.2016.2591513.
%
% Examples:
%
% 1) Constrcut uniform Meyer-type system of kernles, with 5 subbands, 
%     over spectral range: [0,2]:
%
% >> spg_demo_uniformMeyerType('N_subbands',[5],'lambdaMax',2) 
%
% 2) Constrcut uniform Meyer-type system of kernles, with 5,7 & 9 
%     subbands, over spectral range: [0,1]:
%
% >> spg_demo_uniformMeyerType('N_subbands',[5,7,9],'lambdaMax',1) 
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

function spg_demo_uniformMeyerType(varargin)

control_params={'N_subbands',[5,7,9],'lambdaMax',2};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);
fprintf('-----------------------------------\n')

hFig = figure('Name','Uniform Meyer Type system of spectral kernels','NumberTitle','off');
set(hFig, 'Position', [715, 10, 700, 150*numel(N_subbands)]);

% Meyer-type Uniform system of spectral kernels
for i=1:numel(N_subbands)
    g5 = spg_filter_design(lambdaMax,N_subbands(i),...
        'designtype','uniform_meyer_type');
    subplot(numel(N_subbands),1,i)
    spg_view_design(g5,[0,lambdaMax],'plotLineWidth',2);
    title(sprintf('Uniform Meyer-type Frame -- Number of subbands: %d -- lMax: %d',N_subbands(i),lambdaMax))
    set(gca,'Box','off','XLim',[0 lambdaMax],'YLim',[0 1.2])
end

if ~nargin
    fprintf('--------------------------------------------\n')
    message = [...
        sprintf('Uniform Meyer-type (UMT) sytems of spectral kernels with 5, ')...
        sprintf(' 7 and 9 subbabnds, all with a maximal spectral range of 2'),...
        sprintf(' are displayed.\n\n'),...
        sprintf('You may specify a different number of subbands and a different'),...
        sprintf(' maximal spectral range by changing the ''N_subbands'' and'),...
        sprintf(' ''lambdaMax'' parameter when running the demo; for instance:\n\n'),...
        sprintf('>> spg_demo_uniformMeyerType(''N_subbands'',6,''lambdaMax'',1) \n\n'),...
        ];
    uiwait(msgbox(message));
end

end
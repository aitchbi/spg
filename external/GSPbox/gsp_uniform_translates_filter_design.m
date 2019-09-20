function [ filters ] = gsp_uniform_translates_filter_design( num_filters, upper_bound_translates )
% Helper function to generate uniform translate filter banks. Need to have
% it as a separate function to use it in the warped filter bank generation
%
%   Url: http://gspbox.sourceforge.net/doc/filter_design/utils/gsp_uniform_translates_filter_design.php

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

% Design main window
dilation_factor = upper_bound_translates*(3/(num_filters-2));
main_window=@(x) (.5+.5*cos(2*pi*(x/dilation_factor-1/2))).*(x>=0).*(x<=dilation_factor);

% Design filters (uniform translates of the main window, cut-off at the spectrum boundaries)
filters=cell(num_filters,1);
for i=1:num_filters
   filters{i} = @(x) main_window(x-dilation_factor/3*(i-3)); 
end

end



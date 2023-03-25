% This file is part of Biot1D-MATLAB
% Biot1D-MATLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Biot1D-MATLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Biot1D-MATLAB. If not, see <https://www.gnu.org/licenses/>.
% The full text of the license can be found in the file License.md.
%%
global MYCASEFLAG
if MYCASEFLAG ==0
    a = 0; b = 1;
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa  = 1;
elseif MYCASEFLAG ==1 %% Example [NV/1]
    a = 0; b = 1;
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa = 1;
elseif MYCASEFLAG ==2 %% Example [NV/2]
    a = 0; b = 1;
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa = 1;
elseif MYCASEFLAG ==3 %% Example [NV/3]
    a = 0; b = 1;
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa  = 1;
elseif MYCASEFLAG ==4 %% Example [NV/4]
    a = 0; b = 0.1;
    E = 20;
    nu = 0.30;
    COF_c0 = 2.08e-4; 
    COF_lambda = (E*nu)/((1 + nu)*(1-2*nu)); COF_mu = E/(2*(1+nu));
    COF_alpha = 1;
    COF_kappa  = 1e-17/2.7822e-13;
elseif MYCASEFLAG ==5 %% Example [NV/5]
    a = 0; b = 1;
    COF_c0 = 1.5750e-10; 
    COF_lambda = 2.8846e7; 
    COF_mu = 1.9231e7; 
    COF_alpha = 1;
    COF_kappa  = 8e-16/1.1390e-6;
else
    error('Example/scenario not implemented');
end
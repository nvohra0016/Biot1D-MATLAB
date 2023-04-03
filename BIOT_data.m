% This file is part of Biot1D-MATLAB
% Biot1D-MATLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Biot1D-MATLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Biot1D-MATLAB. If not, see <https://www.gnu.org/licenses/>.
% The full text of the license can be found in the file License.md.
%%
global MYCASEFLAG
if MYCASEFLAG ==0 %% custom example
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa  = 1;
elseif MYCASEFLAG ==1 %% Example [NV/1]; manufactured solution
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa = 1;
elseif MYCASEFLAG ==2 %% Example [NV/2]; manufactured solution
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa = 1;
elseif MYCASEFLAG ==3 %% Example [NV/3]; manufactured solution 
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa  = 1;
elseif MYCASEFLAG ==4 %% Example [NV/4]; clay consolidation
    E = 20;
    nu = 0.30;
    COF_c0 = 2.08e-4; 
    COF_lambda = 11.538461538461538;
    COF_mu = 7.692307692307692;
    COF_alpha = 1;
    COF_kappa  = 1e-17/2.7822e-13;
elseif MYCASEFLAG ==5 %% Example [NV/5]; clay-sand consolidation, heterogeneous example
    COF_alpha = 1;
    COF_rhof = 998.21*(1e-6 * (1/3600)*(1/3600)); 
    COF_G =  1.27290528e8 * 1;
    if exist('x','var')
        COF_c0 = 0*x; COF_c0(x <= 0.5) = 0.3 * 4.16e-4; COF_c0(x > 0.5) = 0.5 * 4.16e-4;
        COF_lambda = 0*x; COF_lambda(x <= 0.5) = 6; COF_lambda(x > 0.5) = 11.538461538461538;
        COF_mu = 0*x; COF_mu(x <= 0.5) = 6; COF_mu(x > 0.5) = 7.692307692307692;
        COF_kappa = 0*x; COF_kappa(x <= 0.5)  = 1e-12/2.7822e-13; COF_kappa(x > 0.5) = 1e-17/2.7822e-13;
        COF_rhos = 0*x; COF_rhos(x <= 0.5) = 2650*(1/3600)*(1/3600)*1e-6; COF_rhos(x > 0.5) = 2700*(1/3600)*(1/3600)*1e-6;
    end
else
    error('Example/scenario not implemented');
end
% This file is part of Biot1D-MATLAB
% Biot1D-MATLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Biot1D-MATLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Biot1D-MATLAB. If not, see <https://www.gnu.org/licenses/>.
% The full text of the license can be found in the file License.md.
%% Data used in poroelasticity code.
%% Units: [m, s, Pa]
global MYCASEFLAG
if MYCASEFLAG ==0 %% custom example
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa  = 1;
elseif MYCASEFLAG ==1 %% Example 1; manufactured solution
    COF_c0 = 1; 
    COF_lambda = 1; COF_mu = 1; 
    COF_alpha = 1;
    COF_kappa = 1;
    COF_rhor = 1;
    COF_rhol = 1;
    COF_G = 0;
elseif MYCASEFLAG ==2 %% Example 2; manufactured solution
    COF_c0 = 1;
    COF_lambda = 1; COF_mu = 1;
    COF_alpha = 1;
    COF_kappa = 1;
    COF_rhor = 1;
    COF_rhol = 1;
    COF_G = 0;
elseif MYCASEFLAG ==3 %% Example 3; manufactured solution
    COF_c0 = 1;
    COF_lambda = 1; COF_mu = 1;
    COF_alpha = 1;
    COF_kappa  = 1;
    COF_rhor = 1;
    COF_rhol = 1;
    COF_G = 0;
elseif MYCASEFLAG ==4 %% Example 4; clay consolidation
    E = 20e6;
    nu = 0.30;
    COF_c0 = 0.5 * 4.58e-10; 
    COF_lambda = 1.1538e+07;
    COF_mu = 7.6923e+06;
    COF_alpha = 1;
    COF_viscosity = 1.0005e-3;
    COF_kappa  = 1e-17/COF_viscosity;
    COF_rhor = 2650;
    COF_rhol = 998.21; 
    COF_G = 9.8218*0;
elseif MYCASEFLAG ==5 %% Example 5; clay-sand consolidation, heterogeneous example
    COF_alpha = 1;
    COF_rhol = 998.21; 
    COF_G =  9.8218;
    COF_viscosity = 1.0005e-3;
    if exist('x','var')
        COF_c0 = 0*x; COF_c0(x <= 0.5) = 0.3 * 4.58e-4; COF_c0(x > 0.5) = 0.5 * 4.58e-4;
        COF_lambda = 0*x; COF_lambda(x <= 0.5) = 6e6; COF_lambda(x > 0.5) = 1.1538e+07;
        COF_mu = 0*x; COF_mu(x <= 0.5) = 6e6; COF_mu(x > 0.5) = 7.6923e+06;
        COF_kappa = 0*x; COF_kappa(x <= 0.5)  = 1e-12/COF_viscosity; COF_kappa(x > 0.5) = 1e-17/COF_viscosity;
        COF_rhor = 0*x; COF_rhor(x <= 0.5) = 2650; COF_rhor(x > 0.5) = 2700;
    end
else
    error('Example/scenario not implemented');
end
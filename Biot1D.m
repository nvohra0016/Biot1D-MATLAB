% This file is part of Biot1D-MATLAB
% Biot1D-MATLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Biot1D-MATLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Biot1D-MATLAB. If not, see <https://www.gnu.org/licenses/>.
% The full text of the license can be found in the file License.md.
%%
function   [xfem,usol,xplot,psol]=...
    Biot1D (Tend,nx,dt,bdaryflags,caseflag,ifsave,ifplot) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% <bdaryflags> vector of flags for [MMHH], Dirichlet (0) or Neumann (1) flags 
% (Note: user must code exfun, dexfun to deliver these values)
%% <caseflag> case study number, drawn from BIOT_data, and hard-coded examples 
% only caseflag ==4 has NO exact soln known
%% <ifsave>: flag for saving last time step in a file [->0]
%ifsave =0: no saving; 
%ifsave>0: save to file at the end; 
%ifsave=-1: debug
%ifsave=-2: compute error
%% <ifplot> flag for plotting [->0]
%ifplot>0, then plot at time steps=multiples of ifplot, and at n=1, n=end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLOBALS and files:
% MYCASEFLAG, and file BIOT_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXAMPLES
% case 1, compute error, 
%[xu,u,xp,p]=Biot1d(0.1,10,0.1,[0,0,0,0],1,-2);
%[xu,u,xp,p]=Biot1d(0.1,10,0.1,[1,0,1,1],1,-2);
% case 2, plot, save to file
%[xu,u,xp,p]=Biot1d(0.1,10,0.1,[1,0,1,1],2,1,1);
% case 3, no file, plot every 5
%[xu,u,xp,p]=Biot1d(1,10,0.1,[1,0,1,1],3,0,5);
%% DEFAULT PLOT PROPERTIES
set(groot,'defaultLineLineWidth',4)
set(0,'DefaultaxesLineWidth', 3)
set(0,'DefaultaxesFontSize', 16)
set(0,'DefaultTextFontSize', 16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA
global MYCASEFLAG
MYCASEFLAG = caseflag;
%
if (MYCASEFLAG == 1) a=0; b=1; Tend=1; nx=10; dt=0.1; bdaryflags=[0,0,0,0];
elseif (MYCASEFLAG == 2) a=0; b=1; Tend=1; nx=[0;0.05;0.1;0.15;0.2;0.4;0.6;0.8;0.9;0.95;1.0]; dt=0.1; bdaryflags=[1,0,1,0];
elseif (MYCASEFLAG == 3) a=0; b=1; Tend=1; nx=20; dt=0.1; bdaryflags=[0,1,1,1];
elseif (MYCASEFLAG == 4) a=0; b=0.1; Tend=24; nx=20; dt=2.4; bdaryflags=[1,0,0,1];
elseif (MYCASEFLAG == 5) a=0; b=1; Tend=8760; nx=20; dt=87.6; bdaryflags=[1,0,0,1];
end

%if nargin<5,Tend=1;nx=10;dt=0.1;bdaryflags=[0,0,0,0];caseflag=1;end
%if nargin<6, ifsave=0;end
%if nargin<7, ifplot=0;end
%
if caseflag<4, ifexact=1; else, ifexact=0; end
%
BIOT_data;
%
%%%%%%%%%%%%% grid and data structures 
%%
t1 = 0; t2 = Tend; nt = (t2-t1) / dt; 
%% GRID
% dx is grid size of each CCFD cell;
% x is position of left node of each cell
% xplot is the position of CCFD values
% xfem is position for displacements: x expanded by last node
%% UNIFORM GRID
if (size(nx,1) == 1)
    dx = (b-a)/nx * ones(nx,1); x0 = a;
    x = 0*dx; x(1)=x0; for j=2:nx, x(j)=x(j-1)+dx(j-1); end
    xplot = x + dx/2;
    xfem = x; xfem(nx+1)=x(nx)+dx(nx);
%% NON-UNIFORM GRID    
else
    x = nx;
    nx = size(nx,1) - 1;
    dx = x(2:end) - x(1:end-1);
    xplot = x(1:end-1) + dx/2;
    xfem = x;
 end
%%%% FLOW part
%% permeability coefficient, porosity coefficient, average media density per cell
perm = permfun(xplot);
por = porfun(xplot);
rhoavg = por.*rhosfun(xplot) + (1-por).*COF_rhof;
%% set-up matrices
% compute transmissibilities
tx = zeros(nx+1,1); 
for j=2:nx, tx(j)=2/(dx(j-1)/perm(j-1)+dx(j)/perm(j)); end
% Dirichlet 
if bdaryflags (3) == 0, j = 1;      tx(j)=2/(dx(j)/perm(j)); end
if bdaryflags (4) == 0, j = nx + 1; tx(j)=2/(dx(j-1)/perm(j-1)); end
% for time dependent case
tx = tx*dt;
%%    
stiff = sparse (nx,nx);
for j=2:nx                                       
  gl = j-1; gr = j;
  stiff(gl,gl) = stiff(gl,gl) + tx(j);                                    
  stiff(gl,gr) = stiff(gl,gr) - tx(j);                                                              
  stiff(gr,gl) = stiff(gr,gl) - tx(j);                                    
  stiff(gr,gr) = stiff(gr,gr) + tx(j);     
end
%% Dirichlet contributions to matrix
if bdaryflags(3) == 0, j = 1; gr = 1; stiff(gr,gr) = stiff(gr,gr) + tx(j); end
if bdaryflags(4) == 0, j = nx + 1; gl = nx; stiff(gl,gl) = stiff(gl,gl) + tx(j); end
%%
%full(stiff), pause
%
mass = diag(por.*dx);

%%%% ELASTICITY: 
%% elasticity coefficient
elcof = elcoffun(xplot);
elstiff = sparse(nx+1,nx+1);
elq = zeros(nx+1,1);
%
for j=2:nx %% arithmetic averaging for CG
    elstiff (j,j)  = elcof(j-1)/dx(j-1)+elcof(j)/dx(j);
    elstiff (j,j-1) =-elcof(j-1)/dx(j-1);
    elstiff(j,j+1)  =-elcof(j)/dx(j);
end
%% BCONDITIONS: elasticity, implemented in the entire system. 
if bdaryflags (1) == 0, j=1; multi = elcof(j)/dx(j); elstiff(j,j) = multi; 
else, j=1; multi = elcof(j)/dx(j); elstiff(j,j) = multi; elstiff(j,j+1) =-multi; 
end
%%
if bdaryflags (2) == 0,j=nx+1; multi = elcof(j-1)/dx(j-1); elstiff(j,j) = multi; 
else, j=nx+1; multi = elcof(j-1)/dx(j-1); elstiff(j,j) = multi; elstiff(j,j-1) = -multi; 
end

%% Coupling terms
stiff_pu = sparse(nx,nx+1); %% pressure equation: coupling to displacement
for j=1:nx,  stiff_pu(j,j+1) = COF_alpha; stiff_pu(j,j)=-COF_alpha; end
% 
stiff_up = -stiff_pu';
% remove the entries in the rows of M with Dirichlet b.c. for u
if bdaryflags(1)==0, stiff_up(1,1:end)=0; end
if bdaryflags(2)==0, stiff_up(nx+1,1:end)=0; end

%% Global matrix: assemble
nM = nx+1; nH = nx; nMind = 1:1:nM; nHind = nM+1:1:nM+nH;
mat = sparse(nM + nH, nM + nH);
% MM
mat(nMind,nMind) = elstiff(:,:);
% HH
%same as mat(nM+1:end,nM+1:end) = stiff(:,:) + mass(:,:);
mat (nHind,nHind) = stiff(:,:) + mass(:,:);
% coupling terms
mat(nHind,nMind) = stiff_pu;
mat(nMind,nHind) = stiff_up;
%
if ifsave==-1, full(mat), end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  TIME LOOP AND SOLVER
fcontent = fluid_init(xplot,caseflag);

%% 
%errvec = [];
t = t1;
for n = 1:nt %  time step loop
    t = t + dt;
    fcontent_old = fcontent;
    %% piecewise constant rule on every element for rhs in P
    q = dt * dx.*p_rhs (xplot, t, caseflag); 
    q = q + dx.*fcontent_old;
    
    %% trapezoidal rule on every element for rhs in U
    elq = 0*xfem;
    for j=2:nx, elq(j) = (dx(j-1)+dx(j))/2*u_rhs(xfem(j),t,caseflag); end
    j=1;    elq(j) = dx(j)/2*u_rhs(xfem(j),t,caseflag);
    j=nx+1; elq(j) = dx(j-1)/2*u_rhs(xfem(j),t,caseflag);
    % gravity term
    for j=2:nx, elq(j) = elq(j) + COF_G*(dx(j-1)*rhoavg(j-1)+dx(j)*rhoavg(j))/2; end
    j=1;    elq(j) = elq(j) + COF_G*dx(j)/2*rhoavg(j);
    j=nx+1; elq(j) = elq(j) + COF_G*dx(j-1)/2*rhoavg(j-1);
    %fprintf('rhs\n');elq',
    %% get the values of bconditions or fluxes
    if ifexact == 0 %% no analytical solution
        pval1 = 0; pval2 =0; uval1 = -0.1; uval2 = 0;
    else  %% known analytical solutions or another case
        if bdaryflags(1)==0, uval1 = u_exfun(xfem(1),t,caseflag); else, uval1 = elcof(1)*u_dexfun(xfem(1),t,caseflag)-COF_alpha*p_exfun(xfem(1),t,caseflag); end
        if bdaryflags(2)==0, uval2 = u_exfun(xfem(end),t,caseflag); else, uval2 = elcof(nx)*u_dexfun(xfem(end),t,caseflag)-COF_alpha*p_exfun(xfem(end),t,caseflag); end
        if bdaryflags(3)==0, pval1 = p_exfun(xfem(1),t,caseflag); else, pval1 = p_dexfun(xfem(1),t,caseflag); end
        if bdaryflags(4)==0, pval2 = p_exfun(xfem(end),t,caseflag); else, pval2 = p_dexfun(xfem(end),t,caseflag); end
    end
    % contributions to [M] ELASTICITY from Dirichlet and Neumann    
    if bdaryflags(1) == 0, j=1; multi = elcof(j)/dx(j); elq(j) = uval1*multi; 
    else, j=1; elq(j)=  elq(j) - uval1; end
    %
    if bdaryflags(2) == 0, j=nx+1; multi = elcof(j-1)/dx(j-1); elq(j) = uval2*multi;
    else, j=nx+1; elq(j) =  elq(j) + uval2; end
    %% trapezoidal rule for gravity term in P
    rhsh = 0*xfem;
    for j=2:nx, rhsh(j) = COF_rhof * COF_G * (dx(j-1)+dx(j))/2; end
    j = 1;    rhsh(j) = COF_rhof * COF_G * dx(j)/2;
    j = nx+1; rhsh(j) = COF_rhof * COF_G * dx(j-1)/2;
    rhsh = rhsh .* tx;
    rhshn = rhsh(2:end) - rhsh(1:end-1);
    q = q - rhshn;

    % contributions of [H] FLOW from Dirichlet bdary conditions 
    if bdaryflags(3) == 0, q(1) = q(1) + tx(1) * pval1; else, q(1) = q(1) - dt*pval1 - rhsh(1) ; end
    if bdaryflags(4) == 0,  q(nx) = q(nx) + tx(nx+1) * pval2; else, q(nx) = q(nx) + dt*pval2 + rhsh(nx+1); end
    %debug
    if ifsave==-1 && ifexact ~=0 %% check with exact solution
        fprintf('bcond\n');
       % q',((stiff+mass)*p_exfun(xplot,t,caseflag)+stiff_pu*u_exfun(xfem,t,caseflag))',
        elq', (elstiff*u_exfun(xfem,t,caseflag)+stiff_up*p_exfun(xplot,t,caseflag))',
    end
    
    %% Solve the COUPLED FLOW and ELASTICITY
    rhs = zeros((nx+1) + (nx),1);
    rhs(nx+2:end) = q; rhs(1:nx+1) = elq;
    %
    allsol = mat \ rhs;
    %
    usol = allsol(1:nx+1);psol = allsol(nx+2:end); 
    
    %% postprocess for next time step
    fdilate = (stiff_pu*usol)./dx;
    fcontent = por.*psol + fdilate;
    %% plot etc.
    if ifexact>0,
         pex = p_exfun(xplot,t,caseflag); uex=u_exfun(xfem,t,caseflag);
    end
    if ifplot>0 && (rem(n,ifplot)==0 || n==1 || t==Tend)
        if ifexact == 0
            fprintf('Time step %d at t = %g\n',n,t);
            subplot(2,1,1);plot(xplot,psol,'b--o',xfem,usol,'r--*','linewidth',3,'MarkerSize',12);
            subplot(2,1,2);plot(xplot,psol,'b--o',xfem,usol,'r--*','linewidth',3,'MarkerSize',12);
            if caseflag ==4,
                % displacement
                subplot(2,1,1);plot(xfem,usol,'r--*','linewidth',3,'MarkerSize',12);
                lh = legend('u');
                set(lh,'FontSize',16,'location','northeast');
                legend boxoff;
                xlabel('x [m]');
                ylabel('u [m]')
                xlim([a b]);
                box off;
                % pressure
                subplot(2,1,2);plot(xplot,psol,'b--o','linewidth',3,'MarkerSize',12);
                lh = legend('p');
                set(lh,'FontSize',16,'location','northwest');
                legend boxoff;
                xlabel('x [m]');
                ylabel('p [MPa]')
                xlim([a b]);
                box off;
                sgtitle([sprintf('t=%g M=%d tau=%g. CASE=%d. BDARY=[%g %g %g %g]',t,nx,dt,caseflag,bdaryflags)],'FontSize',15);
            elseif caseflag ==5
                % displacement
                subplot(2,1,1);plot(xfem,usol,'r--*','linewidth',3,'MarkerSize',12);
                lh = legend('u','FontSize',16);
                set(lh,'FontSize',16,'location','northwest');
                legend boxoff;
                xlabel('x [m]');
                ylabel('u [m]')
                xlim([a b]);
                box off;
                % pressure
                subplot(2,1,2);plot(xplot,psol,'b--o','linewidth',3,'MarkerSize',12);
                lh = legend('p');
                set(lh,'FontSize',16,'location','northwest');
                legend boxoff;
                xlabel('x [m]');
                ylabel('p [MPa]')
                xlim([a b]);
                box off;
                sgtitle(sprintf('t=%g M=%d tau=%g. CASE=%d. BDARY=[%g %g %g %g]',t,nx,dt,caseflag,bdaryflags),'FontSize',15);
            end
            pause(.01);
        else
            % displacement
            plot(xplot,pex,'b',xfem,uex,'r','linewidth',3,'MarkerSize',12);
            % pressure
            hold on; 
            plot(xplot,psol,'b--o',xfem,usol,'r--*','linewidth',3,'MarkerSize',12);
            hold off;   
            lh = legend('p','u','pex','uex');
            set(lh,'FontSize',16,'location','northwest');
            legend boxoff;
            xlabel('x');
            xlim([a b]);
            box off;
            title(sprintf('Solution at t=%g M=%d tau=%g. CASE=%d. BDARY=[%g %g %g %g]',t,nx,dt,caseflag,bdaryflags));
            pause(.01);
        end
    end
end
%% save to file at the end of simulation
if ifsave>0
    fname = sprintf('BIOT_step%d_case%d_bd%d%d%d%d_alpha%d_c0%d.txt',n,caseflag,...
        bdaryflags(1),bdaryflags(2),bdaryflags(3),bdaryflags(4),floor(COF_alpha),floor(COF_c0));
    uall = full(usol);pall=full(psol);
    save(fname,'xfem','uall','xplot','pall','-ascii','-double');
end
if ifsave == -2 % compute error
    hmax = max(dx);
    fprintf('[H] dx=%g dt=%g L2err = %g\n',hmax,dt,norm(psol-pex,2)*sqrt(dx(1)));
    %% for FE in [M]
    uhfun = @(x)(interp1(xfem,usol,x,'linear','extrap'));
    %uhfun(0.5),u_exfun(0.5,t,caseflag)
    uherrfun = @(xp,tt,cc)((u_exfun(xp,tt,cc)-uhfun(xp)).^2);
    ul2 = integral(@(xp)uherrfun(xp,t,caseflag),a,b);
    fprintf('[M] dx=%g dt=%g L2err = %g (L2errappr = %g) Linferr = %g\n',...
        hmax,dt,sqrt(ul2),norm(usol-uex,2)*sqrt(dx(1)),norm(usol-uex,inf));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = permfun(x)
BIOT_data;
v = 0*x + COF_kappa;
end

function v = porfun(x)
BIOT_data;
v = 0*x + COF_c0;
end

function v = rhosfun(x)
BIOT_data;
v = 0*x + COF_rhos;
end

function v = elcoffun(x)
BIOT_data;
v = 0*x + COF_lambda + 2*COF_mu;
end

%%%%%%%%%%
function u = u_exfun(x,t,mycase)
BIOT_data;
if mycase ==1
    u = -1/(COF_lambda*pi)*sin(pi*t/2).*cos(pi*x); 
elseif mycase ==2
    u = sin(pi*x/2)*exp(-t);
elseif mycase==3
    u = 2 - x;
else
    error('u_exfun not available');
end
end

function udx = u_dexfun (x,t,mycase)
BIOT_data;
if mycase ==1
    udx = 1/(COF_lambda)*sin(pi*t/2).*sin(pi*x); 
elseif mycase ==2
    udx = pi/2*cos(pi*x/2)*exp(-t);
elseif mycase ==3
    udx = 0*x-1;
else
    error('u_dexfun not implemented');
end
end

function p = p_exfun(x,t,mycase)
BIOT_data;
if mycase ==1 
    p = sin(pi*t/2)*sin(pi*x);
elseif mycase == 2
    p = cos(pi*x/2)*exp(-t);
elseif mycase ==3
    p = x+1;
else
    error("p_exfun not implemented");
end
end

function pdx = p_dexfun (x,t,mycase)
BIOT_data;
if mycase ==1
    pdx = pi * sin(pi*t/2) * cos(pi*x);
elseif mycase ==2
    pdx = -pi/2 * sin(pi*x/2) * exp(-t);
elseif mycase ==3
    pdx = 0*x + 1;
else
    error('p_dexfun not implemented');
end
end

function v = fluid_init(x,mycase)
BIOT_data;
if mycase <4
    v = COF_c0*p_exfun(x,0,mycase) + COF_alpha*u_dexfun(x,0,mycase);
else
    v = 0*x + COF_G*COF_rhof*COF_c0.*x;
end
end

function pfun = p_rhs(x,t,mycase)
BIOT_data;
%
if mycase ==1
    pfun = (COF_c0 + COF_alpha/COF_lambda)*pi/2*cos(pi*t/2)*sin(pi*x)+...
        COF_kappa*pi*pi*sin(pi*t/2).*sin(pi*x);
elseif mycase == 2
    pfun = (-COF_c0 - COF_alpha*pi/2 + COF_kappa*pi^2/4)*cos(pi*x/2)*exp(-t);
elseif mycase == 3
    pfun =0*x;
elseif mycase >=4
    pfun = 0*x;
end
end

function ufun = u_rhs(x,t,mycase)
%
BIOT_data;
if mycase == 1
    ufun = (-(COF_lambda+2*COF_mu)*pi/COF_lambda + COF_alpha*pi)*sin(pi*t/2).*cos(pi*x); 
elseif mycase == 2
    ufun = ((COF_lambda+2*COF_mu)*pi^2/4-COF_alpha*pi/2) * sin(pi*x/2)*exp(-t);
elseif mycase ==3
    ufun = 0*x + COF_alpha;
elseif mycase >=4
    ufun = 0*x;
end
end


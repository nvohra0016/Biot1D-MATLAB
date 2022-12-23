%%
% This file is part of Biot1D-MATLAB
% Biot1D-MATLAB is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
% Biot1D-MATLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with Biot1D-MATLAB. If not, see <https://www.gnu.org/licenses/>.
% The full text of the license can be found in the file License.md.
%%
%% Biot 3 field poroelasticity solver (see Documentation pdf for full implementation details)
%% Input:
%% el_nodes (scalar or vector): 
% if scalar then it is the number of cells and a uniform spatial grid is
% created. If it is a vector, then it contains the element nodes as a column vector (nodal grid), starting
% from x = a and including x = b, i.e., size(M) = #cells + 1. 
%% tsteps (scalar): 
% number of time steps; this determines the time step size tau.
%% BC_flags (4x1 vector):
% boundary condition flags 1x4 vector with values 0 (for Dirichlet) or 1 (for Neumann)
% boundary conditions [M M H H] 
% For example: [1 0 0 1] = [N D D N] for Terzaghi's problem.
% [0 0 0 1] = [D D D N] corresponds to Dirichlet for mechanical and mixed for flow.
%% example (integer):
% parameter that controls the default coded examples.
% example = 1, 2, or 3 corresponds to manufactured solutions
% example = 4 corresponds to Terzaghi's problem of soil consolidation
% using physical parameters.
% example = 5 corresponds to Terzaghi's problem with 2 different soil
% types.
% example >= 6 are custom examples need to be set up by user (i.e., provide
% spatial domain, time period, and physical parameters).
%%
%%
%%
function [xn, xcc, t, U, P, Q] = Biot1D(el_nodes, ntsteps, BC_flags, example)
tic
%% set default plot values
set(groot,'defaultLineLineWidth',4)
set(0,'DefaultaxesLineWidth', 3)
set(0,'DefaultaxesFontSize', 24)
set(0,'DefaultTextFontSize', 18)
%% set boundary conditions. Dirichlet by default, so Neum == 1 for Neumann condition.
Neum_ua = BC_flags(1);
Neum_ub = BC_flags(2);
Neum_pa = BC_flags(3);
Neum_pb = BC_flags(4);
%% spatial and temporal domain
[a, b, Tend] = grid(example);
%% spatial grid size/ nodal/cell centered grid construction
if (size(el_nodes,1) == 1)
    M = el_nodes;
    h = (b-a)/M * ones(M,1);
    %% nodal grid
    xn = a:h:b;
    xn = xn';
    %% cell centered grid
    xcc = a+(h/2):h:b;
    xcc = xcc';
else
    M = size(el_nodes,1)-1;
    h = zeros(M,1);
    for j = 1:1:M
        h(j) = el_nodes(j+1)-el_nodes(j);
    end
    %% nodal grid
    xn = el_nodes(:);
    %% cell centered grid
    xcc = zeros(M,1);
    for j = 1:1:M
        xcc(j) = (el_nodes(j) + el_nodes(j+1))/2;
    end
end
%% time step size and temporal grid
tau = (Tend-0)/ntsteps;
t = 0:tau:Tend;
t = t(:);
%% physical parameters
[lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(xcc, example);
%%
%% compute FE matrices
%%
%% displacement stiffness matrix
Auu = sparse(M+1,M+1);
for j = 2:1:M
    Auu(j,j-1) = -1/h(j-1)*(lambda(j-1) + 2*mu(j-1));
    Auu(j,j) = (1/h(j-1))*(lambda(j-1) + 2*mu(j-1)) + (1/h(j))*(lambda(j) + 2*mu(j));
    Auu(j,j+1) = -1/h(j)*(lambda(j) + 2*mu(j));
end
Auu(1,1) = (1/h(1))*(lambda(1) + 2*mu(1)); Auu(1,2) = -1/h(1)*((lambda(1) + 2*mu(1)));
Auu(M+1,M) = (-1/h(M))*(lambda(M) + 2*mu(M)); Auu(M+1,M+1) = (1/h(M))*(lambda(M) + 2*mu(M));
%% pressure-displacement stiffness matrix
Apu = sparse(M+1,M);
for j = 2:1:M
    Apu(j,j-1) = 1;
    Apu(j,j) = -1;
end
Apu(1,1) = -1;
Apu(M+1,M) = 1;
%% pressure mass matrix
Mpp = sparse(M,M);
for j = 1:1:M
    Mpp(j,j) = betaf * phi(j) * h(j);
end
%% displacement-pressure stiffness matrix
Aup = Apu';
%% transmissibilities
T = zeros(M+1,1);
for j = 2:1:M
    T(j) = ( (h(j-1)/2)*(kappa(j-1)/viscosity)^(-1) + (h(j)/2)*(kappa(j)/viscosity)^(-1) )^(-1);
end
T(1) =  ( (h(1)/2)*(kappa(1)/viscosity)^(-1) )^(-1);
T(M+1) =  ( (h(M)/2)*(kappa(M)/viscosity)^(-1) )^(-1);
%% flux mass matrix using transmissibilities (trapezoidal rule)
Mqfqf = sparse(M+1,M+1);
for j = 1:1:M+1
    Mqfqf(j,j) = T(j)^(-1);
end
%% pressure-flux stiffness matrix
Apqf = sparse(M+1,M);
for j = 2:1:M
    Apqf(j,j-1) = 1;
    Apqf(j,j) = -1;
end
Apqf(1,1) = -1;
Apqf(M+1,M) = 1;
%% flux-pressure stiffness matrix
Aqfp = Apqf';
%% gravity term mechanical (BOM)
rho_avg = rhof * phi + rhos .* (1-phi);
G_vec_M = zeros(M+1,1);
for j = 2:1:M
   G_vec_M(j) = (1/2) * (h(j-1)*rho_avg(j-1) + h(j)*rho_avg(j)) * G;
end
G_vec_M(1) = (1/2) * h(1) * G * rho_avg(1);
G_vec_M(M+1) = (1/2) * h(M) * G * rho_avg(M);
%% gravity term Darcy (COM)
G_vec_H = zeros(M+1,1);
for j = 2:1:M
   G_vec_H(j) = (1/2) * (h(j-1) + h(j)) * G * rhof;
end
G_vec_H(1) = (1/2) * h(1) * G * rhof;
G_vec_H(M+1) = (1/2) * h(M) * G * rhof;
%%
%% boundary conditions, eliminate nodes from matrices
%%
%% displacement boundary conditions
Ua = exact_u(a,t,example);
Ub = exact_u(b,t,example);
%
tNa = -(lambda(1)+2*mu(1))*exact_du(a,t,example) + alpha*exact_p(a,t,example); % after taking product with normal i.e. = \sigma_total \cdot n
tNb = (lambda(M)+2*mu(M))*exact_du(b,t,example) - alpha*exact_p(b,t,example); % after taking product with normal i.e. = \sigma_total \cdot n
%% for example == 4 onwards (physical problem; external stress)
if (example >= 4)
    tNa = 0.1 + 0*t;
end
%% eliminate nodes
if (Neum_ua == 0 && Neum_ub == 0)
    free_nodes_u = 2:1:M;
elseif (Neum_ua == 0 && Neum_ub == 1)
    free_nodes_u = 2:1:M+1;
elseif (Neum_ua == 1 && Neum_ub == 0)
    free_nodes_u = 1:1:M;
else
    error('Cannot have only Neumann BC for displacement');
end
%% eliminate nodes from displacement matrix
Auu = Auu(free_nodes_u, free_nodes_u);     
Apu = Apu(free_nodes_u,:);
Aup = Apu';
G_vec_M = G_vec_M(free_nodes_u);
%% pressure, flux boundary conditions
Pa = exact_p(a,t,example);
Pb = exact_p(b,t,example);
qfa = exact_qf(a,t,example) * -1; % inner product with outward normal
qfb = exact_qf(b,t,example) * 1;
%
if (Neum_pa == 0 && Neum_pb == 0)
    free_nodes_qf = 1:1:M+1;
elseif (Neum_pa == 0 && Neum_pb == 1)
    free_nodes_qf = 1:1:M;
elseif (Neum_pa == 1 && Neum_pb == 0)
    free_nodes_qf = 2:1:M+1;
else
    free_nodes_qf = 2:1:M;
end
%% eliminate nodes from flux matrix
Mqfqf = Mqfqf(free_nodes_qf, free_nodes_qf);
Apqf = Apqf(free_nodes_qf,:);
Aqfp = Apqf';
G_vec_H = G_vec_H(free_nodes_qf);
%% block matrix
A = [Auu -alpha*Apu; -alpha*Aup -Mpp - tau*Aqfp*(Mqfqf \Apqf)];
%%
%% initial conditions
%%
previous_U = exact_u(xn(free_nodes_u),t(1),example);
previous_P = exact_p(xcc,t(1),example);
Q = zeros(M,1);
previous_etaf = initial_fluid_content(xcc,t(1),example);
settlement = zeros(length(t),1);
%% initialize error
linfty_l2_p = 0;
l2_l2_qf = 0;
linfty_H1_u = 0;
%% controls to toggle plotting and print errors
print_displacement = 1;
print_pressure = 1;
print_flux = 1;
if (example == 1 || example == 2 || example == 3)
    print_errors = 1;
elseif (example >= 4)
    print_errors = 0;
end
%%
%% time loop
%%
for n = 2:1:length(t)
    %%
    %% compute rhs_f vector (BOM source)
    %%
    rhs_f_val = rhs_f(xn,t(n),example);
    rhs_f_vector = zeros(M+1,1);
    for j = 2:1:M
        rhs_f_vector(j) = (h(j-1) + h(j))/2 * rhs_f_val(j);
    end
    rhs_f_vector(1) = (h(1)/2)*rhs_f_val(1);
    rhs_f_vector(M+1) = (h(M)/2)*rhs_f_val(M+1);
    %% eliminate boundary nodes
    rhs_f_vector = rhs_f_vector(free_nodes_u);
    %% displacement boundary contribution (Neumann and Dirichlet)
    rhs_f_vector(1) = rhs_f_vector(1) + (1-Neum_ua)*(1/h(1))*(lambda(1)+2*mu(1))*Ua(n) + Neum_ua*tNa(n);
    rhs_f_vector(end) = rhs_f_vector(end) + (1-Neum_ub)*(1/h(M))*(lambda(M)+2*mu(M))*Ub(n) + Neum_ub*tNb(n);
    %% gravity contribution
    rhs_f_vector = rhs_f_vector + G_vec_M;
    %%
    %% compute rhs_h vector (COM source)
    %%
    rhs_h_val = rhs_h(xcc,t(n),example);
    rhs_h_vector = zeros(M,1);
    for j = 1:1:M
        rhs_h_vector(j) = h(j)*rhs_h_val(j);
    end
    %% previous time step contribution
    rhs_h_vector = tau*rhs_h_vector + h .* previous_etaf;
    %% displacement boundary contribution
    rhs_h_vector(1) = rhs_h_vector(1) + alpha*(1-Neum_ua)*(Ua(n));
    rhs_h_vector(M) = rhs_h_vector(M) + alpha*(1-Neum_ub)*(- Ub(n));
    %% pressure boundary contribution and flux boundary contribution from qf = -\kappa \nabla p
    pressure_boundary = zeros(length(free_nodes_qf),1);
    pressure_boundary(1) = (1-Neum_pa)*Pa(n);
    pressure_boundary(end) = -(1-Neum_pb)*Pb(n);
    %
    rhs_h_vector = rhs_h_vector - tau*Aqfp*(Mqfqf \(pressure_boundary + G_vec_H) );
    %% flux Dirichlet boundary condition from COM equation \nabla q_f \cdot \eta_i
    rhs_h_vector(1) = rhs_h_vector(1) - Neum_pa * tau * qfa(n);
    rhs_h_vector(end) = rhs_h_vector(end) - Neum_pb * tau * qfb(n);
    %% compute total rhs vector [rhs_f; rhs_h]
    rhs_vector = zeros(length(free_nodes_u) + M,1);
    rhs_vector(1:length(free_nodes_u)) = rhs_f_vector;
    rhs_vector(length(free_nodes_u)+1:end) = -rhs_h_vector; % to make matrix symmetric
    %%
    %% compute solution at current time step
    %%
    UP = A \ rhs_vector;
    %%
    %% extract solution 
    %%
    U = UP(1:length(free_nodes_u)); 
    previous_U = U;
    P = UP(length(free_nodes_u)+1:end);
    previous_P = P;
    Q = Mqfqf \ (Apqf*P + pressure_boundary + G_vec_H);
    %% add Dirichlet boundary conditions to displacement vector
    U_BC = [Ua(n); zeros(M-1,1); Ub(n)];
    U_BC(free_nodes_u) = U;
    U = U_BC;
    %% add Dirichlet flux (Neumann pressure) boundary conditions to flux vector
    Q_BC = [-qfa(n); zeros(M-1,1); qfb(n)]; % -qfa(n) since qfa(n) is defined as product with normal so to "retrieve"
    % the one dimensional value of exact_qf at x = a
    Q_BC(free_nodes_qf) = Q;
    Q = Q_BC;
    %% update fluid content
    etaf = betaf * phi .* P;
    for j = 1:1:M
        etaf(j) = etaf(j) + alpha * (U(j+1) - U(j))/h(j);
    end
    previous_etaf = etaf;
    %% store settlement value (for Terzaghi's problem)
    settlement(n) = U(1);
    %%
    %% compute error (for example == 1, 2, or 3)
    %%
    if (print_errors == 1)
        %% l^\infty (l^2) error for pressure
        if (l2_err(xn,P,exact_p(xcc,t(n),example)) > linfty_l2_p)
            linfty_l2_p = max(linfty_l2_p,l2_err(xn,P,exact_p(xcc,t(n),example)));
            p_err_time = t(n);
        end
        %% l^2(l^2) error for flux
        Qcc = zeros(M,1);
        for j = 1:1:M
            Qcc(j) = (Q(j) + Q(j+1))/2;
        end
        l2_l2_qf = sqrt( l2_l2_qf^2 + tau*( l2_err(xn,Qcc,exact_qf(xcc,t(n),example)) )^2 );
        %% l^\infty (H1) error for displacement
        Ucc = zeros(M,1);
        for j = 1:1:M
            Ucc(j) = (U(j) + U(j+1))/2;
        end
        dUcc = zeros(M,1);
        for j = 1:1:M
            dUcc(j) = (U(j+1) - U(j))/h(j);
        end
        linfty_H1_u = max(linfty_H1_u, sqrt( l2_err(xn,Ucc,exact_u(xcc,t(n),example))^2 + l2_err(xn,dUcc,exact_du(xcc,t(n),example))^2 ));
        %%
    end
    %%
    %% plot solution
    %%
    if (abs(n - length(t)) < 1e-6)
        %% plot displacement
        if (print_displacement == 1)
        figure(1);
        if (example == 1 || example == 2 || example == 3)
            plot(xn,exact_u(xn,t(n),example),'-k','linewidth',2,'DisplayName','Exact');
        end
        hold on;
        plot(xn,U,'--*r','MarkerSize',12,'DisplayName','Numerical');
        hold off;
        xlim([floor(a) b]);
        xTickLocations = [floor(abs(xn(1))) (xn(1) + xn(end))/2 xn(end)];
        set(gca,'XTick', xTickLocations);
        for k = 1 : length(xTickLocations)
            xTickString{k} = xTickLocations(k);
        end
        set(gca,'XTickLabel', xTickString);
        xlabel('x [m]');
        ylabel('Displacement u [m]');
        if (example == 1 || example == 2 || example == 3)
            title(['t = ',num2str(t(n))]);
        elseif (example >= 4)
            title(['t = ',num2str(t(n)/24),' [day] ']);
        end
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        legend boxoff
        if (example >= 4)
            legend off;
        end
        box off;
        pause(0.1);
        end
        %% plot pressure
        if (print_pressure == 1)
        figure(2);
        if (example == 1 || example == 2 || example == 3)     
            plot(xcc,exact_p(xcc,t(n),example),'-k','linewidth',2,'DisplayName','Exact');
        end
        hold on;
        plot(xcc,P,'ob','MarkerSize',12,'DisplayName','Numerical');
        hold off;
        xlim([floor(a) b]);
        xTickLocations = [floor(abs(xn(1))) (xn(1) + xn(end))/2 xn(end)];
        set(gca,'XTick', xTickLocations);
        for k = 1 : length(xTickLocations)
            xTickString{k} = xTickLocations(k);
        end
        set(gca,'XTickLabel', xTickString);
        xlabel('x [m]');
        ylabel('Pressure p [MPa]'); 
        if (example == 1 || example == 2 || example == 3)
            title(['t = ',num2str(t(n))]);
        elseif (example >= 4)
            title(['t = ',num2str(t(n)/24),' [day] ']);
        end
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        legend boxoff
        if (example >= 4)
            legend off;
        end
        box off;
        pause(0.1);
        end
        %% plot flux 
        if (print_flux == 1)
        figure(3);      
        if (example == 1 || example == 2 || example == 3)
            plot(xn,exact_qf(xn,t(n),example),'-k','linewidth',2,'DisplayName','Exact');
        end
        hold on;
        plot(xn,Q,':*r','MarkerSize',12,'DisplayName','Numerical');
        hold off;
        xlim([floor(a) b]);
        xTickLocations = [floor(abs(xn(1))) (xn(1) + xn(end))/2 xn(end)];
        set(gca,'XTick', xTickLocations);
        for k = 1 : length(xTickLocations)
            xTickString{k} = xTickLocations(k);
        end
        set(gca,'XTickLabel', xTickString);
        xlabel('x [m]');
        ylabel('Flux q_f [m/hr]');
        if (example == 1 || example == 2 || example == 3)
            title(['t = ',num2str(t(n))]);
        elseif (example >= 4)
            title(['t = ',num2str(t(n)/24),' [day] ']);
        end
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        box off;
        if (example >= 4)
            legend off;
        end
        pause(0.1);
        end
    %%
    end
%%
end % time loop ends
%% plot total settlement
if (example >= 4)
    figure(4);
    hold on;
    plot(t,settlement,'-*','MarkerSize',12);
    hold off;
    xlim([0 Tend]);
    xTickLocations = [0 (0 + Tend)/(2) Tend];
    set(gca,'XTick', xTickLocations);
    for k = 1 : length(xTickLocations)
        xTickString{k} = floor(xTickLocations(k)/1);
    end
    set(gca,'XTickLabel', xTickString);
    xlabel('t [hr]');
    ylabel('Settlement [m]');
    title('Total settlement');
    box off;
    %% display maximum settlement
    % format long
     settlement(end)
end
%%
%% print errors
%%
if (print_errors == 1)
    fprintf('l_infty(l_2) error pressure p: %0.5g',linfty_l2_p);
    fprintf('\r');
    fprintf('l_2(l_2) error flux q_f: %0.5g',l2_l2_qf);
    fprintf('\r');
    fprintf('l_infty(H_1) error displacement u: %0.5g',linfty_H1_u);
    fprintf('\r');
end
%%
toc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%% MAIN FUNCTION ENDS
%%
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%% SPATIAL AND TEMPORAL DOMAIN
%% units: domain [m], time [hr]
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [a, b, Tend] = grid(example)
%%
if (example == 1 || example == 2 || example == 3)
    a = 0;
    b = 1;
    Tend = 1;
elseif (example == 4)
    a = 0;
    b = 0.1;
    Tend = 24;
elseif (example == 5)
    a = 0;
    b = 1;
    Tend = 8760;
else
    % custom values for a custom example
    a = 0;
    b = 1;
    Tend = 1;
end
%%
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%% PHYSICAL PARAMETERS
%% units: [m], [hr], [MPa]
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(xcc, example)
%%
if (example == 1 || example == 2 || example == 3)
    lambda = 2*0 + 1 + 0*xcc;
    mu = 3*0 + 1 + 0*xcc;
    kappa = 1e-1*0 + 1 + 0*xcc;
    phi = 0.5*0 + 1 + 0*xcc;
    viscosity = 2*0 + 1;
    alpha = 1;
    betaf = 5*0 + 1;
    rhof = 1;
    rhos = 1 + 0*xcc;
    G = 0.0;
elseif (example == 4)    
    E = 20;
    nu = 0.30;
    lambda = (E.*nu)./((1 + nu).*(1-2*nu)) + 0*xcc;
    mu = E./(2*(1+nu)) + 0*xcc;
    kappa = 1e-17 + 0*xcc;
    phi = 0.50 + 0*xcc;
    viscosity = 2.7822e-13;
    alpha = 1;
    betaf = 4.16e-4;
    rhof = 998.21 * (1e-6 * (1/3600)*(1/3600)); 
    rhos = 2700 * (1/3600) * (1/3600) * 1e-6 + 0*xcc;
    G = 1.27290528e8 * 1;
elseif (example == 5)
    lambda = 0*xcc;
    mu = 0*xcc;
    kappa = 0*xcc;
    phi = 0*xcc;
    rhos = 0*xcc;
    for j = 1:1:length(xcc)
        if (xcc(j) <= 0.5)
            E = 15;
            nu = 0.25;
            lambda(j) = (E*nu)/((1 + nu)*(1-2*nu));
            mu(j) = E/(2*(1+nu));
            kappa(j) = 1e-12;
            phi(j) = 0.30;
            rhos(j) = 2650 * (1/3600) * (1/3600) * 1e-6;
        else
            E = 20;
            nu = 0.30;
            lambda(j) = (E*nu)/((1 + nu)*(1-2*nu));
            mu(j) = E/(2*(1+nu));
            kappa(j) = 1e-17;
            phi(j) = 0.50;
            rhos(j) = 2700 * (1/3600) * (1/3600) * 1e-6;
        end
    end
    viscosity = 2.7822e-13;
    alpha = 1;
    betaf = 4.16e-4;
    rhof = 998.21 * (1e-6 * (1/3600)*(1/3600)); 
    G = 1.27290528e8 * 0; 
else
    %% custom scenario
    lambda = 1 + 0*xcc;
    mu = 1 + 0*xcc;
    kappa = 1 + 0*xcc;
    phi = 1 + 0*xcc;
    viscosity = 1;
    alpha = 1;
    betaf = 1;
    rhof = 1;
    rhos = 1 + 0*xcc;
    G = 0.0;
end
%%
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%% FUNCTION DEFINITONS
%%
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exact solution displacement
function y = exact_u(x,t,example)
%% 1.
if (example == 1)
    y = -sin(pi*t/2)*(1/(pi))*cos(pi*x);
%% 2.
elseif (example == 2)
    y = sin(pi*x/2).*exp(-t); 
%% 3.
elseif (example == 3)
    y = 2 - x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
%%
end
%% exact solution displacement derivative
function y = exact_du(x,t,example)
%% 1.
if (example == 1)
    y = sin(pi*t/2)*sin(pi*x);
%% 2.
elseif (example == 2)
    y = (pi/2)*cos(pi*x/2).*exp(-t);
%% 3.
elseif (example == 3)
    y = -1 + 0*x + 0*t;
%% 4.
else
    [lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
    y =  -(1/(lambda(1) + 2*mu(1))) + 0*x + 0*t; 
end
%%
end
%% exact solution pressure
function y = exact_p(x,t,example)
%% 1.
if (example == 1)
    y = sin(pi*t/2)*sin(pi*x);
%% 2.
elseif (example == 2)
    y = cos(pi*x/2).*exp(-t);
%% 3.
elseif (example == 3)
    y = 1 + x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
%%
end
%% exact solution flux
function y = exact_qf(x,t,example)
%% 1.
if (example == 1)
    [lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
    y = -sin(pi*t/2)*pi*cos(pi*x);
%% 2.
elseif (example == 2)
    [lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
    y =  (kappa(1)/viscosity).*(pi/2)*sin(pi*x/2).*exp(-t);
%% 3.
elseif (example == 3)
    y = -1 + 0*x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
%%
end
%% fluid content function
function y = initial_fluid_content(x,t,example)
%%
[lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
%% 1., 2., 3.
if (example == 1 || example == 2 || example == 3)
    y = alpha*exact_du(x,t,example) + betaf * phi .* exact_p(x,t,example);
%% 4.
else
    y = 0*x + 0*t + betaf*rhof*G*phi.*(x);
end
%%
end
%%
%% compute l2 error using grid norm 
function l2 = l2_err(x,val,exact_val)
l2 = 0;
for j = 1:1:length(x)-1
    h = x(j+1) - x(j);
    l2 = l2 + h*(val(j)-exact_val(j))^2;
end
l2 = sqrt(l2);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%% EXTERNAL SOURCES (BOM AND COM)
%%
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rhs f function (BOM)
function y = rhs_f(x,t,example)
%%
[lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
%%
%% 1.
if (example == 1)
    y = -pi * (lambda + 2*mu).*cos(pi*x)*sin(pi*t/2) + alpha*pi*sin(pi*t/2)*cos(pi*x);
%% 2.
elseif (example == 2)
    y = ( (lambda + 2*mu).*(pi*pi)/4 - alpha*pi/2 ).*sin(pi*x/2).*exp(-t);
%% 3.
elseif (example == 3)
    y = alpha + 0*x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% rhs h function (COM)
function y = rhs_h(x,t,example)
%%
[lambda, mu, alpha, kappa, viscosity, betaf, phi, rhof, rhos, G] = physical_parameters(x, example);
%%
%% 1.
if (example == 1)
    y = (pi/2)*cos(pi*t/2)*sin(pi*x).*(betaf * phi + alpha) + (1/viscosity)*sin(pi*t/2)*pi*pi*kappa.*sin(pi*x);
%% 2. 
elseif (example == 2)
    y = ( -betaf*phi - alpha*pi/2 + (kappa/viscosity)*(pi*pi)/4  ).*cos(pi*x/2).*exp(-t);
%% 3.
elseif (example == 3)
    y = 0*t + 0*x;
%% 4.
else
    y = 0*x + 0*t;
end
%%
end

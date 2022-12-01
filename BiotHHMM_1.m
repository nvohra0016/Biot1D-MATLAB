% Started on 5th April 2022. Biot 3 field formulation using different
% boundary conditions [MM;HH]. For example, [DN;ND] for Terzaghi's problem.
% Using non-uniform cell sizes.
% Using heterogeneous permeabilities.
% Input:
% el_nodes: if scalar then it is the number of cells and a uniform spatial grid is
% created. If it is a vector, then it contains the element nodes as a column vector (nodal grid), starting
% from x = a and including x = b, i.e., size(M) = #cells + 1. 
% tau: time step size.
function [xn, xcc, t, U, P, Q, linfty_H1_u, linfty_l2_p, l2_l2_qf, A, settlement] = BiotHHMM_1(el_nodes, ntsteps, BC_flags, example)
tic
%% set default plot values
set(groot,'defaultLineLineWidth',4)
set(0,'DefaultaxesLineWidth', 3)
set(0,'DefaultaxesFontSize', 24)
set(0,'DefaultTextFontSize', 18)
print_displacement = 1;
print_pressure = 1;
print_flux = 1;
print_errors = 0;
%% set boundary conditions. Dirichlet by default, so Neum == 1 for Neumann condition.
Neum_ua = BC_flags(1);
Neum_ub = BC_flags(2);
Neum_pa = BC_flags(3);
Neum_pb = BC_flags(4);
%% spatial and temporal grid
% spatial/temporal grid end points
a = 0; 
b = 1;
Tend = 1;
%% spatial grid size/ nodal/cell centered grid construction
if (size(el_nodes,1) == 1)
    M = el_nodes;
    h = (b-a)/M * ones(M,1);
    % nodal grid
    xn = a:h:b;
    xn = xn';
    % cell centered grid
    xcc = a+(h/2):h:b;
    xcc = xcc';
else
    M = size(el_nodes,1)-1;
    h = zeros(M,1);
    for j = 1:1:M
        h(j) = el_nodes(j+1)-el_nodes(j);
    end
    % nodal grid
    xn = el_nodes;
    % cell centered grid
    xcc = zeros(M,1);
    for j = 1:1:M
        xcc(j) = (el_nodes(j) + el_nodes(j+1))/2;
    end
end
tau = (Tend-0)/ntsteps
%% temporal grid
t = 0:tau:Tend;
t = t';
%% physical parameters
global lambda mu alpha kappa viscosity c0 betaf phi rhof G example_number
%
example_number = example;
if (example_number == 1 || example_number == 2 || example_number == 3)
    lambda = 1;
    mu = 1;
    kappa = 1;
    viscosity = 1;
    alpha = 1;
    phi = 1.0 * ones(M,1);
    betaf = 1.0;
    rhof = 1.0;
    rhos = 1.0;
    G = 0.0;
    lambda_vec = lambda*ones(M,1);    
    mu_vec = mu*ones(M,1);
    kappa_vec = ones(M,1);
    %
else
    % % sand
    E_sand = 15;
    nu_sand = 0.25;
    kappa_sand = 1e-12;
    phi_sand = 0.30;
    rhos_sand = 2650 * (1/3600) * (1/3600) * 1e-6; 
    %
    % clay
    E_clay = 20;
    nu_clay = 0.30;
    kappa_clay = 1e-17;
    phi_clay = 0.50;
    rhos_clay = 2700 * (1/3600) * (1/3600) * 1e-6;
    %
    % silt
    E_silt = 10;
    nu_silt = 0.35;
    kappa_silt = 1e-14;
    phi_silt = 0.45;
    rhos_silt = 2700 * (1/3600) * (1/3600) * 1e-6;
    %
    % heterogeneous E, nu/ choose soil
    %
    E = E_clay;
    nu = nu_clay;
    phi = phi_clay*ones(M,1);
    kappa = kappa_clay;
    %
    E_vec = E*ones(M,1);
    nu_vec = nu*ones(M,1);
    kappa_vec = kappa*ones(M,1);
    
    rhos = rhos_clay*ones(M,1);
    %
%     for j = 1:1:M
%         if ( (xcc(j) > 0) && (xcc(j) < 6))
%             E_vec(j) = E_silt*1 + 0*400;
%             nu_vec(j) = nu_silt;
%             phi(j) = phi_silt;
%             kappa_vec(j) = kappa_silt;
%             rhos(j) = rhos_silt;
%         else
%             E_vec(j) = E_clay*1000;
%             nu_vec(j) = nu_clay;
%             phi(j) = phi_clay;
%             kappa_vec(j) = kappa_clay;
%             rhos(j) = rhos_clay;
%         end
%     end
    %
    %E_vec(1:20) = 1e-1;
    lambda = (E*nu)/((1 + nu)*(1-2*nu));
    mu = E/(2*(1+nu));
    lambda_vec = (E_vec.*nu_vec)./((1 + nu_vec).*(1-2*nu_vec));
    mu_vec = E_vec./(2*(1+nu_vec));
    %
    % fluid properties
    alpha = 1;
    %viscosity = 1.0016e-9 * (1/3600);
    viscosity = 2.7822e-13;
    betaf = 4.16e-4;
    rhof = 998.21 * (1e-6 * (1/3600)*(1/3600)); % [kg/ m^3] to [MPa hr^2 / m^2]
    %
    G = 1.27290528e8 * 0; % [m / hr^2]
    %
    tNa_mag = 1e-1 * 1;
    
%     water_content = 1.0.*phi;
%     c = 5;
%     kappa_vec = (kappa_clay/viscosity) .* (water_content./phi).^c;
%     viscosity = 1;
    %tNa_mag*(b-a) /(lambda + 2*mu)
    %
%     lambda = 1;
%     mu = 1;
%     kappa = 1;
%     viscosity = 1;
%     alpha = 1;
%     phi = 1.0 * ones(M,1);
%     betaf = 1.0;
%     rhof = 1.0;
%     rhos = 1.0;
%     G = 0.0;
%     lambda_vec = lambda*ones(M,1);    
%     mu_vec = mu*ones(M,1);
%     kappa_vec = ones(M,1);
end
%% compute FE matrices
%% displacement stiffness matrix
Auu = sparse(M+1,M+1);
for j = 2:1:M
    Auu(j,j-1) = -1/h(j-1)*(lambda_vec(j-1) + 2*mu_vec(j-1));
    Auu(j,j) = (1/h(j-1))*(lambda_vec(j-1) + 2*mu_vec(j-1)) + (1/h(j))*(lambda_vec(j) + 2*mu_vec(j));
    Auu(j,j+1) = -1/h(j)*(lambda_vec(j) + 2*mu_vec(j));
end
Auu(1,1) = (1/h(1))*(lambda_vec(1) + 2*mu_vec(1)); Auu(1,2) = -1/h(1)*((lambda_vec(1) + 2*mu_vec(1)));
Auu(M+1,M) = (-1/h(M))*(lambda_vec(M) + 2*mu_vec(M)); Auu(M+1,M+1) = (1/h(M))*(lambda_vec(M) + 2*mu_vec(M));
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
% using variable kappa
% for j = 2:1:M
%     T(j) = (h*(viscosity/kappa(j)))^(-1);
% end
% T(1) =  ((h/2)*(viscosity/kappa(j)))^(-1);
% T(M+1) =  ((h/2)*(viscosity/kappa(j)))^(-1);
% using constant kappa
for j = 2:1:M
    T(j) = ( (h(j-1)/2)*(kappa_vec(j-1)/viscosity)^(-1) + (h(j)/2)*(kappa_vec(j)/viscosity)^(-1) )^(-1);
end
T(1) =  ( (h(1)/2)*(kappa_vec(1)/viscosity)^(-1) )^(-1);
T(M+1) =  ( (h(M)/2)*(kappa_vec(M)/viscosity)^(-1) )^(-1);
%% flux mass matrix using trapezoidal rule
% Mqfqf = sparse(M+1,M+1);
% for j = 2:1:M
%     Mqfqf(j,j) = (h/2)*(kappa/viscosity)^(-1) + (h/2)*(kappa/viscosity)^(-1);
% end
% Mqfqf(1,1) = (h/2)*(kappa/viscosity)^(-1);
% Mqfqf(M+1,M+1) = (h/2)*(kappa/viscosity)^(-1);

% this is same as above i.e. trapezoidal rule
MqfqfCCFD = sparse(M+1,M+1);
for j = 1:1:M+1
    MqfqfCCFD(j,j) = T(j)^(-1);
end
%% flux mass matrix using exact integration
MqfqfSP = sparse(M+1,M+1);
for j =2:1:M
    MqfqfSP(j,j) = (h(j-1)/3) * ( (viscosity/kappa_vec(j-1))) + (h(j)/3) * ( (viscosity/kappa_vec(j)) );
    MqfqfSP(j,j+1) = (h(j)/6) * ( (viscosity/kappa_vec(j)) );
    MqfqfSP(j,j-1) = (h(j-1)/6) * ( (viscosity/kappa_vec(j-1)) );
end
MqfqfSP(1,1) = (h(1)/3) * ((viscosity/kappa_vec(1))); MqfqfSP(1,2) = (h(1)/6) * ( (viscosity/kappa_vec(1)) );
MqfqfSP(M+1,M) = (h(M)/6) * ( (viscosity/kappa_vec(M)) ); MqfqfSP(M+1,M+1) = (h(M)/3) * ( (viscosity/kappa_vec(M)) );
%% choose SP or CCFD
%Mqfqf = MqfqfSP;
Mqfqf = MqfqfCCFD;
%% calculate Mqfqf CCFD inverse while checking if permeability is 0
Mqfqfinv = sparse(M+1,M+1);
for j = 2:1:M
    if (kappa_vec(j) ~= 0 && kappa_vec(j-1) ~= 0)
        Mqfqfinv(j,j) = (1/viscosity) * 2*kappa_vec(j-1) * kappa_vec(j) / (kappa_vec(j) * h(j-1) + kappa_vec(j-1) * h(j));
    else
        Mqfqfinv(j,j) = 0;
    end
end
Mqfqfinv(1,1) = (2/h(1)) * (kappa_vec(1)/viscosity);
Mqfqfinv(M+1,M+1) = (2/h(M)) * (kappa_vec(M)/viscosity);

%Mqfqfinv = inv(Mqfqf);
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
%% gravity term mechanical
rho_avg = rhof * phi + rhos .* (1-phi);
G_vec_M = zeros(M+1,1);
for j = 2:1:M
   G_vec_M(j) = (1/2) * (h(j-1)*rho_avg(j-1) + h(j)*rho_avg(j)) * G;
end
G_vec_M(1) = (1/2) * h(1) * G * rho_avg(1);
G_vec_M(M+1) = (1/2) * h(M) * G * rho_avg(M);
%% gravity term Darcy
G_vec_H = zeros(M+1,1);
for j = 2:1:M
   G_vec_H(j) = (1/2) * (h(j-1) + h(j)) * G * rhof;
end
G_vec_H(1) = (1/2) * h(1) * G * rhof;
G_vec_H(M+1) = (1/2) * h(M) * G * rhof;
%% boundary conditions, eliminate nodes from matrices
%% displacement boundary conditions
Ua = exact_u(a,t);
Ub = exact_u(b,t);

tNa = -(lambda_vec(1)+2*mu_vec(1))*exact_du(a,t) + alpha*exact_p(a,t); % after taking product with normal i.e. = \sigma_total \cdot n
tNb = (lambda_vec(M)+2*mu_vec(M))*exact_du(b,t) - alpha*exact_p(b,t); % after taking product with normal i.e. = \sigma_total \cdot n

%% for example == 4 onwards
%tNa = tNa_mag * ones(length(t),1);

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
%Aup = Aup(free_nodes,free_nodes);
Aup = Apu';
G_vec_M = G_vec_M(free_nodes_u);
%% pressure, flux boundary conditions
Pa = exact_p(a,t);
Pb = exact_p(b,t);
qfa = exact_qf(a,t) * -1; % inner product with outward normal
qfb = exact_qf(b,t) * 1;
%
qfa_bdry = Mqfqf(2,1);
qfb_bdry = Mqfqf(M,M+1);
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
Mqfqfinv = Mqfqfinv(free_nodes_qf, free_nodes_qf);
Apqf = Apqf(free_nodes_qf,:);
Aqfp = Apqf';
G_vec_H = G_vec_H(free_nodes_qf);
%% block matrix
%A = [Auu -alpha*Apu; -alpha*Aup -Mpp - tau*Aqfp*(Mqfqf \Apqf)];
A = [Auu -alpha*Apu; -alpha*Aup -Mpp - tau*Aqfp*(Mqfqfinv * Apqf)];
%% matrices to analyse
%% initial conditions
previous_U = exact_u(xn(free_nodes_u),t(1));
previous_P = exact_p(xcc,t(1));
Q = zeros(M,1);
etaf = initial_fluid_content(xcc,t(1));
settlement = zeros(length(t),1);
%% initialize error
linfty_l2_p = 0;
l2_l2_qf = 0;
linfty_H1_u = 0;
%% time loop
for n = 2:1:length(t)
    %% compute rhs_f vector
    rhs_f_val = rhs_f(xn,t(n));
    rhs_f_vector = zeros(M+1,1);
    for j = 2:1:M
        rhs_f_vector(j) = (h(j-1) + h(j))/2 * rhs_f_val(j);
    end
    rhs_f_vector(1) = (h(1)/2)*rhs_f_val(1);
    rhs_f_vector(M+1) = (h(M)/2)*rhs_f_val(M+1);
    %% eliminate boundary nodes
    rhs_f_vector = rhs_f_vector(free_nodes_u);
    %rhs_f_vector = l2_prod_f(a,b,M,tau,t(n));
    %% displacement boundary contribution (Neumann and Dirichlet)
    rhs_f_vector(1) = rhs_f_vector(1) + (1-Neum_ua)*(1/h(1))*(lambda_vec(1)+2*mu_vec(1))*Ua(n) + Neum_ua*tNa(n);
    rhs_f_vector(end) = rhs_f_vector(end) + (1-Neum_ub)*(1/h(M))*(lambda_vec(M)+2*mu_vec(M))*Ub(n) + Neum_ub*tNb(n);
    %% gravity contribution
    rhs_f_vector = rhs_f_vector + G_vec_M;
    %% compute rhs_h vector
    rhs_h_val = rhs_h(xcc,t(n));
    rhs_h_vector = zeros(M,1);
    for j = 1:1:M
        rhs_h_vector(j) = h(j)*rhs_h_val(j);
    end
    %rhs_h_vector = l2_prod_h(a,b,M,tau,t(n));
    %% previous time step contribution
    rhs_h_vector = tau*rhs_h_vector + h .* etaf;
    %% displacement boundary contribution
    rhs_h_vector(1) = rhs_h_vector(1) + alpha*(1-Neum_ua)*(Ua(n));
    rhs_h_vector(M) = rhs_h_vector(M) + alpha*(1-Neum_ub)*(- Ub(n));
    %% pressure boundary contribution and flux boundary contribution from qf = -\kappa \nabla p
    pressure_boundary = zeros(length(free_nodes_qf),1);
    pressure_boundary(1) = (1-Neum_pa)*Pa(n);
    pressure_boundary(end) = -(1-Neum_pb)*Pb(n);
    %
    % this is primarily for SP formulation
    flux_boundary = zeros(length(free_nodes_qf),1); 
    flux_boundary(1) = -Neum_pa * qfa_bdry * (-qfa(n)); % qfa is flux \cdot n i.e. -qfa is the boundary value 
    flux_boundary(end) = -Neum_pb * qfb_bdry * qfb(n);
    %
    rhs_h_vector = rhs_h_vector - tau*Aqfp*(Mqfqf \(pressure_boundary + flux_boundary + G_vec_H) );
    %% flux Dirichlet boundary condition from COM equation \nabla q_f \cdot \eta_i
    rhs_h_vector(1) = rhs_h_vector(1) - Neum_pa * tau * qfa(n);
    rhs_h_vector(end) = rhs_h_vector(end) - Neum_pb * tau * qfb(n);
    %% compute total rhs vector
    rhs_vector = zeros(length(free_nodes_u) + M,1);
    rhs_vector(1:length(free_nodes_u)) = rhs_f_vector;
    rhs_vector(length(free_nodes_u)+1:end) = -rhs_h_vector; % make matrix symmetric
    %% compute solution at current time step
    UP = A \ rhs_vector;
    %UP = inv(A) * rhs_vector;
    %% extract solution 
    U = UP(1:length(free_nodes_u)); 
    previous_U = U;
    P = UP(length(free_nodes_u)+1:end);
    previous_P = P;
    %Q = Mqfqf \ (Apqf*P + pressure_boundary + flux_boundary + G_vec_H);
    Q = Mqfqfinv * (Apqf*P + pressure_boundary + flux_boundary + G_vec_H);
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
    %% store settlement value
    settlement(n) = U(1);
    %% compute error
    %% l^\infty (l^2) error for pressure
    if (l2_err(xn,P,exact_p(xcc,t(n))) > linfty_l2_p)
      linfty_l2_p = max(linfty_l2_p,l2_err(xn,P,exact_p(xcc,t(n))));
      p_err_time = t(n);
    end
    %% l^2(l^2) error for flux
    Qcc = zeros(M,1);
    for j = 1:1:M
        Qcc(j) = (Q(j) + Q(j+1))/2;
    end
    l2_l2_qf = sqrt( l2_l2_qf^2 + tau*( l2_err(xn,Qcc,exact_qf(xcc,t(n))) )^2 );
    %% l^\infty (H1) error for displacement
%     Ucc = zeros(M,1);
%     for j = 2:1:M
%         Ucc(j) = (U(j) + U(j-1))/2;
%     end
%     dUcc = zeros(M,1);
%     for j = 1:1:M
%         dUcc(j) = (U(j+1) - U(j))/h(j);
%     end
    %linfty_H1_u = max(linfty_H1_u, sqrt( l2_err(xn,Ucc,exact_u(xcc,t(n)))^2 + l2_err(xn,dUcc,exact_du(xcc,t(n)))^2 ));
    %% plot solution
    if (n == length(t) || n == floor(length(t)/2) || n == 2)
        t(n)/24
        %% plot displacement
        if (print_displacement == 1)
        figure(1);
        plot(xn,exact_u(xn,t(n)),'-k','linewidth',2,'DisplayName','Exact');
        hold on;
        plot(xn,U,'--*','linewidth',3,'DisplayName','Numerical');
        hold off;
%         hold on;
%         plot(U,xn,'--*','linewidth',3,'DisplayName','Numerical');
%         hold off;
        %xlim([0 0.012]);
        xlim([floor(a) b]);
        xTickLocations = [floor(abs(xn(1))) (xn(1) + xn(end))/2 xn(end)];
        %yTickLocations = [0 4 6 10];
        set(gca,'XTick', xTickLocations);
        for k = 1 : length(xTickLocations)
            xTickString{k} = xTickLocations(k);
        end
        set(gca,'XTickLabel', xTickString);
        xlabel('Depth x [m]');
        ylabel('Displacement u [m]');
        %set(gca,'YDir','reverse');
        %title(['t = ',num2str(t(n)),', M = ',num2str(M),', tau = ',num2str(tau)],'FontSize',18);
        title(['t = ',num2str(t(n)/24),' [day] '],'FontSize',24);
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        legend boxoff
        legend off;
        box off;
        pause(0.1);
        end
        %% plot pressure
        if (print_pressure == 1)
        figure(2);
        plot(xcc,exact_p(xcc,t(n)),'-k','linewidth',2,'DisplayName','Exact');
        hold on;
        plot(xcc,P,'o','linewidth',3,'MarkerSize',12,'DisplayName','Numerical');
        hold off;
%         hold on;
%         plot(P,xcc,'o','linewidth',3,'MarkerSize',12,'DisplayName','Numerical');
%         hold off;
        xlim([floor(a) b]);
        xTickLocations = [floor(abs(xn(1))) (xn(1) + xn(end))/2 xn(end)];
        %yTickLocations = [0 4 6 10];
        set(gca,'XTick', xTickLocations);
        for k = 1 : length(xTickLocations)
            xTickString{k} = xTickLocations(k);
        end
        set(gca,'XTickLabel', xTickString);
        xlabel('Depth x [m]');
        ylabel('Pressure p [MPa]'); 
        %set(gca, 'YDir','reverse')
        %title(['t = ',num2str(t(n)),', M = ',num2str(M),', tau = ',num2str(tau)],'FontSize',18);
        title(['t = ',num2str(t(n)/24),' [day] '],'FontSize',24);
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        legend boxoff
        legend off;
        box off;
        pause(0.1);
        end
        %% plot flux 
        if (print_flux == 1)
        figure(3);
        plot(xn,exact_qf(xn,t(n)),'-k','linewidth',2,'DisplayName','Exact');
        hold on;
        plot(xn,Q,':*','linewidth',3,'DisplayName','Numerical');
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
        title(['t = ',num2str(t(n)),', M = ',num2str(M),', tau = ',num2str(tau)],'FontSize',18);
        %% legend properties
        lh = legend;
        set(lh,'FontSize',24,'location','southeast');
        box off;
        legend off;
        pause
        end
    %%
    end
    %% print values
%     fprintf('Flux q_f at x = 0: %0.5g',Q(1));
%     fprintf(', ');
%     fprintf('Pressure p at x = h/2: %0.5g',P(1));
%     fprintf('\r');
end
%% plot total settlement
if (example_number == 4)
    figure(4);
    hold on;
    plot(t,settlement,'-','linewidth',2,'DisplayName','With gravity');
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
    title(['M = ',num2str(M),', tau = ',num2str(tau)],'FontSize',18);
    box off;
end

%% print errors
if (print_errors == 1)
    fprintf('l_infty(l2) error pressure p: %0.5g',linfty_l2_p);
    fprintf('\r');
    fprintf(', at time: %0.5g',p_err_time);
    fprintf('\r');
    fprintf('l_2(l2) error flux q_f: %0.5g',l2_l2_qf);
    fprintf('\r');
    fprintf('l_infty(H1) error displacement u: %0.5g',linfty_H1_u);
    fprintf('\r');
else
    linfty_H1_u = 0;
    linfty_l2_p = 0;
    l2_l2_qf = 0;
end
%%
toc
end
%% compute l2 error using grid norm 
function l2 = l2_err(x,val,exact_val)
l2 = 0;
for j = 1:1:length(x)-1
    h = x(j+1) - x(j);
    l2 = l2 + h*(val(j)-exact_val(j))^2;
end
l2 = sqrt(l2);
end
%% exact solution displacement
function y = exact_u(x,t)
global lambda mu alpha kappa kappa_plot viscosity c0 example_number
%% 1.
if (example_number == 1)
    y = -sin(pi*t/2)*(1/(pi))*cos(pi*x);
%% 2.
elseif (example_number == 2)
    y = sin(pi*x/2).*exp(-t); 
%% 3.
elseif (example_number == 3)
    y = 2 - x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% exact solution displacement derivative
function y = exact_du(x,t)
global lambda mu alpha kappa kappa_plot viscosity c0 example_number
%% 1.
if (example_number == 1)
    y = sin(pi*t/2)*sin(pi*x);
%% 2.
elseif (example_number == 2)
    y = (pi/2)*cos(pi*x/2).*exp(-t);
%% 3.
elseif (example_number == 3)
    y = -1 + 0*x + 0*t;
%% 4.
else
    y =  -(1/(lambda + 2*mu)) + 0*x + 0*t; 
end
end
%% exact solution pressure
function y = exact_p(x,t)
global rhof G phi example_number
%% 1.
if (example_number == 1)
    y = sin(pi*t/2)*sin(pi*x);
%% 2.
elseif (example_number == 2)
    y = cos(pi*x/2).*exp(-t);
%% 3.
elseif (example_number == 3)
    y = 1 + x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% exact solution flux
function y = exact_qf(x,t)
global lambda mu alpha kappa kappa_plot viscosity c0 example_number
%% 1.
if (example_number == 1)
    y = -(1/viscosity)*sin(pi*t/2)*pi*kappa*cos(pi*x);
%% 2.
elseif (example_number == 2)
    y =  (pi/2)*(kappa/viscosity)*sin(pi*x/2).*exp(-t);
%% 3.
elseif (example_number == 3)
    y = -kappa/viscosity + 0*x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% fluid content function
function y = initial_fluid_content(x,t)
global alpha phi betaf rhof G example_number
%% 1., 2., 3.
if (example_number == 1 || example_number == 2 || example_number == 3)
    y = alpha*exact_du(x,t) + betaf * phi .* exact_p(x,t);
%% 4.
else
    y = 0*x + 0*t + betaf*rhof*G*phi.*(x);
end
end
%% rhs f function (BOM)
function y = rhs_f(x,t)
global lambda mu alpha kappa kappa_plot viscosity c0 example_number
%% 1.
if (example_number == 1)
    y = -(lambda + 2*mu)*(pi)*sin(pi*t/2)*cos(pi*x) + alpha*pi*sin(pi*t/2)*cos(pi*x);
%% 2.
elseif (example_number == 2)
    y = ( (lambda + 2*mu)*(pi*pi)/4 - alpha*pi/2 )*sin(pi*x/2).*exp(-t);
%% 3.
elseif (example_number == 3)
    y = alpha + 0*x + 0*t;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% rhs h function (COM)
function y = rhs_h(x,t)
global lambda mu alpha kappa kappa_plot viscosity c0 betaf phi example_number
%% 1.
if (example_number == 1)
    y = (pi/2)*cos(pi*t/2)*sin(pi*x).*(betaf * phi + alpha) + (1/viscosity)*sin(pi*t/2)*pi*pi*kappa*sin(pi*x);
%% 2. 
elseif (example_number == 2)
    y = ( -betaf*phi - alpha*pi/2 + (kappa/viscosity)*(pi*pi)/4  ).*cos(pi*x/2).*exp(-t);
%% 3.
elseif (example_number == 3)
    y = 0*t + 0*x;
%% 4.
else
    y = 0*x + 0*t;
end
end
%% L2 projection of exact pressure onto M_h i.e. piecewise constants
function y = L2_proj_exact_p(a,b,M,t)
%% quadrature points
% n = 3;
% points = [-sqrt(3/5) 0 sqrt(3/5)]';
% weights = [5/9 8/9 5/9]';
n = 1;
points = [0]';
weights = [2]';
h = (b-a)/M;
xn = a:h:b;
integral_f_eta = zeros(M,1);
%% create physical quadrature points and compute \int (f \eta_j), 1 \leq j \leq M
for j = 1:1:length(xn)-1
    phys_points = xn(j)*(1/2)*(1 - points) + xn(j+1)*(1/2)*(1 + points);
    integral_f_eta(j) =  sum(exact_p(phys_points,t) .* weights * (h/2));
end
y = (1/h)*integral_f_eta;
end
%% calculate inner product of function and linear basis elements on cell
function y = l2_prod_f(a,b,M,tau,t)
%% quadrature points
% n = 3;
% points = [-sqrt(3/5) 0 sqrt(3/5)]';
% weights = [5/9 8/9 5/9]';
n = 1;
points = 0;
weights = 2;
h = (b-a)/M;
xn = a:h:b;
integral_f_psi = zeros(M,1);
%% create physical quadrature points and compute \int (f \phi_j), 1 \leq j \leq M
for j = 2:1:M
    phys_points = xn(j-1)*(1/2)*(1 - points) + xn(j)*(1/2)*(1 + points);
    psi_j_val = (phys_points - xn(j-1))*(1/h);
    integral_f_psi(j) =  sum(rhs_f(phys_points,t).* psi_j_val .* weights * (h/2));
    
    phys_points = xn(j)*(1/2)*(1 - points) + xn(j+1)*(1/2)*(1 + points);
    psi_j_val = (phys_points - xn(j+1))*(-1/h);
    integral_f_psi(j) = integral_f_psi(j) + sum(rhs_f(phys_points,t).* psi_j_val .* weights * (h/2));
end
%% uncomment for boundary elements
phys_points = xn(1)*(1/2)*(1 - points) + xn(2)*(1/2)*(1 + points);
psi_j_val = (phys_points - xn(2))*(-1/h);
integral_f_psi(1) =  sum(rhs_f(phys_points,t).* psi_j_val .* weights * (h/2));

% phys_points = xn(M)*(1/2)*(1 - points) + xn(M+1)*(1/2)*(1 + points);
% psi_j_val = (phys_points - xn(M))*(1/h);
% integral_f_psi(M+1) =  sum(rhs_f(phys_points,t).* psi_j_val .* weights * (h/2));
%%
y = integral_f_psi;
%%
end
%% calculate inner product of function and piecewise constant basis elements on cell
function y = l2_prod_h(a,b,M,tau,t)
%% quadrature points
% n = 3;
% points = [-sqrt(3/5) 0 sqrt(3/5)]';
% weights = [5/9 8/9 5/9]';
n = 1;
points = [0]';
weights = [2]';
h = (b-a)/M;
xn = a:h:b;
integral_h_eta = zeros(M,1);
%% create physical quadrature points and compute \int (f \eta_j), 1 \leq j \leq M
for j = 1:1:M
    phys_points = xn(j)*(1/2)*(1 - points) + xn(j+1)*(1/2)*(1 + points);
    eta_j_val = 1 + 0*phys_points;
    integral_h_eta(j) =  sum(rhs_h(phys_points,t).* eta_j_val .* weights * (h/2));
end
y = integral_h_eta;
%%
end
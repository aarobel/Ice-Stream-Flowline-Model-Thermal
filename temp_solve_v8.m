function T_out = temp_solve_v8(T,h,x_g,u,v,parameters)
%temp_solve uses a backward euler
%algorithm to find the x-z temperature of an ice stream
%with specified velocity, thickness and grounding line position
%Written by Alex Robel, last modified November 14, 2013

%generate LHS matrix of temp coefficients and RHS vector of constants
[tempLHS,tempRHS] = temp_generate(T,h,x_g,u,v,parameters);

%apply Dirichlet boundary conditions
n_total = parameters.grid.n_elements*parameters.grid.n2_elements;

basal_ind = 1:parameters.grid.n_elements;
sfc_ind = n_total-parameters.grid.n_elements+1:n_total;

tempLHS(basal_ind,:) = 0; %zero out rows for basal layer
tempLHS(sfc_ind,:) = 0; %zero out rows for surface layer

tempLHS(basal_ind,:) = spdiags([1.5*ones(parameters.grid.n_elements,1), -0.5*ones(parameters.grid.n_elements,1)],...
    [0,parameters.grid.n_elements], parameters.grid.n_elements,n_total);
tempLHS(sfc_ind,:) = spdiags([-0.5*ones(parameters.grid.n_elements,1), 1.5*ones(parameters.grid.n_elements,1)],...
    [n_total-2.*parameters.grid.n_elements,n_total-parameters.grid.n_elements], parameters.grid.n_elements,n_total);

tempRHS(basal_ind) = zeros(parameters.grid.n_elements,1); %set dirichlet conditions for upper boundary
tempRHS(sfc_ind) = parameters.T_s; %set dirichlet conditions for lower boundary

%solve
T_out = tempLHS\tempRHS;

%reconstitute
T_out = wrap(T_out',parameters.grid.n_elements);
end

function [foutLHS,foutRHS] = temp_generate(T_old,h,x_g,u,v,parameters)
%note that the input to this function hx_g is the
%ice thickness column vector concatenated with the
%grounding line position in the last row

h_old = parameters.hx_g_old(1:end-1);
x_g_old = parameters.hx_g_old(end);

%% unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1), n_nodes-by-one column vector
sigma_node = repmat(sigma_node,[1 parameters.grid.n2_nodes]);
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1), n_elements-by-one column vector
sigma_element = repmat(sigma_element,[1 parameters.grid.n2_elements]);
dsigma = sigma_node(2:n_nodes,:)-sigma_node(1:n_nodes-1,:); %spacing between nodes

n2_nodes = parameters.grid.n2_nodes;      %number of eta node positions
n2_elements = parameters.grid.n2_elements;  %number of eta finite elements (= n2_nodes-1 in 1-D), h and N have length n_elements
eta_node = parameters.grid.eta_node;     %eta node positions scaled to (0,1), n2_nodes-by-one column vector
eta_node = repmat(eta_node',[parameters.grid.n_nodes 1]);
eta_element = parameters.grid.eta_element;  %eta element centres scaled to (0,1), n2_elements-by-one column vector
eta_element = repmat(eta_element',[parameters.grid.n_elements 1]);
deta = eta_node(:,2:n2_nodes)-eta_node(:,1:n2_nodes-1); %spacing between eta nodes                      %current w ice velocity

dtau = parameters.dtau;

%% indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_n_int = (2:n_nodes-1)';

ind_e_minus = (1:n_elements-1)';
ind_e_plus = (2:n_elements)';

ind_n2_minus = (1:n2_nodes-1)';
ind_n2_plus = (2:n2_nodes)';
ind_n2_int = (2:n2_nodes-1)';

ind_e2_minus = (1:n2_elements-1)';
ind_e2_plus = (2:n2_elements)';

%% unpack other parameters
h_g = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
T_s = parameters.T_s;
x_g_deriv = (x_g - x_g_old)/dtau;

%% calculate nodal interior terms
u_elem = (u(1:n_nodes-1,:) + u(2:n_nodes,:))./2;
u_elem = (u_elem(:,1:n2_nodes-1) + u_elem(:,2:n2_nodes))./2;
h_elem = repmat(h,[1 parameters.grid.n2_elements]);
h_old_elem = repmat(h_old,[1 parameters.grid.n2_elements]);
h_deriv_elem = (h_elem-h_old_elem)./parameters.dtau;
bedHflux = (parameters.bedHflux(1:end-1)+parameters.bedHflux(2:end))./2;

dhdsigma                   = zeros(n_elements,n2_elements);
dhdsigma(2:n_elements-1,:) = (h_elem(3:n_elements,:) - h_elem(1:n_elements-2,:))./(sigma_element(3:n_elements,:)-sigma_element(1:n_elements-2,:));
dhdsigma(1,:)              = (h_elem(2,:) - h_elem(1,:))./(dsigma(1,1:n2_elements));
dhdsigma(n_elements,:)     = (h_g.*ones(1,n2_elements)-h_elem(n_elements-1,:))./(sigma_node(n_nodes,1:n2_elements)-sigma_element(n_elements-1,:));

dbdsigma                   = zeros(n_elements,n2_elements);
dbdsigma(2:n_elements-1,:) = (Base(x_g.*sigma_element(3:n_elements,:),parameters) - Base(x_g.*sigma_element(1:n_elements-2,:),parameters))./(sigma_element(3:n_elements,:)-sigma_element(1:n_elements-2,:));
dbdsigma(1,:)              = (Base(x_g.*sigma_node(2,1:n2_elements),parameters) - Base(x_g.*sigma_node(1,1:n2_elements),parameters))./(sigma_node(2,1:n2_elements)-sigma_node(1,1:n2_elements));
dbdsigma(n_elements,:)     = (Base(x_g.*sigma_node(n_nodes,1:n2_elements),parameters) - Base(x_g.*sigma_node(n_nodes-1,1:n2_elements),parameters))./(sigma_node(n_nodes,1:n2_elements)-sigma_node(n_nodes-1,1:n2_elements));

u_int = (u(ind_n_plus,1:n2_nodes-1)+u(ind_n_plus,2:n2_nodes))./2;
sigma_int = (sigma_node(ind_n_plus,1:n2_nodes-1)+sigma_node(ind_n_plus,2:n2_nodes))./2;
theta_int = theta(x_g,x_g_old,u_int,sigma_int,parameters);

u_plus = (u(ind_n_plus,1:n2_nodes-1)+u(ind_n_plus,2:n2_nodes))./2;
sigma_plus = (sigma_node(ind_n_plus,1:n2_nodes-1)+sigma_node(ind_n_plus,2:n2_nodes))./2;
theta_plus = theta(x_g,x_g_old,u_plus,sigma_plus,parameters);

u_minus = (u(ind_n_minus,1:n2_nodes-1)+u(ind_n_minus,2:n2_nodes))./2;
sigma_minus = (sigma_node(ind_n_minus,1:n2_nodes-1)+sigma_node(ind_n_minus,2:n2_nodes))./2;
theta_minus = theta(x_g,x_g_old,u_minus,sigma_minus,parameters);

v_int = v(:,2:n2_nodes);
eta_int = (eta_node(1:n_nodes-1,2:n2_nodes)+eta_node(2:n_nodes,2:n2_nodes))./2;
theta_z_int = theta_z(v_int - eta_int.*(h_deriv_elem +...
    (dhdsigma.*(u_int-(sigma_element.*x_g_deriv))./x_g)) -...
    ((u_int.*dbdsigma)./x_g),parameters);

% v_plus = (v(1:n_nodes-1,2:n2_nodes)+v(2:n_nodes,2:n2_nodes))./2;
v_plus = v(:,2:n2_nodes);
eta_plus = (eta_node(1:n_nodes-1,2:n2_nodes)+eta_node(2:n_nodes,2:n2_nodes))./2;
theta_z_plus = theta_z(v_plus - eta_plus.*(h_deriv_elem +...
    (dhdsigma.*(u_plus-(sigma_element.*x_g_deriv))./x_g)) -...
    ((u_plus.*dbdsigma)./x_g),parameters);

% v_minus = (v(1:n_nodes-1,1:n2_nodes-1)+v(2:n_nodes,1:n2_nodes-1))./2;
v_minus = v(:,1:n2_nodes-1);
eta_minus = (eta_node(1:n_nodes-1,1:n2_nodes-1)+eta_node(2:n_nodes,1:n2_nodes-1))./2;
theta_z_minus = theta_z(v_minus - eta_minus.*(h_deriv_elem +...
    (dhdsigma.*(u_minus-(sigma_element.*x_g_deriv))./x_g)) -...
    ((u_minus.*dbdsigma)./x_g),parameters);

%%
v_eff_plus = ((v_int(:,1:n2_elements) -...
    eta_int(:,1:n2_elements).*(h_deriv_elem(:,1:n2_elements)+...
    (dhdsigma(:,1:n2_elements)./x_g).*(u_elem(:,1:n2_elements)-sigma_element(:,1:n2_elements).*x_g_deriv)) -...
    (u_elem(:,1:n2_elements).*dbdsigma(:,1:n2_elements)./x_g))./...
    h_elem(:,1:n2_elements))./deta(1:n_elements,1:n2_elements);

Vadv_minus = zeros(n_elements,n2_elements); %ensure that first term is zero (though it does go to zero anyway)
v_eff_minus(:,2:n2_elements) = ((v_int(:,1:n2_elements-1) - ...
    eta_int(:,1:n2_elements-1).*(h_deriv_elem(:,2:n2_elements)+...
    (dhdsigma(:,2:n2_elements)./x_g).*(u_elem(:,2:n2_elements)-sigma_element(:,2:n2_elements).*x_g_deriv)) -...
    (u_elem(:,2:n2_elements).*dbdsigma(:,2:n2_elements)./x_g))./...
    h_elem(:,2:n2_elements))./deta(1:n_elements,2:n2_elements);

%% (i,j+1) temperature term
dFdT_plus = zeros(n_elements,n2_elements); %ensure that last term is zero
%u plus advection term
dFdT_plus(2:n_elements,:) = ((u_int(1:n_elements-1,:) - sigma_int(1:n_elements-1,:).*x_g_deriv)./x_g).*(1-theta_int(1:n_elements-1,:))./...
    dsigma(ind_e_minus,1:n2_elements);

%% (i,j-1) temperature term
dFdT_minus = zeros(n_elements,n2_elements); %ensure that first term is zero
%u minus advection term
dFdT_minus(1:n_elements-1,:) = -((u_int(1:n_elements-1,:) - sigma_int(1:n_elements-1,:).*x_g_deriv)./x_g).*theta_int(1:n_elements-1,:)./...
    dsigma(ind_e_plus,1:n2_elements);

%% (i,j) temperature term
dFdT_zero = zeros(n_elements,n2_elements);

%advection terms

dFdT_zero = ((u_plus - sigma_plus.*x_g_deriv)./x_g).*theta_plus./dsigma(:,1:n2_elements) -...
            ((u_minus - sigma_minus.*x_g_deriv)./x_g).*(1-theta_minus)./dsigma(:,1:n2_elements) +...
            v_eff_plus.*theta_z(v_eff_plus,parameters) -...
            v_eff_minus.*(1-theta_z(v_eff_minus,parameters));

dFdT_zero(1,2:n2_elements) = ((u_plus(1,2:n2_elements) - sigma_plus(1,2:n2_elements).*x_g_deriv)./x_g).*theta_plus(1,2:n2_elements)./dsigma(1,2:n2_elements) -...
            (u(1)./x_g)./dsigma(1,2:n2_elements) +...
            v_eff_plus(1,2:n2_elements).*theta_z(v_eff_plus(1,2:n2_elements),parameters) -...
            v_eff_minus(1,2:n2_elements).*(1-theta_z(v_eff_minus(1,2:n2_elements),parameters));                        
     
dFdT_zero(1:n_elements-1,1) = ((u_plus(1:n_elements-1,1) - sigma_plus(1:n_elements-1,1).*x_g_deriv)./x_g).*theta_plus(1:n_elements-1,1)./dsigma(1:n_elements-1,1) -...
            ((u_minus(1:n_elements-1,1) - sigma_minus(1:n_elements-1,1).*x_g_deriv)./x_g).*(1-theta_minus(1:n_elements-1,1))./dsigma(1:n_elements-1,1) +...
            v_eff_plus(1:n_elements-1,1).*theta_z(v_eff_plus(1:n_elements-1,1),parameters) -...
            (((v(1:n_elements-1,1) - ((u(1:n_nodes-2,1)+u(2:n_nodes-1,1))./2).*dbdsigma(1:n_elements-1,1)./x_g)./h_elem(1:n_elements-1,1))./deta(1:n_elements-1,1));

dFdT_zero(n_elements,1:n2_elements) = ((u_plus(n_elements,1:n2_elements) - sigma_plus(n_elements,1:n2_elements).*x_g_deriv)./x_g)./dsigma(n_elements,1:n2_elements) -...
            ((u_minus(n_elements,1:n2_elements) - sigma_minus(n_elements,1:n2_elements).*x_g_deriv)./x_g).*(1-theta_minus(n_elements,1:n2_elements))./dsigma(n_elements,1:n2_elements) +...
            ((v_plus(n_elements,1:n2_elements) - eta_plus(n_elements,1:n2_elements).*(h_deriv_elem(n_elements,1:n2_elements) + (dhdsigma(n_elements,1:n2_elements)./x_g).*(u_elem(n_elements,1:n2_elements) - sigma_element(n_elements,1:n2_elements).*x_g_deriv)) - (u_elem(n_elements,1:n2_elements).*dbdsigma(n_elements,1:n2_elements)./x_g))./h_elem(n_elements,1:n2_elements)).*theta_z_plus(n_elements,1:n2_elements)./deta(n_elements,1:n2_elements) -...
            ((v_minus(n_elements,1:n2_elements) - eta_minus(n_elements,1:n2_elements).*(h_deriv_elem(n_elements,1:n2_elements) + (dhdsigma(n_elements,1:n2_elements)./x_g).*(u_elem(n_elements,1:n2_elements) - sigma_element(n_elements,1:n2_elements).*x_g_deriv)) - (u_elem(n_elements,1:n2_elements).*dbdsigma(n_elements,1:n2_elements)./x_g))./h_elem(n_elements,1:n2_elements)).*(1-theta_z_minus(n_elements,1:n2_elements))./deta(n_elements,1:n2_elements);                        

dFdT_zero(1,1) = ((u_plus(1,1) - sigma_plus(1,1).*x_g_deriv)./x_g).*theta_plus(1,1)./dsigma(1,1) -...
            (u(1)./x_g)./dsigma(1,1) +...
            v_eff_plus(1,1).*theta_z(v_eff_plus(1,1),parameters) -...
            (((v(1,1) - ((u(1,1)+u(2,1))./2).*dbdsigma(1,1)./x_g)./h_elem(1,1))./deta(1,1));
        
dFdT_zero(n_elements,1)   = ((u_plus(n_elements,1) - sigma_plus(n_elements,1).*x_g_deriv)./x_g)./dsigma(n_elements,1) -...
            ((u_minus(n_elements,1) - sigma_minus(n_elements,1).*x_g_deriv)./x_g).*(1-theta_minus(n_elements,1))./dsigma(n_elements,1) +...
            ((v_plus(n_elements,1) - eta_plus(n_elements,1).*(h_deriv_elem(n_elements,1) + (dhdsigma(n_elements,1)./x_g).*(u_elem(n_elements,1) - sigma_element(n_elements,1).*x_g_deriv)) - (u_elem(n_elements,1).*dbdsigma(n_elements,1)./x_g))./h_elem(n_elements,1)).*theta_z_plus(n_elements,1)./deta(n_elements,1) -...
            (((v(n_elements,1) - ((u(n_nodes-1,1)+u(n_nodes,1))./2).*dbdsigma(n_elements,1)./x_g)./h_elem(n_elements,1))./deta(n_elements,1));                        

%fder term
dFdT_zero = dFdT_zero + (1/dtau);

%stretch term 1
dFdT_zero = dFdT_zero + (x_g_deriv./x_g) +...
    (h_deriv_elem./h_elem) +...
    dhdsigma.*(u_elem - sigma_element.*x_g_deriv)./(h_elem.*x_g);

%diffusion ee term
dFdT_zero(:,2:n2_elements-1) = dFdT_zero(:,2:n2_elements-1) + 2.*parameters.Kappa.*...
    (((eta_element(:,2:n2_elements-1).^2)./((x_g.^2).*h_elem(:,2:n2_elements-1).^2)) +...
    (h_elem(:,2:n2_elements-1).^(-2)))./(deta(1:n_elements,2:n2_elements-1).^2);

dFdT_zero(:,1) = dFdT_zero(:,1) + parameters.Kappa.*...
    (((eta_element(:,1).^2)./((x_g.^2).*h_elem(:,1).^2)) +...
    (h_elem(:,1).^(-2)))./(deta(1:n_elements,1).^2);

dFdT_zero(:,n2_elements) = dFdT_zero(:,n2_elements) + 3.*parameters.Kappa.*...
    (((eta_element(:,n2_elements).^2)./((x_g.^2).*h_elem(:,n2_elements).^2)) +...
    (h_elem(:,n2_elements).^(-2)))./(deta(1:n_elements,n2_elements).^2);


%% (i-1,j) temperature term
dFdT_iminus = zeros(n_elements,n2_elements);

dFdT_iminus(:,1:n2_elements-1) = -v_eff_minus(:,2:n2_elements).*...
    theta_z(v_eff_minus(:,2:n2_elements),parameters);

%Diffusion ee term
dFdT_iminus(:,1:n2_elements-1) = dFdT_iminus(:,1:n2_elements-1) -...
    parameters.Kappa.*(((eta_element(:,2:n2_elements).^2)./((x_g.^2).*h_elem(:,2:n2_elements).^2)) +...
    (h_elem(:,2:n2_elements).^(-2)))./(deta(1:n_elements,2:n2_elements).^2);

%% (i+1,j) temperature term
dFdT_iplus = zeros(n_elements,n2_elements);

dFdT_iplus(:,2:n2_elements) = v_eff_plus(:,1:n2_elements-1).*...
    (1-theta_z(v_eff_plus(:,1:n2_elements-1),parameters));

%Diffusion ee term
dFdT_iplus(:,2:n2_elements) = dFdT_iplus(:,2:n2_elements) -...
    parameters.Kappa.*(((eta_element(:,1:n2_elements-1).^2)./((x_g.^2).*h_elem(:,1:n2_elements-1).^2)) +...
    (h_elem(:,1:n2_elements-1).^(-2)))./(deta(1:n_elements,1:n2_elements-1).^2);

%% assemble Jacobian

foutLHS = spdiags([dFdT_iminus(:) dFdT_minus(:) dFdT_zero(:) dFdT_plus(:) dFdT_iplus(:)], [-n_elements,-1:1,n_elements], n_elements*n2_elements, n_elements*n2_elements);

foutRHS = zeros(n_elements,n2_elements);
foutRHS(:,n2_elements) = -v_eff_plus(:,n2_elements).*(T_s.*(1-theta_z(v_eff_plus(:,n2_elements),parameters))) +...
    parameters.Kappa.*(((eta_element(:,n2_elements).^2)./((x_g.^2).*h_elem(:,n2_elements).^2)) +(h_elem(:,n2_elements).^(-2))).*(2.*parameters.T_s./(deta(1:n_elements,n2_elements).^2));

% foutRHS(:,1) = -parameters.Kappa.*(((eta_element(:,1).^2)./((x_g.^2).*h_elem(:,1).^2)) +(h_elem(:,1).^(-2))).*...
%     ((h.*bedHflux.*deta(1:n_elements,1))./(deta(1:n_elements,1).^2));

foutRHS = foutRHS + (T_old./parameters.dtau);

foutRHS=foutRHS(:);
end

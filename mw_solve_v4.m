function [e_new,h_till_new,HF,parameters] = mw_solve_v4(T,h,x_g,u,e,h_till,parameters)
%mw_solve uses an explicit Euler solver to solve the
%meltwater production budget and till water content
%evolution of an ice stream with specified temperature,
%velocity, thickness and grounding line position
%Written by Alex Robel, last modified December 19, 2013

%unpack grid
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
deta = eta_node(:,2:n2_nodes)-eta_node(:,1:n2_nodes-1); %spacing between eta nodes

u_elem = (u(1:n_nodes-1)+u(2:n_nodes))./2; %Put velocity onto elements

%% Calculate MW production
FrictionalHF = parameters.taub.*u_elem;
FrictionalHF = max([FrictionalHF';zeros(size(e))'])';

parameters.bedHflux = -(parameters.G + FrictionalHF);

HF = (parameters.k_i.*((T(:,2)-T(:,1))./...
    (h.*(parameters.grid.eta_element(2)-parameters.grid.eta_element(1)))) -...
    parameters.bedHflux); %Net Basal Heat Flux (in first eta layer)

frozen_zone = (e==parameters.e_min);
temperate_zone = (e>parameters.e_min);

becomes_frozen_zone = (e + temperate_zone.*parameters.dtau.*HF./...
    (parameters.rho.*parameters.L_f.*h_till))<parameters.e_min;
becomes_temperate_zone = (h_till + frozen_zone.*parameters.dtau.*HF./...
    (parameters.rho.*parameters.L_f))>parameters.h_till_max;

e_new = e + temperate_zone.*parameters.dtau.*HF./(parameters.rho.*parameters.L_f.*h_till); 
h_till_new = h_till + becomes_frozen_zone.*...
    parameters.h_till_max.*(e_new-parameters.e_min)./parameters.dtau;

% h_till_new = h_till + frozen_zone.*parameters.dtau.*HF./(parameters.rho.*parameters.L_f); 
h_till_new = h_till_new + frozen_zone.*parameters.dtau.*HF./(parameters.rho.*parameters.L_f); 
e_new = e_new + becomes_temperate_zone.*...
    (h_till_new-parameters.h_till_max)./(parameters.h_till_max.*parameters.dtau);

% e_new = min([e_new';parameters.e_max.*ones(size(e_new))'])';
e_new = max([e_new';parameters.e_min.*ones(size(e_new))'])';
h_till_new = min([h_till_new';parameters.h_till_max.*ones(size(h_till_new))'])';

end
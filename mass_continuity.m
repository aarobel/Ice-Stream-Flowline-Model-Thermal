function v = mass_continuity(h,x_g,u,parameters)
%mass_continuity finds the vertical velocity, v from
%horizontal velocity, u, using the mass continuity equation
%Written by Alex Robel, last modified August 30, 2013

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
deta = eta_node(:,2:n2_nodes)-eta_node(:,1:n2_nodes-1); %spacing between eta nodes                      %current w ice velocity

%calculate dhdt
% h_old = parameters.hx_g_old(1:end-1);
% dtau = parameters.dtau;
% h_deriv = (h-h_old)./dtau;
% h_deriv_eta_nodes = repmat(h_deriv,[1 parameters.grid.n2_nodes]);

sigma_elem_nodes = repmat(parameters.grid.sigma_element,[1 parameters.grid.n2_nodes]);
eta_elem_nodes = repmat(parameters.grid.eta_node,[1 parameters.grid.n_elements])';

dbdsigma                   = zeros(n_elements,n2_nodes);
dbdsigma(2:n_elements-1,:) = (Base(x_g.*sigma_elem_nodes(3:n_elements,:),parameters) - Base(x_g.*sigma_elem_nodes(1:n_elements-2,:),parameters))./(sigma_elem_nodes(3:n_elements,:)-sigma_elem_nodes(1:n_elements-2,:));
dbdsigma(1,:)              = (Base(x_g.*sigma_node(2,:),parameters) - Base(x_g.*sigma_node(1,:),parameters))./(sigma_node(2,:)-sigma_node(1,:));
dbdsigma(n_elements,:)     = (Base(x_g.*sigma_node(n_nodes,:),parameters) - Base(x_g.*sigma_node(n_nodes-1,:),parameters))./(sigma_node(n_nodes,:)-sigma_node(n_nodes-1,:));

urep = repmat((u(1:end-1,1)+u(2:end,1))./2,[1 parameters.grid.n2_nodes]);

%calculate dudx
dudx = diff(u,1,1)./(x_g.*dsigma);

%integrate to calculate w
h_eta_nodes = repmat(h,[1 parameters.grid.n2_nodes]);
% dz = deta(1:n_elements,:).*repmat(h,[1 parameters.grid.n2_elements]);
% v = h_deriv_eta_nodes - h_eta_nodes.*fliplr(cumtrapz(flipud(parameters.grid.eta_node),fliplr(dudx),2));

v = urep.*(dbdsigma./x_g) - h_eta_nodes.*eta_elem_nodes.*dudx;
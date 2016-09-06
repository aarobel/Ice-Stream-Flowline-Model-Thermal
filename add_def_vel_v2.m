function [u_old_full,u_old_mean] = add_def_vel_v2(u_old,h,x_g,parameters)

%add velocity deformation

h_nodes = zeros(parameters.grid.n_nodes,1);
h_nodes(2:end-1) = (h(1:end-1)+h(2:end))/2;
h_nodes(1) = 1.5*h(1) - 0.5*h(2);
h_nodes(end) = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
h_nodes_full = repmat(h_nodes,[1 parameters.grid.n2_nodes]);

b = Base(x_g.*parameters.grid.sigma_element,parameters);
dsdx_nodes = zeros(parameters.grid.n_nodes,1);
dsdx_nodes(2:end-1) = diff(h+b)./(x_g.*diff(parameters.grid.sigma_element));
dsdx_nodes(end) = dsdx_nodes(end-1);
taud_nodes = -parameters.rho.*parameters.g.*h_nodes.*dsdx_nodes;
taud_nodes_full = repmat(taud_nodes,[1 parameters.grid.n2_nodes]);

B_Glen_nodes_full = zeros(parameters.grid.n_nodes,parameters.grid.n2_nodes);
B_Glen_nodes_full(2:parameters.grid.n_nodes-1,1:parameters.grid.n2_nodes-1) = (parameters.B_Glen_full(1:parameters.grid.n_elements-1,:)+parameters.B_Glen_full(2:parameters.grid.n_elements,:))./2;
B_Glen_nodes_full(1,1:parameters.grid.n2_nodes-1) = 1.5*parameters.B_Glen_full(1,:) - 0.5*parameters.B_Glen_full(2,:);
B_Glen_nodes_full(parameters.grid.n_nodes,1:parameters.grid.n2_nodes-1) = 1.5*parameters.B_Glen_full(parameters.grid.n_nodes-1,:) - 0.5*parameters.B_Glen_full(parameters.grid.n_nodes-2,:);

B_Glen_nodes_full(:,2:parameters.grid.n2_nodes-1) = (B_Glen_nodes_full(:,1:parameters.grid.n2_elements-1)+B_Glen_nodes_full(:,2:parameters.grid.n2_elements))./2;
B_Glen_nodes_full(:,1) = 1.5*B_Glen_nodes_full(:,2) - 0.5*B_Glen_nodes_full(:,3);
B_Glen_nodes_full(:,parameters.grid.n2_nodes) = 1.5*B_Glen_nodes_full(:,parameters.grid.n2_nodes-1) - 0.5*B_Glen_nodes_full(:,parameters.grid.n2_nodes-2);

B_Glen_nodes_mean = mean(B_Glen_nodes_full,2);
B_Glen_nodes_full = repmat(B_Glen_nodes_mean,[1 parameters.grid.n2_nodes]);

eta_full = repmat(parameters.grid.eta_node',[parameters.grid.n_nodes,1]);

%calculate mean column velocity with deformation
u_old_mean = u_old + ((2.*B_Glen_nodes_mean.^(-parameters.n_Glen)).*...
    h_nodes.*(taud_nodes.^parameters.n_Glen))/(parameters.n_Glen+2);
u_old_mean(1) = 0;

%calculate full column velocity
u_old_full = repmat(u_old,[1 parameters.grid.n2_nodes]) + ((2.*B_Glen_nodes_full.^(-parameters.n_Glen)).*...
    h_nodes_full.*(1-(1-eta_full).^(parameters.n_Glen+1)).*(taud_nodes_full.^parameters.n_Glen))/(parameters.n_Glen+1);
u_old_full(1,:) = 0;
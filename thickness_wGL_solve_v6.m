function [h_out,x_g_out,iter] = thickness_wGL_wdiff_solve(h,x_g,u,parameters,flags)
%thickness_wGL_solve uses a Newton's method minimization
%algorithm to find the flowline thickness and 
%grounding line positions of an ice stream
%with specified velocity
%Written by Alex Robel, last modified Oct. 10

if nargin < 5 || ~isfield(flags,'test')
    flags.test = false;
end

if flags.test
    if parameters.grid.n_nodes > 500, error('Run derivative test only with small number of degrees of freedom (< 40)'), end
    parameters.u=u;
    h_out = Jacobian_test_v2(@thickness_function,@thickness_Jacobian,[h;x_g],parameters,sqrt(eps));
    x_g_out=[];
    iter=[];
    return
end
       
%tolerances for Newton solver
srchparams.toldelta = 1e5*(parameters.grid.n_nodes)*sqrt(eps);
srchparams.tolgrad = 50*(parameters.grid.n_nodes)*sqrt(eps);%/min(parameters.u_eps,parameters.D_eps);
srchparams.tolF = 50*(parameters.grid.n_nodes)*sqrt(eps);
srchparams.verbose = 0;
srchparams.itmax=parameters.hiter_max;

%solve
hx_g = [h;x_g]; %initial guess
parameters.u=u;

[hx_g_out,error_flag,iter] =  Newton(@thickness_function,@thickness_Jacobian,hx_g,parameters,srchparams);
h_out = hx_g_out(1:end-1);
x_g_out = hx_g_out(end);
end

function fout = thickness_function(hx_g,parameters)
%note that the input to this function hx_g is the
%ice thickness column vector concatenated with the
%grounding line position in the last row

%unpack variables
h = hx_g(1:end-1);
x_g = hx_g(end);

h_old = parameters.hx_g_old(1:end-1);
x_g_old = parameters.hx_g_old(end);

%unpack parameters
rho = parameters.rho; %ice density
rho_w = parameters.rho_w; %water density

%unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1), n_nodes-by-one column vector
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1), n_elements-by-one column vector
dsigma = sigma_node(2:n_nodes)-sigma_node(1:n_nodes-1); %spacing between nodes
dsigma_elem = sigma_element(2:n_elements)-sigma_element(1:n_elements-1);

%unpack other things from parameters structure
u = parameters.u;                       %current ice velocity
width_nodes = width_fun(x_g.*sigma_node,parameters);
width_elems = width_fun(x_g.*sigma_element,parameters);
width_Nplus_ghost = 2*width_elems(n_elements)-width_elems(n_elements-1);
width_Nminus_ghost = 2*width_elems(1)-width_elems(2);

%unpack time step
dtau = parameters.dtau;

%indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_n_int = (2:n_nodes-1)';
ind_e_minus = (1:n_elements-1)';
ind_e_plus = (2:n_elements)';

%Compute derivative term
Fder = dsigma.*width_elems.*(h-h_old)./dtau;

%Compute advection terms

%calculate nodal interior terms
u_int = u(ind_n_int);
sigma_int = sigma_node(ind_n_int);
theta_int = theta(x_g,x_g_old,u_int,sigma_int,parameters);

x_g_deriv = (x_g - x_g_old)/dtau;
h_g = -(rho_w/rho).*Base(x_g,parameters);

% u_eff_int = (u_int - sigma_int.*x_g_deriv)./(x_g.*width_nodes(ind_n_int));
u_eff_int = (u_int - sigma_int.*x_g_deriv)./x_g;

% u_eff_GL = ((u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./(x_g.*width_nodes(n_nodes)));
u_eff_GL = (u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./x_g;
theta_GL = theta(x_g,x_g_old,u(n_nodes),sigma_node(n_nodes),parameters);
h_Nplus_ghost = 2*h(n_elements)-h(n_elements-1);
h_Nminus_ghost = 2*h(1)-h(2);

%Plus advection term
Fadv_plus = zeros(n_elements,1); %specify last term separately
Fadv_plus(1:n_elements-1) = u_eff_int.*(width_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int) + width_elems(ind_e_minus).*h(ind_e_minus).*theta_int);
% Fadv_plus(1:n_elements-1) = u_eff_int.*(h(ind_e_plus).*(1-theta_int) + h(ind_e_minus).*theta_int);

%Last term can be spcified differently because we know thickness at the GL
% Fadv_plus(n_elements) = ((u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./x_g).*h_g;
Fadv_plus(n_elements) = u_eff_GL.*(width_Nplus_ghost.*h_Nplus_ghost.*(1-theta_GL) + width_elems(n_elements).*h(n_elements).*theta_GL);% +...

%Minus advection term
Fadv_minus = zeros(n_elements,1); %ensure that first term is zero
Fadv_minus(2:n_elements) = -u_eff_int.*(width_elems(ind_e_minus).*h(ind_e_minus).*theta_int + width_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int));% - ...
% Fadv_minus(2:n_elements) = -u_eff_int.*(h(ind_e_minus).*theta_int + h(ind_e_plus).*(1-theta_int));

Fadv_minus(1) = -(u(1)./x_g).*h(1).*width_elems(1);

%Compute stretching term
Fstretch = dsigma.*h_old.*width_elems.*(1-(x_g_old./x_g))./dtau;

%Compute source term
Fsource = -accumulation(parameters).*dsigma.*width_elems;

%assemble thickness minimization
fout = zeros(n_nodes,1);
fout(1:n_elements) = Fder+Fadv_plus+Fadv_minus+Fstretch+Fsource;
% fout(1:n_elements) = Fadv_plus+Fadv_minus;

%assemble GL position minimization
fout(n_nodes) = 1.5.*h(n_elements) - 0.5*h(n_elements-1) - h_g;
% fout(n_nodes) = h_g - h(n_elements);

end

function fout = thickness_Jacobian(hx_g,parameters)
%note that the input to this function hx_g is the
%ice thickness column vector concatenated with the
%grounding line position in the last row

%unpack variables
h = hx_g(1:end-1);
x_g = hx_g(end);

h_old = parameters.hx_g_old(1:end-1);
x_g_old = parameters.hx_g_old(end);

%unpack parameters
rho = parameters.rho; %ice density
rho_w = parameters.rho_w; %water density

%unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1), n_nodes-by-one column vector
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1), n_elements-by-one column vector
dsigma = sigma_node(2:n_nodes)-sigma_node(1:n_nodes-1); %spacing between nodes
dsigma_elem = sigma_element(2:n_elements)-sigma_element(1:n_elements-1);

%unpack other things from parameters structure
u = parameters.u;                       %current ice velocity
width_nodes = width_fun(x_g.*sigma_node,parameters);
width_elems = width_fun(x_g.*sigma_element,parameters);
width_Nplus_ghost = 2*width_elems(n_elements)-width_elems(n_elements-1);
width_Nminus_ghost = 2*width_elems(1)-width_elems(2);

dwdx_nodes = dwidthdx(x_g.*sigma_node,parameters);
dwdx_elems = dwidthdx(x_g.*sigma_element,parameters);

%unpack time step
dtau = parameters.dtau;

%indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_n_int = (2:n_nodes-1)';
ind_e_minus = (1:n_elements-1)';
ind_e_plus = (2:n_elements)';

%Set Dirichlet condition
% uD = u(1);
% u(1) = parameters.u_in;

%calculate thickness dertivatives

%calculate nodal interior terms
u_int = u(ind_n_int);
sigma_int = sigma_node(ind_n_int);
theta_int = theta(x_g,x_g_old,u_int,sigma_int,parameters);
dthetadx_g_int = dthetadx_g(x_g,x_g_old,u_int,sigma_int,parameters);

h_plus = h(ind_e_plus);
u_plus = u(ind_n_plus);
sigma_plus = sigma_node(ind_n_plus);
theta_plus = theta(x_g,x_g_old,u_plus,sigma_plus,parameters);
dthetadx_g_plus = dthetadx_g(x_g,x_g_old,u_plus,sigma_plus,parameters);
width_n_plus = width_nodes(ind_n_plus);
width_e_plus = width_elems(ind_e_plus);

dwdx_n_plus = dwdx_nodes(ind_n_plus);
dwdx_e_plus = dwdx_elems(ind_e_plus);

h_minus = h(ind_e_minus);
u_minus = u(ind_n_minus);
sigma_minus = sigma_node(ind_n_minus);
theta_minus = theta(x_g,x_g_old,u_minus,sigma_minus,parameters);
dthetadx_g_minus = dthetadx_g(x_g,x_g_old,u_minus,sigma_minus,parameters);
width_n_minus = width_nodes(ind_n_minus);
width_e_minus = width_elems(ind_e_minus);

dwdx_n_minus = dwdx_nodes(ind_n_minus);
dwdx_e_minus = dwdx_elems(ind_e_minus);

x_g_deriv = (x_g - x_g_old)/dtau;
h_g = -(rho_w/rho).*Base(x_g,parameters);
dh_gdx_g = -(rho_w/rho).*dBasedx(x_g,parameters);

u_eff_GL = ((u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./x_g);
theta_GL = theta(x_g,x_g_old,u(n_nodes),sigma_node(n_nodes),parameters);
dtheta_dx_g_GL = dthetadx_g(x_g,x_g_old,u(n_nodes),sigma_node(n_nodes),parameters);
h_Nplus_ghost = 2*h(n_elements)-h(n_elements-1);
h_Nminus_ghost = 2*h(1)-h(2);

% theta_int = ones(size(theta_int));
% theta_plus = ones(size(theta_plus));
% theta_minus = ones(size(theta_minus));
% theta_GL=1;
% dthetadx_g_plus = zeros(size(dthetadx_g_plus));
% dthetadx_g_minus = zeros(size(dthetadx_g_minus));

%j+1 thickness term (on thickness)
dFdh_plus = zeros(n_elements,1); %ensure that last term is zero
dFdh_plus(2:n_elements) = ((u_int - sigma_int.*x_g_deriv)./(x_g)).*(1-theta_int).*width_elems(ind_e_plus);

%j-1 thickness term (on thickness)
dFdh_minus = zeros(n_elements,1); %ensure that first term is zero
dFdh_minus(1:n_elements-1) = -((u_int - sigma_int.*x_g_deriv)./(x_g)).*theta_int.*width_elems(ind_e_minus);

%j thickness term (on thickness)
dFdh_zero = zeros(n_elements,1);
dFdh_zero(2:n_elements-1) = (width_elems(2:n_elements-1).*dsigma(2:n_elements-1)/dtau) + ((u_plus(2:n_elements-1) - sigma_plus(2:n_elements-1).*x_g_deriv)./(x_g)).*theta_plus(2:n_elements-1).*width_elems(2:n_elements-1) -...
                            ((u_minus(2:n_elements-1) - sigma_minus(2:n_elements-1).*x_g_deriv)./(x_g)).*(1-theta_minus(2:n_elements-1)).*width_elems(2:n_elements-1);

dFdh_zero(1) = (width_elems(1).*dsigma(1)/dtau) + ((u_plus(1) - sigma_plus(1).*x_g_deriv)./(x_g)).*theta_plus(1).*width_elems(1) -...
                            ((u_minus(1))./(x_g)).*width_elems(1);

% dFdh_zero(1) = (dsigma(1)/dtau) + ((u_plus(1) - sigma_plus(1).*x_g_deriv)./x_g).*theta_plus(1) +...
%                             dsigma(1).*(1-(x_g_old./x_g))./dtau +...
%                             parameters.HKappa.*(1./(x_g^2))./dsigma_elem(1);


dFdh_zero(n_elements) = (width_elems(n_elements).*dsigma(n_elements)/dtau) + ((u_plus(n_elements) - x_g_deriv)./(x_g)).*width_elems(n_elements) -...
                            ((u_minus(n_elements) - sigma_minus(n_elements).*x_g_deriv)./(x_g)).*(1-theta_minus(n_elements)).*width_elems(n_elements);

%x_g term (on thickness)
dFdx_g = zeros(n_elements,1);

%Fadv_plus term
dFdx_g(1:n_elements-1) = -(1/(x_g.^2)).*(u_int + (sigma_int.*x_g_old./dtau)).*...
    (width_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int) + width_elems(ind_e_minus).*h(ind_e_minus).*theta_int) +...
                          ((u_int - sigma_int.*x_g_deriv)./x_g).*...
    (dwdx_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int) + dwdx_elems(ind_e_minus).*h(ind_e_minus).*theta_int) +...
                          ((u_int - sigma_int.*x_g_deriv)./x_g).*...
    (width_elems(ind_e_plus).*h(ind_e_plus).*(-dthetadx_g_int) + width_elems(ind_e_minus).*h(ind_e_minus).*dthetadx_g_int);

dFdx_g(n_elements) = -(1/(x_g.^2)).*(u(n_nodes) + (sigma_node(n_nodes).*x_g_old./dtau)).*...
    (width_Nplus_ghost.*h_Nplus_ghost.*(1-theta_GL) + width_elems(n_elements).*h(n_elements).*theta_GL) +...
                          ((u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./x_g).*...
    ((2*dwdx_elems(n_elements)-dwdx_elems(n_elements-1)).*h_Nplus_ghost.*(1-theta_GL) + dwdx_elems(n_elements).*h(n_elements).*theta_GL) +...
                          ((u(n_nodes) - sigma_node(n_nodes).*x_g_deriv)./x_g).*...
    (width_Nplus_ghost.*h_Nplus_ghost.*(-dtheta_dx_g_GL) + width_elems(n_elements).*h(n_elements).*dtheta_dx_g_GL);

%Fadv_minus term
dFdx_g(2:n_elements) = dFdx_g(2:n_elements) + (1/(x_g.^2)).*(u_int + (sigma_int.*x_g_old./dtau)).*...
    (width_elems(ind_e_minus).*h(ind_e_minus).*theta_int + width_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int)) -...
                          ((u_int - sigma_int.*x_g_deriv)./x_g).*...
    (dwdx_elems(ind_e_minus).*h(ind_e_minus).*theta_int + dwdx_elems(ind_e_plus).*h(ind_e_plus).*(1-theta_int)) -...
                          ((u_int - sigma_int.*x_g_deriv)./x_g).*...
    (width_elems(ind_e_minus).*h(ind_e_minus).*dthetadx_g_int + width_elems(ind_e_plus).*h(ind_e_plus).*(-dthetadx_g_int));

%Stretch and source terms
dFdx_g(1:n_elements) = dFdx_g(1:n_elements) + dsigma.*h_old.*width_elems.*(x_g_old./(x_g^2))./dtau +...
    dsigma.*h_old.*dwdx_elems.*(1-(x_g_old./x_g))./dtau - accumulation(parameters).*dsigma.*dwdx_elems;

%upstream term
% dFdx_g(1) =              -(1/(x_g.^2)).*...
%     ((u_plus(1) + (sigma_plus(1).*x_g_old./dtau))).*...
%     (width_e_plus(1).*h_plus(1).*(1-theta_plus(1)) + width_elems(1).*h(1).*theta_plus(1)) +...
%                          (h(1).*x_g_old.*dsigma(1)./(dtau.*(x_g.^2)));

% dFdx_g(1) =              -(1/(x_g.^2)).*...
%     (u_plus(1) + (sigma_plus(1).*x_g_old./dtau)).*...
%     (h_plus(1).*(1-theta_plus(1)) + h(1).*theta_plus(1)) +...
%                          ((u_plus(1) - sigma_plus(1).*x_g_deriv)./x_g).*...
%     (-h_plus(1).*dthetadx_g_plus(1) + h(1).*dthetadx_g_plus(1)) +...
%                          (h(1).*x_g_old.*dsigma(1)./(dtau.*(x_g.^2)));
                      
%downstream term
% dFdx_g(n_elements) =    -(1/(x_g.^2)).*((u(n_nodes) + (x_g_old./dtau))./width_nodes(n_nodes)).*h(n_elements).*width_elems(n_elements) +...
%                          (1/(x_g.^2)).*...
%     ((u_minus(n_elements) + (sigma_minus(n_elements).*x_g_old./dtau))./width_n_minus(n_elements)).*...
%     (width_e_minus(n_elements-1).*h_minus(n_elements-1).*theta_minus(n_elements) + width_elems(n_elements).*h(n_elements).*(1-theta_minus(n_elements))) -...
%                          ((u_minus(n_elements) - sigma_minus(n_elements).*x_g_deriv)./(x_g.*width_n_minus(n_elements))).*...
%     (width_e_minus(n_elements-1).*h_minus(n_elements-1).*dthetadx_g_minus(n_elements) - width_elems(n_elements).*h(n_elements).*dthetadx_g_minus(n_elements)) +...
%                          (h(n_elements).*x_g_old.*dsigma(n_elements)./(dtau.*(x_g.^2)));

%boundary condition derivative
dGdh_minus = -0.5;
dGdh_zero = 1.5;
dGdx_g = -dh_gdx_g;

% dGdh_minus = 0;
% dGdh_zero = -1;
% dGdx_g = dh_gdx_g;

%assemble Jacobian
dF = spdiags([dFdh_minus dFdh_zero dFdh_plus], -1:1, n_elements, n_elements);

fout = spalloc(n_nodes,n_nodes,4*n_nodes+3);
fout(1:n_elements,1:n_elements) = dF;
fout(1:n_elements,n_nodes) = dFdx_g;
fout(n_nodes,n_elements-1) = dGdh_minus;
fout(n_nodes,n_elements) = dGdh_zero;
fout(n_nodes,n_nodes) = dGdx_g;

end

function fout = accumulation(parameters)
%accumulation on elements
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements

fout = parameters.accumrate.*ones(n_elements,1);
end

function fout = theta(x_g,x_g_old,u,sigma_node,parameters)
%Heavyside for ensuring good behavior of grounding line

HS_sensitivity = parameters.HS_sensitivity;

%unpack time step
dtau = parameters.dtau;
%grounding line migration rate
x_g_deriv = (x_g - x_g_old)/dtau;

heavy_var = u - sigma_node.*x_g_deriv;

fout = 0.5.*(1+tanh(HS_sensitivity.*heavy_var)); % here we use a hyperbolic tan function as a regularized hevyside function
% fout=ones(size(u));
end

function fout = dthetadx_g(x_g,x_g_old,u,sigma_node,parameters)
%Derivative of regularized Heavyside function for use in Jacobian
HS_sensitivity = parameters.HS_sensitivity;

%unpack time step
dtau = parameters.dtau;
%grounding line migration rate
x_g_deriv = (x_g - x_g_old)/dtau;

heavy_var = u - sigma_node.*x_g_deriv;

fout = -HS_sensitivity.*(sigma_node./dtau).*(0.5.*sech(HS_sensitivity.*heavy_var).^2);
% fout = zeros(size(u));
end

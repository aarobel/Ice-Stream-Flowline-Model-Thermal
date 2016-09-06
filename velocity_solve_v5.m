function [uout,error_flag,iter] = velocity_solve_v4(u,h,x_g,parameters,flags)
%velocity_solve uses a Newton's method minimization
%algorithm to find the flowline velocity of an ice stream
%with specified thickness, tau_b and tau_l
%Written by Christian/Alex/James
%Also need to set default parameters values, e.g.
% if nargin = 0
%   u = blag
%end
%if nargin <2  || ~isfield(parameters,'B_Glen')
%   B_Glen = blah
%end

%pack thickness and grounding line position
parameters.h=h;
parameters.x_g=x_g;

if nargin < 5 || ~isfield(flags,'test')
    flags.test = false;
end

if flags.test
    if parameters.grid.n_nodes > 40, error('Run derivative test only with small number of degrees of freedom (< 40)'), end
    uout = difftest_v3(@stream_Lagrangian,@stream_gradient,@stream_Hessian,rand(size(u)),parameters,sqrt(eps));
    return
end

%tolerances for Newton solver
srchparams.toldelta = (parameters.grid.n_nodes)*sqrt(eps);
srchparams.tolgrad = 100.*max(parameters.B_Glen)*(parameters.grid.n_nodes)*sqrt(eps);
srchparams.tolF = max(parameters.B_Glen)*(parameters.grid.n_nodes)*sqrt(eps);
srchparams.verbose = parameters.uverbose;
srchparams.itmax=parameters.uiter_max;

%solve
% [uout,error_flag] =  Newton_min(@stream_Lagrangian,@stream_gradient,@stream_Hessian,u,parameters,srchparams);
[uout,error_flag,iter] =  Newton_min(@stream_Lagrangian,@stream_gradient,@stream_Hessian,u,parameters,srchparams);
%set Dirichlet condition
%uout(1)=parameters.u_in;

end

function Jout = stream_Lagrangian(u,parameters)

%unpack parameters
B_Glen = parameters.B_Glen; %B in Glen's law (vertically averaged if necessary)
n_Glen = parameters.n_Glen; %n in Glen's law
rho = parameters.rho; %ice density
rho_w = parameters.rho_w; %water density
g = parameters.g;   %acceleration due to gravity
D_eps = parameters.D_eps; %strain rate regularizer
u_eps = parameters.u_eps; %velocity regualrizer

%unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1), n_nodes-by-one column vector
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1), n_elements-by-one column vector
dsigma = sigma_node(2:n_nodes)-sigma_node(1:n_nodes-1); %spacing between nodes
%dsigma = parameters.grid.dsigma %n_elements-by-one list of element lengths

%unpack other dynamical variables from parameters structure
h = parameters.h;                       %current ice thickness
x_g = parameters.x_g;                   %current grounding line position

%compute edge lengths and grid point positions
dx = x_g*dsigma;
x_node = x_g*sigma_node;
x_element = x_g*sigma_element;

%indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_e_plus = [2:n_elements n_elements]';
ind_e_minus = [1 1:n_elements-1]';

%Set Dirichlet condition
uD = u(1);
u(1) = parameters.u_in;

%Compute viscous term
dudx = (u(ind_n_plus)-u(ind_n_minus))./dx;
Jvisc = (n_Glen/(n_Glen+1))*sum(2*B_Glen.*h.*sqrt(dudx.^2+D_eps^2).^(1/n_Glen+1).*dx);

%Compute lateral shear term
Jlateral_integrand = n_Glen/(n_Glen+1).*DupontG(x_element,parameters).*h.*(sqrt(u(ind_n_minus).^2+u_eps.^2).^(1/n_Glen+1)+sqrt(u(ind_n_plus).^2+u_eps^2).^(1/n_Glen+1)).*dx/2;
Jlateral = sum(Jlateral_integrand);

%Compute basal shear stress term
Jbase_integrand = (IntTauB(u(ind_n_minus),parameters)+IntTauB(u(ind_n_plus),parameters)).*dx/2;
Jbase = sum(Jbase_integrand);

%Compute driving term
b = Base(x_element,parameters);
slope_plus = ((h+b)-(h(ind_e_plus)+b(ind_e_plus)))./dx;
slope_plus(n_elements) = slope_plus(n_elements-1);
slope_minus =  ((h(ind_e_minus)+b(ind_e_minus))-(h+b))./dx;
Jdrive_integrand = zeros(n_nodes,1);
Jdrive_integrand(ind_e_plus) = Jdrive_integrand(ind_e_plus) + rho.*g.*h.*slope_plus.*u(ind_e_plus).*dx/2;
Jdrive_integrand(ind_e_minus) = Jdrive_integrand(ind_e_minus) + rho.*g.*h.*slope_minus.*u(ind_e_minus).*dx/2;
Jdrive_integrand(n_nodes) = 0;
Jdrive = sum(Jdrive_integrand);

%Boundary term
if ~parameters.float
    h_g = parameters.h_g;
    T_f = 0.5*rho*g*h_g^2;   
else
    h_g = -(rho_w/rho)*Base(x_g,parameters);
    T_f = (1-parameters.buttress).*0.5*rho*g*(1-rho/rho_w).*(h_g^2);
end
Jbdy = T_f.*u(n_nodes);

%assemble
Jout = Jvisc+Jlateral+Jbase-Jdrive-Jbdy+uD^2/2;

end

function fout = stream_gradient(u,parameters)

%unpack parameters
B_Glen = parameters.B_Glen; %B in Glen's law (vertically averaged if necessary)
n_Glen = parameters.n_Glen; %n in Glen's law
rho = parameters.rho; %ice density
rho_w = parameters.rho_w; %water density
g = parameters.g; %acceleration due to gravity
u_eps = parameters.u_eps; %velocity regularizer
D_eps = parameters.D_eps; %strain rate regularizer

%unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1)
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1)
dsigma = sigma_node(2:n_nodes)-sigma_node(1:n_nodes-1); %spacing between nodes
%dsigma = parameters.grid.dsigma %n_elements-by-one list of element lengths

%unpack other dynamical variables from parameters structure
h = parameters.h;                       %current ice thickness
x_g = parameters.x_g;                   %current grounding line position

%compute edge lengths and grid point positions
dx = x_g*dsigma;
x_node = x_g*sigma_node;
x_element = x_g*sigma_element;

%indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_e_plus = [2:n_elements n_elements]';
ind_e_minus = [1 1:n_elements-1]';

%Dirichlet condition
uD = u(1);
u(1) = parameters.u_in;

%Differentiate Jvisc term: Jvisc = n_Glen/(n_Glen+1)*sum(h.*sqrt(dudx.^2+D_eps^2).^(1/n_Glen+1).*dx);
dudx = (u(ind_n_plus)-u(ind_n_minus))./dx;
dJviscdD = 2*B_Glen.*h.*((dudx.^2+D_eps^2).^((1-n_Glen)/(2*n_Glen))).*dudx.*dx;
dDdu_minus = -1./dx;
dDdu_plus = 1./dx;
fvisc = zeros(n_nodes,1);
fvisc(ind_n_minus) = fvisc(ind_n_minus) + dJviscdD.*dDdu_minus;
fvisc(ind_n_plus) = fvisc(ind_n_plus) + dJviscdD.*dDdu_plus;

%Differentiate Jlateral term: Jlateral = n_Glen/(n_Glen+1)*sum(G.*h.*(sqrt(u(ind_n_minus).^2+u_eps^2).^(1/n_Glen+1)+sqrt(u(ind_n_plus).^2+u_eps^2).^(1/n_Glen+1)).*dx/2);
flateral = zeros(n_nodes,1);
flateral(ind_n_minus) = flateral(ind_n_minus) + DupontG(x_element,parameters).*h.*(u(ind_n_minus).^2+u_eps^2).^((1-n_Glen)/(2*n_Glen)).*u(ind_n_minus).*dx/2;
flateral(ind_n_plus) =  flateral(ind_n_plus) + DupontG(x_element,parameters).*h.*(u(ind_n_plus).^2+u_eps^2).^((1-n_Glen)/(2*n_Glen)).*u(ind_n_plus).*dx/2;

%Differentiate Jbase term: Jbase = sum((IntTauB(u(ind_n_minus),parameters)+IntTauB(u(ind_n_minus),parameters)).*dx)
fbase = zeros(n_nodes,1);
fbase(ind_n_minus) = fbase(ind_n_minus) + TauB(u(ind_n_minus),parameters).*dx/2;
fbase(ind_n_plus) = fbase(ind_n_plus) + TauB(u(ind_n_plus),parameters).*dx/2;

%Differentiate driving term
b = Base(x_element,parameters);
slope_plus = ((h+b)-(h(ind_e_plus)+b(ind_e_plus)))./dx;
slope_plus(n_elements) = slope_plus(n_elements-1);
slope_minus =  ((h(ind_e_minus)+b(ind_e_minus))-(h+b))./dx;
fdrive = zeros(n_nodes,1);
fdrive(ind_e_plus) = fdrive(ind_e_plus) + rho.*g.*h.*slope_plus.*dx/2;
fdrive(ind_e_minus) = fdrive(ind_e_minus) + rho.*g.*h.*slope_minus.*dx/2;
fdrive(n_nodes) = 0;

%Differentiate boundary term
fbdy = zeros(n_nodes,1);
if ~parameters.float
    h_g = parameters.h_g;
    T_f = 0.5*rho*g*h_g^2;
else
    h_g = -(rho_w/rho)*Base(x_g,parameters);
    T_f = (1-parameters.buttress).*0.5*rho*g*(1-rho/rho_w).*(h_g^2);
end

fbdy(n_nodes) = T_f;

%Assemble
fout = fvisc+flateral+fbase-fdrive-fbdy;

%Deal with Dirichlet
fout(1) = uD;
end

function Hout = stream_Hessian(u,parameters)

%unpack parameters
B_Glen = parameters.B_Glen; %B in Glen's law (vertically averaged if necessary)
n_Glen = parameters.n_Glen; %n in Glen's law
rho = parameters.rho; %ice density
rho_w = parameters.rho_w; %water density
D_eps = parameters.D_eps; %strain rate regularizer
u_eps = parameters.u_eps; %velocity regularizer

%unpack grid
n_nodes = parameters.grid.n_nodes;      %number of node positions, u has length n_nodes
n_elements = parameters.grid.n_elements;  %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
sigma_node = parameters.grid.sigma_node;     %node positions scaled to (0,1)
sigma_element = parameters.grid.sigma_element;  %element centres scaled to (0,1)
dsigma = sigma_node(2:n_nodes)-sigma_node(1:n_nodes-1); %spacing between nodes
%dsigma = parameters.grid.dsigma %n_elements-by-one list of element lengths

%unpack other dynamical variables from parameters structure
h = parameters.h;                       %current ice thickness
x_g = parameters.x_g;                   %current grounding line position

%compute edge lengths and grid point positions
dx = x_g*dsigma;
x_node = x_g*sigma_node;
x_element = x_g*sigma_element;

%indexing
ind_n_minus = (1:n_nodes-1)';
ind_n_plus = (2:n_nodes)';
ind_e_plus = [2:n_elements n_elements]';
ind_e_minus = [1 1:n_elements-1]';

%set Dirichlet condition
u(1) = parameters.u_in;

%Differentiate fvisc term: fvisc_minus = dJviscdD.*dDdu_minus; 
%fvisc_plus = dJviscdD.*dDdu_plus;
dudx = (u(ind_n_plus)-u(ind_n_minus))./dx;
dJviscdD2 = 2*B_Glen.*h.*(dudx.^2+D_eps^2).^((1-n_Glen)/(2*n_Glen)).*dx +...
    ((1-n_Glen)/n_Glen)*2*B_Glen.*h.*(dudx.^2+D_eps^2).^((1-3*n_Glen)/(2*n_Glen)).*(dudx.^2).*dx;
dDdu_minus = -1./dx;
dDdu_plus = 1./dx;
Hvisc_minusminus = dJviscdD2.*dDdu_minus.^2;
Hvisc_minusplus = dJviscdD2.*dDdu_plus.*dDdu_minus;
Hvisc_plusminus = dJviscdD2.*dDdu_minus.*dDdu_plus;
Hvisc_plusplus = dJviscdD2.*dDdu_plus.^2;
Hvisc = sparse([ind_n_minus; ind_n_minus; ind_n_plus; ind_n_plus],[ind_n_minus; ind_n_plus; ind_n_minus; ind_n_plus],[Hvisc_minusminus; Hvisc_minusplus; Hvisc_plusminus; Hvisc_plusplus],n_nodes,n_nodes);

%Differentiate flateral term
Hlateral_minusminus = DupontG(x_element,parameters).*h.*((1-n_Glen)/n_Glen*(u(ind_n_minus).^2+u_eps^2).^((1-3*n_Glen)/(2*n_Glen)).*u(ind_n_minus).^2  +  (u(ind_n_minus).^2+u_eps^2).^((1-n_Glen)/(2*n_Glen)) ).*dx/2;
Hlateral_plusplus = DupontG(x_element,parameters).*h.*((1-n_Glen)/n_Glen*(u(ind_n_plus).^2+u_eps^2).^((1-3*n_Glen)/(2*n_Glen)).*u(ind_n_plus).^2  +  (u(ind_n_plus).^2+u_eps^2).^((1-n_Glen)/(2*n_Glen)) ).*dx/2;
Hlateral = sparse([ind_n_minus; ind_n_plus],[ind_n_minus; ind_n_plus],[Hlateral_minusminus; Hlateral_plusplus],n_nodes,n_nodes);
Hlateral(n_nodes,:) = 0;

%Differentiate fbase term
Hbase_minusminus = DTauB(u(ind_n_minus),parameters).*dx/2;
Hbase_plusplus = DTauB(u(ind_n_plus),parameters).*dx/2;
Hbase = sparse([ind_n_minus; ind_n_plus],[ind_n_minus; ind_n_plus],[Hbase_minusminus; Hbase_plusplus],n_nodes,n_nodes);

%No need to differentiate driving or boundary terms as they are linear
Hout = Hvisc + Hlateral + Hbase;
Hout(:,1)=0;
Hout(1,:)=0;
Hout(1,1)=1;
end

function fout = IntTauB(u,parameters)
%antiderivative of TauB with respect to u

switch parameters.frictionlaw
  case 'Coulomb'
    %unpack parameters
%     taub = (parameters.taub(1:end-1)+parameters.taub(2:end))./2;
    taub = parameters.taub;
    u_eps = parameters.u_eps;   %regularizing velocity
    fout = taub.*sqrt(u.^2+u_eps.^2);
  case 'Weertman'
    %unpack
    C_schoof = parameters.C_schoof;
    m_schoof = parameters.m_schoof;
    u_eps = parameters.u_eps;
    fout = 1/(m_schoof+1)*C_schoof.*sqrt(u.^2+u_eps.^2).^(m_schoof+1);
  otherwise
    fout = u.^2/2;	%default linear law
end

end

function fout = TauB(u,parameters)
%basal shear stress in an element

switch parameters.frictionlaw
  case 'Coulomb'
    %unpack parameters
    taub = parameters.taub;
    u_eps = parameters.u_eps;   %regularizing velocity
    fout = taub.*u./sqrt(u.^2+u_eps.^2);
  case 'Weertman'
    %unpack
    C_schoof = parameters.C_schoof;
    m_schoof = parameters.m_schoof;
    u_eps = parameters.u_eps;
    fout = C_schoof.*sqrt(u.^2+u_eps.^2).^(m_schoof-1).*u;
  otherwise
    fout = u;	%default linear law
end

end

function fout = DTauB(u,parameters)
%derivative w.r.t u of basal shear stress in an element

switch parameters.frictionlaw
  case 'Coulomb'
    %unpack parameters
    taub = parameters.taub;
    u_eps = parameters.u_eps;   %regularizing velocity
    fout = taub./sqrt(u.^2+u_eps.^2) - (taub.*u.^2)./((u.^2+u_eps.^2).^(3/2));
  case 'Weertman'
    %unpack
    C_schoof = parameters.C_schoof;
    m_schoof = parameters.m_schoof;
    u_eps = parameters.u_eps;
    fout = C_schoof.*(sqrt(u.^2+u_eps.^2).^(m_schoof-1) + (m_schoof-1)*sqrt(u.^2+u_eps.^2).^(m_schoof-3).*u.^2);
  otherwise
    fout = 0;	%default linear law
end

end

function Gout = DupontG(x,parameters)
%Dupont lateral shear stress coefficient, may depend on local width of
%stream (and hence on x)
halfwidth = width_fun(x,parameters)./2;

Gout = parameters.B_shear./(halfwidth.*(parameters.width_shear^(1/parameters.n_Glen)));
end
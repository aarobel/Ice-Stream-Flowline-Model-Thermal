function parameters = setparams(parameters)
%list of parameters

%% Time step parameters
parameters.dtau = parameters.tfinal/parameters.nsteps; %length of time steps
parameters.dtau_max = parameters.dtau;

%% Newton Parameters
parameters.HS_sensitivity = pi*parameters.year/10;     %sensitivity of the HS function (as this gets larger, theta approaches the actual HS function)
parameters.uverbose = 1;
parameters.iteration_threshold = 1e-3;
parameters.hiter_max=1e3;
parameters.uiter_max=5e2;
parameters.titer_max=4e1;
parameters.CFL_max = 4;

%% Grid Parameters
%USE THIS FUNCTION ALONG WITH ANTICIPATED X_G TO FIND EXPONENT FOR GRID
%REFINEMENT:
%fzero(@(x) (1/parameters.grid.n_nodes).^(x) - (250/x_g),1) 

parameters.grid.n_elements = parameters.grid.n_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
% parameters.grid.sigma_node = linspace(0,1,parameters.grid.n_nodes)';  %node positions scaled to (0,1)
% parameters.grid.sigma_node = flipud(1-linspace(0,1,parameters.grid.n_nodes)'.^parameters.grid.n_exponent); %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_node = [linspace(0,0.97,parameters.grid.n_nodes-parameters.grid.gz_nodes),linspace(0.97+(.03/parameters.grid.gz_nodes),1,parameters.grid.gz_nodes)]'; %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_element =...
    (parameters.grid.sigma_node(1:parameters.grid.n_nodes-1)+...
    parameters.grid.sigma_node(2:parameters.grid.n_nodes))/2;     %element centres scaled to (0,1)

parameters.grid.n2_elements = parameters.grid.n2_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
parameters.grid.eta_node = linspace(0,1,parameters.grid.n2_nodes)';  %eta node positions scaled to (0,1)
parameters.grid.eta_element =...
    (parameters.grid.eta_node(1:parameters.grid.n2_nodes-1)+...
    parameters.grid.eta_node(2:parameters.grid.n2_nodes))/2;     %eta element centres scaled to (0,1)

parameters.regrid_on = 0;

%% Glen's Law parameters
% parameters.B_Glen = (1e-24)^(-1/3) .* ones(parameters.grid.n_elements,1);                     %B in Glen's law (vertically averaged if necessary)
parameters.n_Glen = 3;

parameters.R = 8.3145;
parameters.Qminus = 60e3;
parameters.Qplus =139e3;
parameters.Aminus = 4e-13;
parameters.Aplus = 1.734e3;
parameters.MP = 273.15;

%% Heat budget parameters

parameters.G = 0.07;
parameters.T_s = -15.*ones(parameters.grid.n_nodes-1,1);

%% Physical parameters
parameters.rho = 917;  %917                                 %ice density
parameters.rho_w = 1028;  %1028                               %water density
parameters.g = 9.81;                                    %acceleration due to gravity
parameters.D_eps = 1e-10;                               %strain rate regularizer
parameters.u_eps = 1e-9;                %velocity regularizer
parameters.u_in = 0./parameters.year; 

%% Sliding Law Parameters
parameters.frictionlaw = 'Coulomb';

parameters.C_schoof = 0; %7.624e6;      %See Schoof (2007)
parameters.m_schoof = 1/3;          %See Schoof (2007)

parameters.B_shear = (1e-24)^(-1/parameters.n_Glen);
parameters.width_shear = 1e3;

parameters.float = 1;
parameters.buttress = 0;

%% Till parameters
parameters.h_till_max=4;        %max till thickness
parameters.e_min=0.5;           %void ratio consolidation threshold

parameters.k_i = 2.1;           %thermal conductivity of ice

parameters.a_till = 9.44e8;     %empirical coefficient a from Tulaczyk
parameters.b_till = -21.7;      %empirical coefficient b from Tulaczyk
parameters.L_f = 3.35e5;        %specific latent heat of ice
parameters.C_i = 1.938e6;       %volumetric heat capacity of ice

parameters.HKappa = 0;
parameters.MWKappa = 0;     %MW diffusion turned off
parameters.Kappa = 1.41e-6; %See p.400 of Cuffey & Patterson
parameters.bedHflux = -parameters.G/parameters.k_i; %-3e-2 - 0.1.*exp(-40.*(parameters.grid.sigma_element-0.5).^2);
parameters.dsHflux = 0;
parameters.usHflux = 0;

parameters.e_max=1.0;

end
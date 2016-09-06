%% Ice stream flowline model
% This flowline model solves for ice thickness, grounding line
% position, velocity, void ratio, unfrozen till thickness and
% internal ice temperature using an operator splitting approach.
% Thickness and grounding line position are calculated simultaneously
% using a backward Euler method. Velocity and temperature are
% calculated separately also using implicit methods. Void ratio and
% till thickness are calculated explicitly using forward Euler
% to allow for ad-hoc corrections at the consolidation threshold.
%
% See Robel, Schoof, Tziperman (2014) in JGR-Earth Surface for details
%
% This code has been written by Alex Robel and Christian Schoof
% over the period 2011-2014.
% This version of the code has been prepared for public distribution

warning('off')
plotting_on=1;
%% Set Parameters
parameters.grid.n_nodes = 200;                      %Horizontal Resolution
% parameters.grid.n_exponent = 1.34;                %Horizontal Resolution refinement exponent
parameters.grid.n2_nodes = 20;                      %Vertical Resolution
parameters.grid.gz_nodes = 50;                      %Horizontal Resolution in grounding zone

parameters.year = 3600*24*365;                      %length of a year in seconds
parameters.tfinal = 1e3.*parameters.year;          %total time of integration
parameters.nsteps = 1e3;                           %number of time steps

parameters.icedivide = 100;                                 %bed elevation at upstream domain boundary
parameters.bedslope = -5e-4;                                %bed slope
parameters.width = 50e3;                                    %ice stream width (used for calc lateral shear stress)
parameters.accumrate = (0.3/parameters.year);               %accumulation rate
% parameters.sill_slope = sill_slope;
% parameters.sill_min = sill_min;
% parameters.sill_max = sill_min+sill_length;

parameters=setparams_init(parameters);              %set all other parameters
n_plots = 1e2;                                      %number of plots to make

%% Pre-allocate storage
e_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
htill_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
MWHF_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
h_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
xg_all = nan*ones(1,round(parameters.nsteps/5));
time_all = nan*ones(1,round(parameters.nsteps/5));
ub_all = nan*ones(parameters.grid.n_nodes,round(parameters.nsteps/5));

%% Initialize variables
disp('Initializing Variables...')
x_g = 800e3;                                                    %initial grounding line position  
h_g = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);      %initial grounding line thickness
h = 900 - (900-h_g).*parameters.grid.sigma_element.^2;               %initial thickness   

T_s = parameters.T_s;                                                %initial temperature  
T=zeros(parameters.grid.n_elements,parameters.grid.n2_elements);
for i=1:parameters.grid.n_elements
    T(i,:) = T_s(i)-(T_s(i)).*(1-parameters.grid.eta_element').^2;
end

e=0.6.*ones(size(h));                                                %initial void ratio
h_till = parameters.h_till_max.*ones(size(e));                       %initial unfrozen till thickness

parameters.taub = parameters.a_till.*exp(parameters.b_till.*e); 
parameters.B_Glen_full = set_B_Glen(T,parameters);
parameters.B_Glen = mean(parameters.B_Glen_full,2);

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_v5(linspace(0,10,parameters.grid.n_nodes)'./parameters.year,h,x_g,parameters); %initialize basal velocity
[u_old_full,u_old_mean] = add_def_vel_v2(u,h,x_g,parameters);        %initialize deformational velocity
v = mass_continuity(h,x_g,u_old_full,parameters);                    %initialize vertical velocities
parameters.hx_g_old = [h;x_g];

%% Time integration
parameters.uverbose = 0;
time = 0;
t=0;
check=1;
q=0;

while(parameters.year*time <= parameters.tfinal)
    
    t=t+1;

    %calculate velocity, thickness and GL position 
    %using a time step that is reduced if either of
    %the solvers reaches max iterations or CFL condition isn't met
    parameters.hx_g_old = [h;x_g];
    u_old = u;
    e_old = e;
    h_till_old = h_till;
    T_original = T;
    
    hiter=parameters.hiter_max;
    parameters.dtau = parameters.dtau_max; %start with max time step length
    courant = parameters.CFL_max;
    
    while(hiter==parameters.hiter_max || courant >= parameters.CFL_max)
          
        %calculate MW production and till properties (from previous time step)
        [e,h_till,MW_HF,parameters] = mw_solve_v4(T_original,parameters.hx_g_old(1:end-1),parameters.hx_g_old(end),u_old,e_old,h_till_old,parameters);
        parameters.taub = parameters.a_till.*exp(parameters.b_till.*e);

        %calculate temperature-dependent values of Glen's law parameters
        parameters.B_Glen_full = set_B_Glen(T_original,parameters);
        parameters.B_Glen = mean(parameters.B_Glen_full,2);
        
        %add in velocity due to deformation and calculate vertical velocity
        %with mass continuity
        [u_old_full,u_old_mean] = add_def_vel_v2(u_old,parameters.hx_g_old(1:end-1),parameters.hx_g_old(end),parameters);
        v = mass_continuity(parameters.hx_g_old(1:end-1),parameters.hx_g_old(end),u_old_full,parameters);
        
        %calculate thickness and GL position (with velocity from previous
        %time step)
        [h,x_g,hiter] = thickness_wGL_solve_v6(parameters.hx_g_old(1:end-1),parameters.hx_g_old(end),u_old_mean,parameters);
        
        %calculate temperature field (with current thickness and GL pos)
        T = temp_solve_v8(T_original,h,x_g,u_old_full,v,parameters);  
        
        %calculate velocity (with current thickness and GL pos)
        [u,error,uiter] = velocity_solve_v5(u_old,h,x_g,parameters);
      
        %Calculate CFL number
        u_eff = 0.5.*(u_old_full(1:end-1,:)+u_old_full(2:end,:)) - repmat(parameters.grid.sigma_element,[1 parameters.grid.n2_nodes]).*(x_g-parameters.hx_g_old(end))./parameters.dtau;
        v_eff = v - repmat(parameters.grid.eta_node',[parameters.grid.n_elements 1]).*repmat((h-parameters.hx_g_old(1:end-1))./parameters.dtau,[1 parameters.grid.n2_nodes]);

        courant1 = max(max(abs(u_eff./(repmat(x_g.*diff(parameters.grid.sigma_node)./parameters.dtau,[1 parameters.grid.n2_nodes])))));
        courant2 = max(max(abs(v_eff./(repmat(h.*diff(parameters.grid.eta_node(1:2))./parameters.dtau,[1 parameters.grid.n2_nodes])))));
        courant = max([courant1 courant2]);
        
        %If necessary, halve time step so that hx_g solver converges or CFL
        %is satisfied
        if(hiter==parameters.hiter_max || courant >= parameters.CFL_max)
            parameters.dtau = parameters.dtau/2; 
            disp(['Adapting time step...',int2str(uiter),...
                ' u iterations; ',int2str(hiter),' h iterations; Courant Number: ',num2str(courant)]);
        end
        
    end
    
    %plot solution on occasional time steps
    if((mod(t,parameters.nsteps/n_plots)==0 || t==1) && plotting_on==1)
        FlowlinePlot
        clc
    end
    
    %save some solutions
    if(mod(t,1)==0) %save solution on some steps
        
        q=q+1;
        
        e_all(:,q) = e;
        htill_all(:,q) = h_till;
        MWHF_all(:,q) = MW_HF;
        h_all(:,q) = h;
        xg_all(q) = x_g;
        time_all(q) = time;
    	ub_all(:,q) = u;      
        
    end
    
    time = time + (parameters.dtau./parameters.year);
    disp(['Time = ',num2str(time),' years; ',...
        int2str(uiter),' u iterations; ',int2str(hiter),' h iterations; Courant Number: ',num2str(courant)])
    
    if(sum(isnan(u))>0);disp('Integration FAIL');break;end
    if(parameters.dtau <= (2^(-20)).*parameters.dtau_max);disp('Integration FAIL');break;end
    
end

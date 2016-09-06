function fout = theta(x_g,x_g_old,u,sigma_node,parameters)
%Heavyside for ensuring good behavior of grounding line

HS_sensitivity = parameters.HS_sensitivity;

%unpack time step
dtau = parameters.dtau;
%grounding line migration rate
x_g_deriv = (x_g - x_g_old)./dtau;

heavy_var = u - sigma_node.*x_g_deriv;

fout = 0.5.*(1+tanh(HS_sensitivity.*heavy_var)); % here we use a hyperbolic tan function as a regularized hevyside function
% fout = 1.*ones(size(u)); %set theta to be a true upwind scheme
end
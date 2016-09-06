function fout = theta_z(heavy_var,parameters)
%Heavyside for ensuring good behavior of thickness stretch

HS_sensitivity = parameters.HS_sensitivity*1000;

fout = 0.5.*(1+tanh(HS_sensitivity.*heavy_var)); % here we use a hyperbolic tan function as a regularized hevyside function
% fout = 1.*ones(size(v)); %set theta to be a true upwind scheme
end
function fout = dwidthdx(x,parameters)
%derivative of ice stream width; needs to be consistent with that used
%in evolution problem for h


% width_int = width_fun(x,parameters);
% 
% fout = zeros(size(width_int));
% fout(2:end-1) = (width_int(3:end)-width_int(1:end-2))./(x(3:end)-x(1:end-2));
% fout(end) = (width_int(end)-width_int(end-1))./(x(end)-x(end-1));
% fout(1) = fout(2);

% bed_x = parameters.bed_x;
% width = parameters.width;
% if(x(end)>bed_x(end))
%     width = [width(1:end-1);width(end) + parameters.outer_widthslope.*linspace(0,x(end)-bed_x(end),100)'];
%     bed_x = [bed_x(1:end-1);linspace(bed_x(end),x(end),100)'];
%     
% end
% 
% width_interpolant = interp1(bed_x,width,'cubic','pp'); 
% dwdx_interpolant = fnder(width_interpolant,1);
% fout = ppval(dwdx_interpolant,x);
fout = zeros(size(x));
end
function fout = width_fun(x,parameters)
%width of ice stream; needs to be consistent with that used
%in evolution problem for h

% fout=nan.*ones(size(x));
% bed_x = parameters.bed_x;
% width = parameters.width;

% if(x(end)>bed_x(end))
%     width = [width(1:end-1);width(end) + parameters.outer_widthslope.*linspace(0,x(end)-bed_x(end),100)'];
%     bed_x = [bed_x(1:end-1);linspace(bed_x(end),x(end),100)'];
%     
% end
% 
% fout = interp1(bed_x,width,x);

%pw cubic interpolation of available width data
% width_interpolant = interp1(bed_x,width,'cubic','pp'); 
% fout = ppval(width_interpolant,x);

fout = parameters.width.*ones(size(x));

end
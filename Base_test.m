function fout = Base_test(x,parameters)
%elevation of base above sea level; needs to be consistent with that used
%in evolution problem for h

fout = parameters.icedivide + (parameters.bedslope*x);     %linear bed profile

% fout = (729 - (2184.8.*(x/750e3).^2) + (1031.72.*(x/750e3).^4) -(151.72.*(x/750e3).^6));
% fout = (parameters.icedivide - (parameters.bedparam1.*(x/750e3).^2) + (1031.72.*(x/750e3).^4) -(151.72.*(x/750e3).^6));
end
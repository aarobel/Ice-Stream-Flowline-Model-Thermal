function fout = Base(x,parameters)
%elevation of base w.r.t sea level; needs to be consistent with that used
%in evolution problem for h

fout = parameters.icedivide + (parameters.bedslope*x);

end
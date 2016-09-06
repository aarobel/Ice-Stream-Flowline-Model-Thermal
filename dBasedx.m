function fout = dBasedx(x,parameters)
%derivative of elevation of base w.r.t sea level; needs to be consistent with that used
%in evolution problem for h

fout=parameters.bedslope.*ones(length(x),1);

end
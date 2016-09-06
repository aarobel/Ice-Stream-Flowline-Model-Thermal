function B_Glen = set_B_Glen(T,parameters)
%Sets Glen's law parameter based on mean column ice temperature

B_Glen = zeros(size(T));

% B_Glen = B_Glen + (parameters.Aminus.*exp(-parameters.Qminus./(parameters.R.*(T+parameters.MP)))).^(-1/parameters.n_Glen);
i = find(T<-10);
i2 = find(T>=-10);
B_Glen(i) = (parameters.Aminus.*exp(-parameters.Qminus./(parameters.R.*(T(i)+parameters.MP)))).^(-1/parameters.n_Glen);
B_Glen(i2) = (parameters.Aplus.*exp(-parameters.Qplus./(parameters.R.*(T(i2)+parameters.MP)))).^(-1/parameters.n_Glen);

% B_Glen = 1e8.*ones(size(T));
B_Glen = wrap(min([B_Glen(:)';3e8.*ones(size(B_Glen(:)))']),parameters.grid.n_elements);
% B_Glen = wrap(max([B_Glen(:)';9e7.*ones(size(B_Glen(:)))']),parameters.grid.n_elements);
end
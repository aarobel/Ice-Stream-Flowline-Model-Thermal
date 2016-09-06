function z=wrap(x,m)

% Function z=wrap(x,m) ***
%   Where <x> is M rows by N columns, wraps each row of <x> into
%   m successive new rows within a new matrix <z>, such that <z>
%   will have dimensions M*m rows by N/m columns. 
%   E.g., for Z=WRAP(X,2):
%
%   X =  1  6  2  7  3  8  4  9  5 10  ==>  Z =  1  2  3  4  5
%       11 16 12 17 13 18 14 19 15 20            6  7  8  9 10
%                                               11 12 13 14 15
%                                               16 17 18 19 20
%
%   Uses RESHAPE. See also DEWRAP. 

[nr,nc]=size(x);
n=nc/m;z=zeros(nr*m,n);
for j=1:nr
  j1=(j-1)*m+1;j2=j*m;
  z(j1:j2,:)=reshape(x(j,:),m,n);
end

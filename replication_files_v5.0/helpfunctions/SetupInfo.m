function info = SetupInfo(nvar, npredetermined, Z, h)
%
%
%

% computes the dimension of space in which the spheres live that map into
% the orthogonal matrices satisfying the restriction.  
dim=0;
for i=1:nvar
    dim=dim+nvar-(i-1+size(Z{i},1));
end

% number zero restrictions
nzeros=0;
for i=1:nvar
    nzeros=nzeros+size(Z{i},1);
end

% gets random W
W=cell(nvar,1);
for j=1:nvar
    W{j}=randn(nvar-(j-1+size(Z{j},1)),nvar);
end

% info
info.nvar=nvar;
info.npredetermined=npredetermined;
info.Z=Z;
info.h=h;
info.dim=dim;
info.nzeros=nzeros;
info.W=W;

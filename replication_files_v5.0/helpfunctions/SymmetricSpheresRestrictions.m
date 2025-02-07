function r = SymmetricSpheresRestrictions(x, info)
%
%
%

nvar=info.nvar;
npredetermined=info.npredetermined;
W=info.W;

Sigma=reshape(x(nvar*npredetermined+1:nvar*(npredetermined+nvar)),nvar,nvar);
w=x(nvar*(npredetermined+nvar)+1:end);

r=zeros(nvar*(nvar-1)/2 + nvar,1);

m=1;

% Sigma must be symmetric
for i=1:nvar
    for j=i+1:nvar
        r(m)=Sigma(i,j) - Sigma(j,i);
        m=m+1;
    end
end

% w(i) must be of norm one
k=1;
for i=1:nvar
    s=size(W{i},1);
    r(m)=norm(w(k:k+s-1))^2 - 1.0;
    k=k+s;
    m=m+1;
end
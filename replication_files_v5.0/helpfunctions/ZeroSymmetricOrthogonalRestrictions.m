function y = ZeroSymmetricOrthogonalRestrictions(x, info)

nvar=info.nvar;
npredetermined=info.npredetermined;
nzeros=info.nzeros;
ZF=info.ZF(f_h_inv(x,info),info);

y=zeros(nvar*nvar+nzeros,1);

% Sigma must be symmetric
Sigma = reshape(x(nvar*npredetermined+1:nvar*npredetermined+nvar*nvar),nvar,nvar);
k=1;
for i=1:nvar
    for j=i+1:nvar
        y(k)=Sigma(i,j) - Sigma(j,i);
        k=k+1;
    end
end

% Q must be orthogonal
Q = reshape(x(nvar*(npredetermined+nvar)+1:end),nvar,nvar);
for i=1:nvar
    y(k)=Q(:,i)'*Q(:,i) - 1.0;
    k=k+1;
    for j=i+1:nvar
        y(k)=Q(:,i)'*Q(:,j);
        k=k+1;
    end
end

% zero restrictions
for i=1:nvar
    s=size(ZF{i},1);
    y(k:k+s-1)=ZF{i}(:,i);
    k=k+s;
end
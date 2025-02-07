function z = ZeroOrthogonalRestrictions(x, info, B, Sigma)
%Z_IRF, nvar)
%
% zero restrictions on Q=reshape(x,nvar,nar), which are of the form
%
%   Z_IRF{i}*Q(:,i) = 0
%
% and the restrictions that Q be orthogonal.
%

nvar=info.nvar;
total_zeros=nvar*(nvar+1)/2 + info.nzeros;
ZF=info.ZF(f_h_inv([vec(B); vec(Sigma); vec(eye(nvar))],info),info);

% restrictions
Q = reshape(x,nvar,nvar);
z=zeros(total_zeros,1);
k=1;
for i=1:nvar
    s=size(ZF{i},1);
    if s > 0
        z(k:k+s-1)=ZF{i}*Q(:,i);
        k=k+s;
    end

    for j=1:i-1
        z(k)=dot(Q(:,i),Q(:,j));
        k=k+1;
    end
    
    z(k)=dot(Q(:,i),Q(:,i)) - 1.0;
    k=k+1;
end
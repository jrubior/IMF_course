function q = SpheresToQ(w, info, B, Sigma)
%Z_IRF, nvar, W)
%
%  w=[w(1); ... ; w(n)] where the dimension of w(j) is n-(j-1+z(j)) > 0. 
%
%  Z_IRF{j} - z(j) x n matrix of full row rank. 
%      Z_IRF{j} = Z{j} * F(f_h^{-1}(B,Sigma,I_n)).
%
%  W{j} - (n-(j-1+z(j))) x n matrix.
%
%  q=[q(1); ... ; q(n)] where the dimension of q(j) is n.
%
%  norm(q(j)) == 1 if and only if norm(x(j)) == 1.
%

nvar=info.nvar;
W=info.W;
ZF=info.ZF(f_h_inv([vec(B); vec(Sigma); vec(eye(nvar))],info),info);

Q=zeros(nvar,nvar);
k=0;
for j=1:nvar
    s=size(W{j},1);
    wj=w(k+1:k+s);
    Mj_tilde=[Q(:,1:j-1)'; ZF{j}; W{j}];
    [K,R]=qr(Mj_tilde');    
    for i=nvar-s+1:nvar
        if (R(i,i) < 0) 
            K(:,i)=-K(:,i);
        end
    end 
    Kj=K(:,nvar-s+1:nvar);
    Q(:,j)=Kj*wj;
    k=k+s;
end
q=reshape(Q,nvar*nvar,1);

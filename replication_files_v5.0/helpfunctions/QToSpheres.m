function w = QToSpheres(q, info, B, Sigma)
%
%  q - n^2 dimensional vector
%
%  Z_IRF{j} - z(j) x n matrix of full row rank. 
%      Z_IRF{j} = Z{j} * F(f_h^{-1}(B,Sigma,I_n)).
%
%  W{j} - (n-(j-1+z(j))) x n matrix.
%
%  w=[w(1); ... w(n)] where the dimension of w(j) is n-(j-1+z(j)) > 0. 
%
%  norm(w(j)) == 1 if Q is orthogonal and satisfies the zero restrictions.
%

nvar=info.nvar;
dim=info.dim;
W=info.W;

ZF=info.ZF(f_h_inv([vec(B); vec(Sigma); vec(eye(nvar))],info),info);

Q=reshape(q,nvar,nvar);
w=zeros(dim,1);
k=0;
for j=1:nvar
    s=size(W{j},1);
    Mj_tilde=[Q(:,1:j-1)'; ZF{j}; W{j}];
    [K,R]=qr(Mj_tilde');    
    for i=nvar-s+1:nvar
        if (R(i,i) < 0) 
            K(:,i)=-K(:,i);
        end
    end 
    Kj=K(:,nvar-s+1:nvar);
    w(k+1:k+s)=Kj'*Q(:,j);
    k=k+s;
end
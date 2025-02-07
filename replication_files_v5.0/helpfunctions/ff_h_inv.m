function y = ff_h_inv(x,info)
%
% Mapping from (B, Sigma, W) to (A0,A+).  Is Sigma is symmetric and
% positive definite, then (A0,A+) will satisfy the zero restrictions and
% ff_h(ff_h_inv(B,Sigma,W)) = (B,Sigma,W)
%

nvar=info.nvar;
npredetermined=info.npredetermined;
z=f_h_inv([x(1:(npredetermined+nvar)*nvar); vec(eye(nvar))],info);
ZF=info.ZF(z,info);
W=info.W;

w=x((npredetermined+nvar)*nvar+1:end);
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

A0=reshape(z(1:nvar*nvar),nvar,nvar);
Aplus=reshape(z(nvar*nvar+1:end),npredetermined,nvar);

y=[vec(A0*Q); vec(Aplus*Q)];



function y = ff_h(x,info)
%
% Mapping from (A0,A+) to (B,Sigma,W).  If (A0,A+) satisfies the zero
% restrictions, then ff_h_inv(ff_h(A0,A+)) = (A0,A+)
%

nvar=info.nvar;
npredetermined=info.npredetermined;
W=info.W;
dim=info.dim;
z=f_h(x,info);
ZF=info.ZF(f_h_inv([z(1:nvar*(npredetermined+nvar)); vec(eye(nvar))],info),info);

Q=reshape(z(nvar*(npredetermined+nvar)+1:end),nvar,nvar);

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

 y=[z(1:nvar*(npredetermined+nvar)); w];


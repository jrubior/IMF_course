function w = DrawW(info)
%
%
%

nvar=info.nvar;
dim=info.dim;
W=info.W;

w=zeros(dim,1);
k=1;
for j=1:nvar
    s=size(W{j},1);
    wj=randn(s,1);
    w(k:k+s-1)=wj/norm(wj);
    k=k+s;
end
    
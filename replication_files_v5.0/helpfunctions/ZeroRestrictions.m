function z=ZeroRestrictions(a, info)
%
%
%

nvar=info.nvar;
nzeros=info.nzeros;
ZF=info.ZF(a,info);

z=zeros(nzeros,1);
k=1;
for i=1:nvar
    s=size(ZF{i},1);
    z(k:k+s-1)=ZF{i}(:,i);
    k=k+s;
end
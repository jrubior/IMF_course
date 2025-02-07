function y = StructuralRestrictions(x,Z)
%
%  Z - n x 1 cell of z_j x k matrices
%  f(x) - stacked A0, impulse responses, or both
%  
%  restrictions:
%
%     Z{j}*f(x)*e_j = 0
%

n=size(Z,1);
total_zeros=0;
for j=1:n
    total_zeros=total_zeros+size(Z{j},1);
end


A0           = reshape(x(1:n*n),n,n);
L0           = inv(A0)';
f            = L0;

y=zeros(total_zeros,1);
ib=1;
for j=1:n
    ie=ib+size(Z{j},1);  
    y(ib:ie-1)=Z{j}*f(:,j);
    ib=ie;
end
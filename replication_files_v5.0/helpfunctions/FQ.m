function y = FQ(x)

n=3;
A0=reshape(x(1:n*n),n,n);
Sigma=inv(A0*A0');
Q=chol(Sigma)*A0;
y=vec(Q);
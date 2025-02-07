function y = f_h(x,info)
%
% mapping from vectorized structural parameters to vectorized othogonal
% reduced-form parameters.
%

% h_tilde() also uses info.h
nvar=info.nvar;
npredetermined=info.npredetermined;

A0=reshape(x(1:nvar*nvar),nvar,nvar);
Aplus=reshape(x(nvar*nvar+1:end),npredetermined,nvar);

B=Aplus/A0;
Sigma=inv(A0*A0');
Q=h_tilde(Sigma,info)*A0;

y=[vec(B); vec(Sigma); vec(Q)];

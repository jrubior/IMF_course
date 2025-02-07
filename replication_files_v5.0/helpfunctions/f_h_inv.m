function y = f_h_inv(x,info)
%
% Mapping from an open set containing the orthogonal reduced-form
% parameters to the structural parameters.  Is the inverse of f_h over the
% orthogonal reduced-form parameters.
%

nvar=info.nvar;
npredetermined=info.npredetermined;

B=reshape(x(1:npredetermined*nvar),npredetermined,nvar);
Sigma=reshape(x(npredetermined*nvar+1:(npredetermined+nvar)*nvar),nvar,nvar);
Q=reshape(x((npredetermined+nvar)*nvar+1:end),nvar,nvar);

A0=h_tilde(Sigma,info)\Q;
Aplus=B*A0;

y=[vec(A0); vec(Aplus);];
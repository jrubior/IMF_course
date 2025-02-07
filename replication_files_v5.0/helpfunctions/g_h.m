function y = g_h(x, info)
%
% Mapping from (B, Sigma, Q) to (B, Sigma, W)
%

nvar=info.nvar;
npredetermined=info.npredetermined;

B=reshape(x(1:nvar*npredetermined),npredetermined,nvar);
Sigma=reshape(x(nvar*npredetermined+1:nvar*(npredetermined+nvar)),nvar,nvar);

w=QToSpheres(x(nvar*(npredetermined+nvar)+1:end),info,B,Sigma);

y=[x(1:nvar*(npredetermined+nvar)); w];
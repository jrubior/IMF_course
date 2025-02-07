function Y = h_tilde(X, info)
%
% X - n x n matrix
%
% Y - n x n matrix such that Y'*Y = 0.5*(X+X')
%

Y=info.h(0.5*(X+X'));

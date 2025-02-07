function [c,ceq] = mycon(x)
c = [];     % Compute nonlinear inequalities at x.
ceq = x'*x-1;   % Compute nonlinear equalities at x.

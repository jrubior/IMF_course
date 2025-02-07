function x = vec(X)
%
% stacks the columns of X into a column vector
%

x=reshape(X,numel(X),1);
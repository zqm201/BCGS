function [V, D] = cs(X)
% Column scaling
% input:
%   X : matrix of shape (m, n)
% output:
%   V : matrix with normalized columns
%   D : diagonal matrix

%%
[m, s] = size(X);
D = zeros(s,s);
V = zeros(m,s);

D(1,1) = norm(X(:,1));
V(:,1) = X(:,1) / D(1,1);

for k = 1:s-1
    D(k+1,k+1) = norm(X(:,k+1));
    V(:,k+1) = X(:,k+1) / D(k+1,k+1);
end
end
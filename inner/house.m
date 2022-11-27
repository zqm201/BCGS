function [Q, R] = house(X)
% Householder orthogonalization
% references:
%   Golub, Van Loan. 2013.
% input:
%   X : matrix of shape (m, n)
% output:
%   Q : orthogonal matrix
%   R : upper triangular matrix

%%
[m, n] = size(X);
Q = eye(m, n);
beta = zeros(n, 1);

for k = 1:n
    v = X(k:end,k);
    sigma = v(2:end)' * v(2:end);
    mu = -sign(v(1)) * sqrt(v(1)^2 + sigma);
    v(1) = v(1) - mu;
    beta(k) = 2 * v(1)^2 / (sigma + v(1)^2);
    v = v / v(1);  % P = I - beta vv', v(1) = 1, beta = 2/(vv')
    X(k:m,k:n) = X(k:m,k:n) - (beta(k) * v) * (v' * X(k:m,k:n));
    X(k+1:m,k) = v(2:end);

    Q(k:m,k) = Q(k:m,k) - (beta(k) * v) * (v' * Q(k:m,k));
    for j = k-1:-1:1
        v = [1; X(j+1:m,j)];
        Q(j:m,k) = Q(j:m,k) - (beta(j) * v) * (v' * Q(j:m,k));
    end
end
R = triu(X(1:n,:));
end
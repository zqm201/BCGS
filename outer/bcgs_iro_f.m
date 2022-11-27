function [QQ, RR] = bcgs_iro_f(XX, s, fun1, fun2)
% Reorthogonalized block classical Gram-Schmidt
% references:
%   Barlow, Smoktunowicz. 2013.
% input:
%   XX   : block matrix of shape (m, n)
%   s    : block size
%   fun1 : first factorization function
%   fun2 : second factorization function
% output:
%   QQ   : block orthogonal matrix
%   RR   : block upper triangular matrix

%%
addpath(genpath('../inner/'))

[m, n] = size(XX);
QQ = zeros(m, n);
RR = zeros(n, n);
p = n / s;

kk = 1:s;
sk = 0;
[QQ(:,kk), RR(kk,kk)] = house(XX(:,kk));

for k = 2:p
    kk = kk + s;
    sk = sk + s;

    RR(1:sk,kk) = QQ(:,1:sk)' * XX(:,kk);
    W = XX(:,kk) - QQ(:,1:sk) * RR(1:sk,kk);
    [QQ(:,kk), RR(kk,kk)] = fun1(W);

    RR2 = QQ(:,1:sk)' * QQ(:,kk);
    W = QQ(:,kk) - QQ(:,1:sk) * RR2;
    [QQ(:,kk), R2] = fun2(W);

    RR(1:sk,kk) = RR(1:sk,kk) + RR2 * RR(kk,kk);
    RR(kk,kk) = R2 * RR(kk,kk);
end

end
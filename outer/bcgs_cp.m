function [QQ, RR] = bcgs_cp(XX, s, fun)
% Block classical Gram-Schmidt with continuous projection
% input:
%   XX  : block matrix of shape (m, n)
%   s   : block size
%   fun : local factorization function
% output:
%   QQ  : block orthogonal matrix
%   RR  : block upper triangular matrix

%%
addpath(genpath('../inner/'))

[m, n] = size(XX);
QQ = zeros(m, n);
RR = zeros(n, n);
p = n / s;

kk = 1:s;
sk = 0;
[QQ(:,kk), RR(kk,kk)] = fun(XX(:,kk));

for k = 2:p
    kk = kk + s;
    sk = sk + s;

    RR(1:sk,kk) = QQ(:,1:sk)' * XX(:,kk);
    W = XX(:,kk) - QQ(:,1:sk) * RR(1:sk,kk);

    RR2 = QQ(:,1:sk)' * W;
    W = W - QQ(:,1:sk) * RR2;
    [QQ(:,kk), RR(kk,kk)] = fun(W);

    RR(1:sk,kk) = RR(1:sk,kk) + RR2;
end

end
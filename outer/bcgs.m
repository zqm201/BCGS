function [QQ, RR] = bcgs(XX, s)
% Block classical Gram-Schmidt
% input:
%   XX : block matrix of shape (m, n)
%   s  : block size
% output:
%   QQ : block orthogonal matrix
%   RR : block upper triangular matrix

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
    [QQ(:,kk), RR(kk,kk)] = house(W);
end

end
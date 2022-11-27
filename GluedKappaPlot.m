clear
close all
format long

addpath(genpath('inner/'))
addpath(genpath('outer/'))
addpath(genpath('matrices/'))

%%
nalg = 6;  % BCGS, BCGSI+, BCGSI+F(MGS), BCGSI+F(CholQR), BCGS-CS, BCGS-CP
logcondXX = 1:8;
nmat = length(logcondXX);
loss_ortho = zeros(nmat, nalg);
XXcond = zeros(1, nmat);

rng(0)
XXdim = [10000, 50, 10];  % [m, p, s]
m = XXdim(1); p = XXdim(2); s = XXdim(3);
n = p * s;
I = eye(n);

% Fix glued dimensions
factors = factor(n);
mid_ind = round(length(n)/2);
r = prod(factors(1:mid_ind));
t = prod(factors(mid_ind+1:end));

for i = 1:nmat
    XX = CreateGluedMatrix(m, r, t, .5 * logcondXX(i), logcondXX(i));
    XXcond(i) = cond(XX);

    [QQ, ~] = bcgs(XX, XXdim(3));  % BCGSI+
    loss_ortho(i, 1) = norm(I - QQ' * QQ, 2);
    [QQ, ~] = bcgs_iro_f(XX, XXdim(3), @house, @house);  % BCGSI+
    loss_ortho(i, 2) = norm(I - QQ' * QQ, 2);
    [QQ, ~] = bcgs_iro_f(XX, XXdim(3), @mgs, @house);  % BCGSI+F(MGS)
    loss_ortho(i, 3) = norm(I - QQ' * QQ, 2);
    [QQ, ~] = bcgs_iro_f(XX, XXdim(3), @cholqr, @house);  % BCGSI+F(CholQR)
    loss_ortho(i, 4) = norm(I - QQ' * QQ, 2);
    [QQ, ~] = bcgs_iro_f(XX, XXdim(3), @cs, @house);  % BCGS-CS
    loss_ortho(i, 5) = norm(I - QQ' * QQ, 2);
    [QQ, ~] = bcgs_cp(XX, XXdim(3), @house);  % BCGS-CP
    loss_ortho(i, 6) = norm(I - QQ' * QQ, 2);
end

alg = {'BCGS(HouseQR)', 'BCGSI+(HouseQR)', 'BCGSI+F(MGS, HouseQR)', ...
       'BCGSI+F(CholQR, HouseQR)', 'BCGS-CS(HouseQR)', 'BCGS-CP(HouseQR)'};
color = {[228,26,28]/255, [55,126,184]/255, [77,175,74]/255, ...
         [255,127,0]/255, [152,78,163]/255, [166,86,40]/255};
lbl = {'o-', '*--', 's:', 'd-.', '+-', 'x--'};
hold on
for j = 1:nalg
    plot(XXcond, loss_ortho(:,j), lbl{j}, 'Color', color{j}, ...
         'LineWidth', 1.5, 'MarkerSize', 9);
end
grid on
set(gca, 'Yscale', 'log', 'Xscale', 'log', 'XMinorGrid', 'off', ...
         'YMinorGrid', 'off', 'FontSize', 13);
xlabel('\(\kappa(\mathcal{X})\)', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('\(||I-\widehat{\mathcal{Q}}^T\widehat{\mathcal{Q}}||\)', ...
       'Interpreter', 'latex', 'FontSize', 16)
lgd = legend(alg, 'FontSize', 14, 'Location', 'east');
ylim([1e-16,1e1])
saveas(gca, 'fig1b.eps', 'epsc')

clear
close all
format long

addpath(genpath('inner/'))
addpath(genpath('outer/'))
addpath(genpath('matrices/'))

%%
%  rand_uniform, rand_normal, rank_def, laeuchli, monomial,
%  stewart, stewart_extreme, hilbert, s-step, newton
XXdim = [10000, 50, 10];  % [m, p, s]
n = XXdim(2) * XXdim(3);
mat = {'rand_uniform', 'rand_normal', 'monomial', 's-step', 'newton'};
nummat = length(mat);

rng(0)
XX = cell(nummat, 1);

fID = fopen('tab.txt', 'w');

fprintf(fID, '------------------------------------------\n');
fprintf(fID, 'matrix ID | sigma_min | sigma_max | kappa \n');
fprintf(fID, '------------------------------------------\n');
for i = 1:nummat
    [XX{i}, ~, XXprops] = MatGen(mat{i}, XXdim);
    fprintf(fID, '%s\t (%.2e %.2e %.2e)\n', mat{i}, XXprops.sv(1), ...
            XXprops.sv(end), XXprops.cond);
end

fprintf(fID, '-----------------------\n');
fprintf(fID, 'BCGS:     LOO  |  Res  \n');
fprintf(fID, '-----------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs(XX{i}, XXdim(3));  % BCGS
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end
fprintf(fID, '-----------------------\n');
fprintf(fID, 'BCGSI+:   LOO  |  Res  \n');
fprintf(fID, '-----------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs_iro_f(XX{i}, XXdim(3), @house, @house);  % BCGSI+
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end
fprintf(fID, '----------------------------\n');
fprintf(fID, 'BCGSI+F(MGS):  LOO  |  Res  \n');
fprintf(fID, '----------------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs_iro_f(XX{i}, XXdim(3), @mgs, @house);  % BCGSI+F(MGS)
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end
fprintf(fID, '-------------------------------\n');
fprintf(fID, 'BCGSI+F(CholQR):  LOO  |  Res  \n');
fprintf(fID, '-------------------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs_iro_f(XX{i}, XXdim(3), @cholqr, @house);  % BCGSI+F(CholQR)
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end
fprintf(fID, '-----------------------\n');
fprintf(fID, 'BCGS-CS:  LOO  |  Res  \n');
fprintf(fID, '-----------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs_iro_f(XX{i}, XXdim(3), @cs, @house);  % BCGS-CS
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end
fprintf(fID, '-----------------------\n');
fprintf(fID, 'BCGS-CP:  LOO  |  Res  \n');
fprintf(fID, '-----------------------\n');
for i = 1:nummat
    [QQ, RR] = bcgs_cp(XX{i}, XXdim(3), @house);  % BCGS-CP
    norm_XX = norm(XX{i}, 2);

    loo = norm(eye(n) - QQ'*QQ, 2);
    res = norm(XX{i} - QQ*RR, 2) / norm_XX;
    fprintf(fID, '%.2e\t %.2e\t %s\n', loo, res, mat{i});
end

fclose(fID);

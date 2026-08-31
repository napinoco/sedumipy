% generate_symbcholden_oracle.m
%
% Oracle for this port's symbcholden.py (dense-columns optimization,
% Stage 3 of the plan). Builds a real L = symbchol(X) struct via the
% vendored symbchol.mex machinery (global ADA_sedumi_), then hand-builds
% small dense/DAt structs and calls the real vendored symbcholden.m
% (itself calling the real symbfwblk/incorder/finsymbden mex binaries)
% directly -- symbcholden.m only uses dense.l/length(dense.q)/
% length(dense.cols)/dense.A/DAt.denq, so no full getdense/pretransfo
% pipeline is needed to exercise it faithfully.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_symbcholden_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'symbcholden');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

global ADA_sedumi_

function S = run_case(X, dense, DAt)
    global ADA_sedumi_
    ADA_sedumi_ = X;
    L = symbchol();
    Lden = symbcholden(L, dense, DAt);
    S.X = X;
    S.L_perm = L.perm;
    S.L_L = L.L;
    S.L_xsuper = L.xsuper;
    S.dense_l = dense.l;
    S.dense_q = dense.q;
    S.dense_cols = dense.cols;
    S.dense_A = dense.A;
    S.DAt_denq = DAt.denq;
    S.Lden_LAD = Lden.LAD;
    S.Lden_perm = Lden.perm;
    S.Lden_dz = Lden.dz;
    S.Lden_first = Lden.first;
end

rand('seed', 19);

%% Case 1: m=8, nl=2 LP-dense cols, nq=1 Lorentz block, nden=2 norm-bound cols.
m1 = 8;
T1 = tril(sprandn(m1, m1, 0.3), -1);
X1 = T1 + T1' + m1 * speye(m1);
dense1.l = 2;
dense1.q = 1;
dense1.cols = zeros(5,1);  % only length used by symbcholden.m
dense1.A = sprand(m1, 5, 0.5);   % [2 LP][1 trace][2 normbound]
DAt1.denq = sprand(m1, 1, 0.6);
S1 = run_case(X1, dense1, DAt1);
save('-v7', fullfile(out_dir, 'case1.mat'), '-struct', 'S1');
fprintf('case1 Lden.perm: %s\n', mat2str(S1.Lden_perm'));
fprintf('case1 Lden.first: %s\n', mat2str(S1.Lden_first'));

%% Case 2: m=10, nl=0, nq=2 Lorentz blocks, nden=3 norm-bound cols.
m2 = 10;
T2 = tril(sprandn(m2, m2, 0.25), -1);
X2 = T2 + T2' + (m2+1) * speye(m2);
dense2.l = 0;
dense2.q = [1;2];
dense2.cols = zeros(5,1);  % nq(2) + nden(3)
dense2.A = sprand(m2, 5, 0.4);   % [0 LP][2 trace][3 normbound]
DAt2.denq = sprand(m2, 2, 0.5);
S2 = run_case(X2, dense2, DAt2);
save('-v7', fullfile(out_dir, 'case2.mat'), '-struct', 'S2');
fprintf('case2 Lden.perm: %s\n', mat2str(S2.Lden_perm'));
fprintf('case2 Lden.first: %s\n', mat2str(S2.Lden_first'));

%% Case 3: m=12, nl=3 LP-dense cols, nq=0 (no dense Lorentz blocks at all).
m3 = 12;
T3 = tril(sprandn(m3, m3, 0.2), -1);
X3 = T3 + T3' + (m3+1) * speye(m3);
dense3.l = 3;
dense3.q = zeros(0,1);
dense3.cols = zeros(3,1);  % nl(3) only
dense3.A = sprand(m3, 3, 0.5);
DAt3.denq = sparse(m3, 0);
S3 = run_case(X3, dense3, DAt3);
save('-v7', fullfile(out_dir, 'case3.mat'), '-struct', 'S3');
fprintf('case3 Lden.perm: %s\n', mat2str(S3.Lden_perm'));
fprintf('case3 Lden.first: %s\n', mat2str(S3.Lden_first'));

fprintf('symbcholden oracle written to %s\n', out_dir);

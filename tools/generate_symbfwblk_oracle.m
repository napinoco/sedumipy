% generate_symbfwblk_oracle.m
%
% Oracle for this port's _native.symbfwblk (dense-columns optimization,
% Stage 2 of the plan). Builds a real L = symbchol(X) struct (perm, L,
% xsuper) via the vendored symbchol.mex machinery (global ADA_sedumi_,
% exactly as symbchol.m expects), then calls the real vendored
% symbfwblk.mex directly on a random sparse B, matching how
% symbcholden.m uses it.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_symbfwblk_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'symbfwblk');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

global ADA_sedumi_

function S = run_case(X, B)
    global ADA_sedumi_
    ADA_sedumi_ = X;
    L = symbchol();
    Xpat = symbfwblk(L, B);
    S.X = X;
    S.B = B;
    S.L_perm = L.perm;
    S.L_L = L.L;
    S.L_xsuper = L.xsuper;
    S.Xpat = Xpat;
end

rand('seed', 17);

%% Case 1: chain-like sparse SPD pattern (m=8), multi-column sparse B.
m1 = 8;
T1 = tril(sprandn(m1, m1, 0.3), -1);
X1 = speye(m1) + T1 + T1' ;
X1 = X1 + m1 * speye(m1);  % ensure diag dominance / real symbchol path
B1 = sprand(m1, 3, 0.3);
S1 = run_case(X1, B1);
save('-v7', fullfile(out_dir, 'case1.mat'), '-struct', 'S1');
fprintf('case1 Xpat nnz=%d shape=%dx%d\n', nnz(S1.Xpat), size(S1.Xpat,1), size(S1.Xpat,2));

%% Case 2: block-diagonal-ish disconnected pattern (m=10), single-column B's.
m2 = 10;
X2 = speye(m2);
X2(1:4,1:4) = X2(1:4,1:4) + tril(ones(4),-1) + triu(ones(4),1);
X2(5:10,5:10) = X2(5:10,5:10) + tril(ones(6),-1) + triu(ones(6),1);
X2 = X2 + m2*speye(m2);
B2 = sparse(m2, 2);
B2(2,1) = 1;
B2(9,2) = 1;
S2 = run_case(X2, B2);
save('-v7', fullfile(out_dir, 'case2.mat'), '-struct', 'S2');
fprintf('case2 Xpat nnz=%d shape=%dx%d\n', nnz(S2.Xpat), size(S2.Xpat,1), size(S2.Xpat,2));

%% Case 3: larger random sparse SPD pattern (m=15), denser B.
m3 = 15;
T3 = tril(sprandn(m3, m3, 0.15), -1);
X3 = T3 + T3' + (m3+2) * speye(m3);
B3 = sprand(m3, 5, 0.4);
S3 = run_case(X3, B3);
save('-v7', fullfile(out_dir, 'case3.mat'), '-struct', 'S3');
fprintf('case3 Xpat nnz=%d shape=%dx%d\n', nnz(S3.Xpat), size(S3.Xpat,1), size(S3.Xpat,2));

fprintf('symbfwblk oracle written to %s\n', out_dir);

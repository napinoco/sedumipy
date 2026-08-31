% generate_deninfac_oracle.m
%
% Oracle for this port's deninfac.py dense-column branch (Stage 5 of the
% dense-columns optimization plan). Builds a real mixed L+Lorentz problem
% (same K/At/b/c shape as generate_getdatm_oracle.m), flags Lorentz block 1
% dense (same dense.cols=[3;5;6] layout generate_getdatm_oracle.m's own
% "dense2" case already established and verified: K2.mainblks=[3 5 10] =>
% trace row 3, K2.qblkstart=[5 7 10] => block 1's norm-bound rows 5:6),
% zeroes those rows out of At2 (matching sedumi.m's own
% `A(dense.cols,:)=0.0` pre-loop step) before building a real symbchol/
% blkchol L and symbcholden symLden, then calls the real vendored
% deninfac.m directly.
%
% Case 1 uses the default (large) pars.chol.maxuden -- dpr1fact.c
% shouldn't need to reorder for numerical stability here. Case 2 uses a
% deliberately tiny maxuden (1.01) to force dpr1fact's pivoting branch
% (dopiv=1), which test_cluster2.py's own docstring flags as otherwise
% untested anywhere in this port, plus a hand-injected L.skip entry with
% a huge pars.chol.canceltol/abstol to force deninfac.m's tail
% `Ld(skip(...))=1` assignment branch deterministically (that logic is
% pure elementwise post-processing of L.skip/absd/L.perm -- it doesn't
% care whether L.skip came from a genuine blkchol pivot decision).
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_deninfac_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'sdinit');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

global ADA_sedumi_

rand('seed', 202);

K.f = 0; K.l = 2; K.q = [3;4]; K.r = zeros(1,0); K.s = [2;3];
n = K.l + sum(K.q) + K.s(1)^2 + K.s(2)^2;
m = 4;
At = rand(n,m) - 0.5;
b = rand(m,1) - 0.5;
c = rand(n,1) - 0.5;
pars = checkpars(struct());

[At2, b2, c2, K2, prep] = pretransfo(At, b, c, K, pars);

dense0.cols = zeros(0,1);
dense0.A = sparse(length(b2), 0);
dense0.q = zeros(0,1);
dense0.l = 0;

[d, v, vfrm, y, y0, R] = sdinit(At2, b2, c2, dense0, K2, pars);
d.q1 = d.q1 .* (1 + 0.1*rand(size(d.q1)));
d.q2 = d.q2 + 0.05*(rand(size(d.q2))-0.5);

dense.l = 0;
dense.q = 1;
dense.cols = [3; 5; 6];
dense.A = sparse(full(At2(dense.cols, :))');

At2z = At2;
At2z(dense.cols, :) = 0.0;
Ablkjc = partitA(At2z, K2.mainblks);

DAtdenq_placeholder = sparse(ones(length(b2), 1));
DAt = getDAtm(At2z, Ablkjc, dense, DAtdenq_placeholder, d, K2);

ADA = zeros(length(b2), length(b2));
for j = 1:length(b2)
    ej = zeros(length(b2),1); ej(j) = 1;
    Aej = Amul(At2z, dense, ej, 1);
    PAej = PopK(d, Aej, K2);
    ADA(:,j) = Amul(At2z, dense, PAej, 0);
end
ADA = sparse((ADA+ADA')/2);
ADA_sedumi_ = ADA;

L = symbchol();
Ltmpsiz = L.tmpsiz;
[L.L, L.d, L.skip, L.add] = blkchol(L, ADA);
L.tmpsiz = Ltmpsiz;

symLden = symbcholden(L, dense, DAt);

absd = full(diag(ADA));

%% Case 1: default (large) maxuden -- no forced pivoting.
pars_chol1 = pars.chol;
pars_chol1.maxuden = 500;

[Lden1, Ld1] = deninfac(symLden, L, dense, DAt, d, absd, K2.qblkstart, pars_chol1);

S1.At2z = At2z; S1.K2 = K2; S1.d = d; S1.dense = dense; S1.DAt = DAt;
S1.L_L = L.L; S1.L_d = L.d; S1.L_skip = L.skip; S1.L_perm = L.perm; S1.L_xsuper = L.xsuper;
S1.symLden_LAD = symLden.LAD; S1.symLden_perm = symLden.perm;
S1.symLden_dz = symLden.dz; S1.symLden_first = symLden.first;
S1.absd = absd; S1.qblkstart = K2.qblkstart;
S1.pars_chol = pars_chol1;
S1.Lden_betajc = Lden1.betajc;
S1.Lden_beta = Lden1.beta;
S1.Lden_p = Lden1.p;
S1.Lden_dopiv = Lden1.dopiv;
S1.Lden_pivperm = Lden1.pivperm;
S1.Ld = Ld1;
save('-v7', fullfile(out_dir, 'deninfac_case1.mat'), '-struct', 'S1');
fprintf('case1: dopiv nnz = %d, betajc(end) = %d\n', nnz(Lden1.dopiv), Lden1.betajc(end));

%% Case 2: tiny maxuden to force dpr1fact pivoting, plus a hand-injected
%% L.skip entry (with huge canceltol/abstol) to force deninfac.m's tail
%% Ld(skip(...))=1 assignment branch.
L2 = L;
L2.skip = zeros(size(L2.skip));
L2.skip(1) = 1;

% Scale the dense-column DATA (not its sparsity pattern, which symLden
% was already built from and stays valid) way up, so dpr1fact.c's
% |p(i,k)*beta(j,k)| > maxu stability check actually fires.
dense2 = dense;
dense2.A = dense.A * 1e6;
DAt2 = DAt;
DAt2.denq = DAt.denq * 1e6;

pars_chol2 = pars.chol;
pars_chol2.maxuden = 1.01;
pars_chol2.canceltol = 1e10;
pars_chol2.abstol = 1e10;

[Lden2, Ld2] = deninfac(symLden, L2, dense2, DAt2, d, absd, K2.qblkstart, pars_chol2);

S2.At2z = At2z; S2.K2 = K2; S2.d = d; S2.dense = dense2; S2.DAt = DAt2;
S2.L_L = L2.L; S2.L_d = L2.d; S2.L_skip = L2.skip; S2.L_perm = L2.perm; S2.L_xsuper = L2.xsuper;
S2.symLden_LAD = symLden.LAD; S2.symLden_perm = symLden.perm;
S2.symLden_dz = symLden.dz; S2.symLden_first = symLden.first;
S2.absd = absd; S2.qblkstart = K2.qblkstart;
S2.pars_chol = pars_chol2;
S2.Lden_betajc = Lden2.betajc;
S2.Lden_beta = Lden2.beta;
S2.Lden_p = Lden2.p;
S2.Lden_dopiv = Lden2.dopiv;
S2.Lden_pivperm = Lden2.pivperm;
S2.Ld = Ld2;
save('-v7', fullfile(out_dir, 'deninfac_case2.mat'), '-struct', 'S2');
fprintf('case2: dopiv nnz = %d, betajc(end) = %d, Ld(skip)=%g\n', ...
    nnz(Lden2.dopiv), Lden2.betajc(end), Ld2(1));

fprintf('deninfac oracle written to %s\n', out_dir);

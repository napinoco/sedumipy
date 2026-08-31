% generate_pcg_dense_oracle.m
%
% Oracle for pcg.py's dense-column path (Stage 6 of the dense-columns
% optimization plan): wrapPcg.m/loopPcg.m driven with a genuine nonempty
% dense.cols/dense.q and a real (non-trivial) Lden from deninfac.m,
% rather than the trivial no-dense-columns Lden.betajc=1 that
% generate_pcg_oracle.m uses. Reuses the same K/At/dense-column setup as
% generate_deninfac_oracle.m (K2.mainblks=[3 5 10], K2.qblkstart=[5 7 10],
% dense.cols=[3;5;6], dense.q=1), then calls deninfac.m for a real Lden
% (default/large maxuden -- no forced pivoting, since the point here is
% to exercise loopPcg/wrapPcg's OWN new dense-column term, already
% covered by test_deninfac.py's own pivoting case), then drives
% wrapPcg.m/loopPcg.m exactly as generate_pcg_oracle.m already does.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_pcg_dense_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'sdinit');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

global ADA_sedumi_

rand('seed', 303);

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

ADA = zeros(m,m);
for j = 1:m
    ej = zeros(m,1); ej(j) = 1;
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

pars_chol = pars.chol;
pars_chol.maxuden = 500;

[Lden, Ld] = deninfac(symLden, L, dense, DAt, d, absd, K2.qblkstart, pars_chol);
L.d = Ld;

pars.cg = checkpars(struct()).cg;

% ---- wrapPcg ----
rb = rand(m,1) - 0.5;
rv = rand(K2.N,1) - 0.5;
[wy, wdx, wk, wr] = wrapPcg(L, Lden, At2z, dense, d, DAt, K2, rb, rv, pars.cg, min(1,y0)*R.maxRb);

% ---- loopPcg (fresh start, p=[]) ----
bvec = rand(m,1) - 0.5;
restol = 1e-10;
[ly, lk, lDAy] = loopPcg(L, Lden, At2z, dense, d, DAt, K2, bvec, [], 0, pars.cg, restol);

save('-v7', fullfile(out_dir, 'pcg_dense.mat'), ...
     'At2z', 'K2', 'd', 'dense', 'DAt', 'ADA', ...
     'rb', 'rv', 'y0', 'R', 'wy', 'wdx', 'wk', 'wr', ...
     'bvec', 'restol', 'ly', 'lk', 'lDAy');

Lstruct.L = L.L; Lstruct.d = L.d; Lstruct.xsuper = L.xsuper; Lstruct.perm = L.perm;
save('-v7', '-append', fullfile(out_dir, 'pcg_dense.mat'), 'Lstruct');

Ldenstruct.betajc = Lden.betajc; Ldenstruct.beta = Lden.beta; Ldenstruct.p = Lden.p;
Ldenstruct.dopiv = Lden.dopiv; Ldenstruct.pivperm = Lden.pivperm;
Ldenstruct.dz = Lden.dz; Ldenstruct.first = Lden.first; Ldenstruct.perm = Lden.perm;
save('-v7', '-append', fullfile(out_dir, 'pcg_dense.mat'), 'Ldenstruct');

fprintf('dense-column PCG oracle written to %s (wk=%d, lk=%d)\n', out_dir, wk, lk);

% generate_sedumi_dense_oracle.m
%
% Stage 8 of the dense-columns optimization plan: end-to-end verification
% that sedumi.py's full pipeline -- getdense() detection, the one-time
% A(dense.cols,:)=0.0 zeroing, symbcholden(), and deninfac.py/pcg.py's
% product-form dense-column correction -- produces the exact same (x,y,
% info) as the real Octave sedumi.m, on a problem where getdense.m's
% heuristic genuinely fires (unlike every existing tests/fixtures/sedumi/
% *.mat case, which are either small-and-fully-dense or -- case 7's
% sprand LP -- too uniformly sparse to trigger it).
%
% Background A is sparse (sprand); a single LP row and one whole Lorentz
% block (trace + norm-bound rows) are then overwritten to be fully
% dense, mirroring generate_getdense_oracle.m's case1/case2 construction.
% pars.denf is lowered from checkpars.m's default (10) to 3 (same value
% generate_getdense_oracle.m uses) so detection doesn't depend on a
% precariously-tuned random seed; pars.denf is a normal user-settable
% pars field so passing it explicitly to both the Octave and Python
% sedumi() calls is a faithful, symmetric test setup, not a special case.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_sedumi_dense_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'sedumi');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

pars.fid = 0;   % quiet
pars.denf = 3;  % aggressive dense-column threshold (checkpars.m default is 10)

function save_case(out_dir, name, At, b, c, K, pars)
    [x, y, info] = sedumi(At, b, c, K, pars);
    save('-v7', fullfile(out_dir, [name '.mat']), 'At', 'b', 'c', 'K', 'pars', 'x', 'y', 'info');
    fprintf('  %s: iter=%d numerr=%d pinf=%d dinf=%d cx=%g by=%g\n', ...
        name, info.iter, info.numerr, info.pinf, info.dinf, full(sum(c.*x)), full(sum(b.*y)));
end

% LP (15) + one Lorentz block (4: trace + 3 norm-bound) + one 2x2 PSD
% block (4), m=20 < n=23 (matching the m<n ratio the other _feasible
% cases in generate_sedumi_oracle.m use, which leaves the dual enough
% room to be strictly feasible instead of just barely closed). Row 2
% (an LP variable) and rows 16-19 (the whole Lorentz block) are
% overwritten to be fully dense; the LP/Lorentz background is
% sprand(0.06) and the Lorentz+PSD background is an extra-sparse
% sprand(0.02) so that getdense.m's `h` floor (max(NORMDEN=5, total
% PSD-block-row nnz)) stays pinned at NORMDEN=5 instead of being pushed
% up by the PSD block's own row density -- with h=5 fixed, the dense
% rows' full m-column density (colnz=20) clears denf*spquant=15
% comfortably. (A first attempt with a uniform 0.12 background, m=25,
% nl=10 did NOT trigger detection at all: the PSD block's row-nnz sum
% alone pushed h to 10, so denf*h=30 exceeded even the fully-dense
% rows' colnz=25. A second attempt, nl=10/m=30/seed=701, did trigger
% detection but landed on a dual-infeasible outcome (dinf=1) for every
% seed tried except a handful; seed=701 with nl=15/m=20 below was
% swept over seeds 1-30 and seed=21 is the one that lands on a clean
% pinf=0/dinf=0/numerr=0 solve with cx==by, confirmed by manually
% re-running getdense()+sedumi() across that seed range.)
rand('seed', 21); randn('seed', 21);
nl = 15; nq = 4; ns = 2; m = 20;
n = nl + nq + ns^2;
K = struct('f', 0, 'l', nl, 'q', nq, 'r', zeros(1,0), 's', ns);

At = full(sprand(n, m, 0.06));
At(nl+1:nl+nq+ns^2, :) = full(sprand(nq+ns^2, m, 0.02)); % extra-sparse Lorentz+PSD rows
At(2, :) = 1 + rand(1, m);          % dense LP variable
At(nl+1:nl+nq, :) = 1 + rand(nq, m); % dense Lorentz block (trace + norm-bound)

xq = randn(nq-1,1);
Xfeas = [2 0.3; 0.3 2];             % symmetric, strictly PD
xfeas = [rand(nl,1)+0.1; norm(xq)+0.5; xq; Xfeas(:)];
b = At' * xfeas;
c = rand(n,1) - 0.5;

fprintf('nnz(dense LP row) = %d/%d, nnz(dense Lorentz rows) = %d/%d\n', ...
    nnz(At(2,:)), m, nnz(At(nl+1:nl+nq,:)), nq*m);
save_case(out_dir, 'lp_socp_sdp_dense_feasible', At, b, c, K, pars);

fprintf('sedumi dense-columns oracle written to %s\n', out_dir);

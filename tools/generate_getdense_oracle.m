% generate_getdense_oracle.m
%
% Oracle for this port's getdense.py, the dense-column *detection*
% heuristic used by the dense-columns optimization. Follows
% generate_incorder_oracle.m's style of hand-building At/K directly in
% already-normalized (post-pretransfo) internal form -- mainblks/
% qblkstart/sblkstart/lq built by hand per pretransfo.py's own
% construction formulas -- rather than running the full pretransfo
% pipeline, so injected dense rows land at exactly the intended internal
% row positions. Exercises: a plain LP dense-column case, a dense-
% Lorentz-block case (dense.q), a mixed LP+Lorentz+PSD case (also checks
% the PSD-driven "h" floor), a case with nothing flagged dense, and the
% ">m/2" safety-valve reset case.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_getdense_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'getdense');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

function S = run_case(At, K, denq, denf)
    m = size(At, 2);
    Ablkjc = partitA(At, K.mainblks);
    pars.denq = denq;
    pars.denf = denf;
    [dense, Adotdden] = getdense(At, Ablkjc, K, pars);

    S.A = At;
    S.mainblks = K.mainblks;
    S.qblkstart = K.qblkstart;
    S.sblkstart = K.sblkstart;
    S.lq = K.lq;
    S.Kl = K.l;
    S.Kq = K.q;
    S.Ks = K.s;
    S.Ablkjc = Ablkjc;
    S.denq = denq;
    S.denf = denf;
    S.dense_l = dense.l;
    S.dense_cols = dense.cols;
    S.dense_q = dense.q;
    S.Adotdden = Adotdden;
    fprintf('  dense.l=%d dense.cols=%d dense.q=%d (m=%d)\n', ...
        dense.l, numel(dense.cols), numel(dense.q), m);
end

rand('seed', 7);

%% Case 1: LP-only, one artificially dense row (variable).
K1.l = 30; K1.q = []; K1.s = [];
K1.mainblks = [31; 31; 31];
K1.qblkstart = 31;
K1.sblkstart = 31;
K1.lq = 30;
m1 = 40;
At1 = sprandn(30, m1, 0.05);
At1(3, :) = 1 + rand(1, m1);            % dense LP variable
fprintf('Case 1 (LP dense column):\n');
S1 = run_case(At1, K1, 0.75, 3);
save('-v7', fullfile(out_dir, 'case1.mat'), '-struct', 'S1');

%% Case 2: LP + Lorentz, one whole Lorentz block (trace + norm-bound
%% rows) made dense.
K2.l = 20; K2.q = [5;5;5;5]; K2.s = [];
K2.mainblks = [21; 25; 41];
K2.qblkstart = [25; 29; 33; 37; 41];
K2.sblkstart = 41;
K2.lq = 40;
m2 = 40;
At2 = sprandn(40, m2, 0.05);
At2([21 25 26 27 28], :) = 1 + rand(5, m2);   % dense first Lorentz block
fprintf('Case 2 (dense Lorentz block):\n');
S2 = run_case(At2, K2, 0.75, 3);
save('-v7', fullfile(out_dir, 'case2.mat'), '-struct', 'S2');

%% Case 3: LP + Lorentz + PSD, dense LP row + dense Lorentz block together.
K3.l = 15; K3.q = [4;4;4]; K3.s = [3;2];
K3.mainblks = [16; 19; 28];
K3.qblkstart = [19; 22; 25; 28];
K3.sblkstart = [28; 37; 41];
K3.lq = 27;
m3 = 45;
At3 = [sprandn(27, m3, 0.05); sprandn(13, m3, 0.01)];  % lower PSD density keeps h near its floor
At3(2, :) = 1 + rand(1, m3);                  % dense LP variable
At3([16 19 20 21], :) = 1 + rand(4, m3);      % dense first Lorentz block
fprintf('Case 3 (mixed LP+Lorentz+PSD):\n');
S3 = run_case(At3, K3, 0.75, 3);
save('-v7', fullfile(out_dir, 'case3.mat'), '-struct', 'S3');

%% Case 4: nothing flagged dense (uniform sparsity, generous threshold).
K4.l = 20; K4.q = [4;4]; K4.s = 3;
K4.mainblks = [21; 23; 29];
K4.qblkstart = [23; 26; 29];
K4.sblkstart = [29; 38];
K4.lq = 28;
m4 = 30;
At4 = sprandn(37, m4, 0.1);
fprintf('Case 4 (nothing dense):\n');
S4 = run_case(At4, K4, 0.75, 10);
save('-v7', fullfile(out_dir, 'case4.mat'), '-struct', 'S4');

%% Case 5: safety valve -- low thresholds flag more than m/2 columns dense.
K5.l = 30; K5.q = []; K5.s = [];
K5.mainblks = [31; 31; 31];
K5.qblkstart = 31;
K5.sblkstart = 31;
K5.lq = 30;
m5 = 10;
At5 = sprandn(30, m5, 0.5);
fprintf('Case 5 (safety valve, low thresholds):\n');
S5 = run_case(At5, K5, 0.01, 1.01);
save('-v7', fullfile(out_dir, 'case5.mat'), '-struct', 'S5');

fprintf('getdense oracle written to %s\n', out_dir);

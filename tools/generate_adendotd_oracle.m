% generate_adendotd_oracle.m
%
% Oracle for this port's _native.adendotd (dense-columns optimization,
% Stage 2 of the plan). Calls the real vendored adendotd.mex directly
% with small hand-built dense/d/sparAd/Ablk/blkstart structs -- a pure,
% isolated kernel test, no pretransfo/solve pipeline needed.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_adendotd_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'adendotd');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

function S = run_case(dense, d, sparAd, Ablk, qblkstart)
    Ad = adendotd(dense, d, sparAd, Ablk, qblkstart);
    S.dense_l = dense.l;
    S.dense_q = dense.q;
    S.dense_cols = dense.cols;
    S.dense_A = dense.A;
    S.d_q1 = d.q1;
    S.d_q2 = d.q2;
    S.sparAd = sparAd;
    S.Ablk = Ablk;
    S.qblkstart = qblkstart;
    S.Ad = Ad;
end

rand('seed', 13);

%% Case 1: single dense Lorentz block (nq=1), one dense norm-bound column,
%% m=5 constraints, full Ablk pattern.
m1 = 5;
dense1.l = 0;
dense1.q = 1;
dense1.cols = [999; 6];  % [trace placeholder, dense normbound row (1-idx)]
dense1.A = sparse(m1, 2);
dense1.A(:,1) = rand(m1,1);   % trace column
dense1.A([1 3],2) = rand(2,1); % normbound column, sparse
d1.q1 = 3.2;
d1.q2 = [1.5; 2.5; 3.5];  % covers global rows firstQ..firstQ+2 = 4,5,6 (qblkstart=[5,8])
qblkstart1 = [5; 8];
sparAd1 = sparse(m1, 1);
sparAd1(2) = 4.4;
Ablk1 = sparse(ones(m1,1));  % full pattern
S1 = run_case(dense1, d1, sparAd1, Ablk1, qblkstart1);
save('-v7', fullfile(out_dir, 'case1.mat'), '-struct', 'S1');
fprintf('case1 Ad: %s\n', mat2str(full(S1.Ad)'));

%% Case 2: two dense Lorentz blocks (nq=2), multiple dense norm-bound
%% columns split across them, partial Ablk pattern (not every row stored).
m2 = 6;
dense2.l = 0;
dense2.q = [1; 2];
dense2.cols = [901; 902; 7; 9; 14];  % nq=2 placeholders, then 3 normbound rows
dense2.A = sparse(m2, 5);  % nq(2)+nden(3) = 5 columns
dense2.A(:,1) = rand(m2,1);       % trace col block1
dense2.A(:,2) = rand(m2,1);       % trace col block2
dense2.A([1 2],3) = rand(2,1);    % normbound col (block1, row7)
dense2.A([2 4],4) = rand(2,1);    % normbound col (block1, row9)
dense2.A([3 5],5) = rand(2,1);    % normbound col (block2, row14)
d2.q1 = [1.1; 2.2];
d2.q2 = zeros(11,1);
d2.q2([7 9 14] - 5) = [0.7; 0.9; 1.4];  % firstQ = qblkstart(1)-1 = 5
qblkstart2 = [6; 11; 16];
sparAd2 = sparse(m2, 2);
sparAd2([1 3], 1) = [0.3; 0.6];
sparAd2([2 5], 2) = [0.2; 0.8];
Ablk2 = sparse(m2, 2);
Ablk2([1 2 3 4], 1) = 1;   % row 5 (0-indexed) excluded on purpose
Ablk2([2 3 5 6], 2) = 1;
S2 = run_case(dense2, d2, sparAd2, Ablk2, qblkstart2);
save('-v7', fullfile(out_dir, 'case2.mat'), '-struct', 'S2');
fprintf('case2 Ad:\n'); disp(full(S2.Ad));

fprintf('adendotd oracle written to %s\n', out_dir);

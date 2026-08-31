% generate_adenscale_oracle.m
%
% Oracle for this port's _native.adenscale (dense-columns optimization,
% Stage 2 of the plan). Calls the real vendored adenscale.mex directly
% with small hand-built dense/d/qblkstart structs -- no pretransfo/solve
% pipeline needed since this is a pure, isolated kernel.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_adenscale_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'adenscale');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

function S = run_case(dense, d, qblkstart)
    smult = adenscale(dense, d, qblkstart);
    S.dense_l = dense.l;
    S.dense_q = dense.q;
    S.dense_cols = dense.cols;
    S.d_det = d.det;
    S.qblkstart = qblkstart;
    S.smult = smult;
end

rand('seed', 11);

%% Case 1: two dense Lorentz blocks, several dense norm-bound columns each.
dense1.l = 0;
dense1.q = [1; 2];
% cols layout: [trace placeholders (nq=2)][dense norm-bound rows (global,
% 1-indexed, within block1=[5,10), block2=[10,17))]
dense1.cols = [901; 902; 6; 8; 9; 12; 15; 16];
d1.det = [3.5; 7.25; 1.1];
qblkstart1 = [5; 10; 17; 20];
S1 = run_case(dense1, d1, qblkstart1);
save('-v7', fullfile(out_dir, 'case1.mat'), '-struct', 'S1');
fprintf('case1 smult: %s\n', mat2str(S1.smult'));

%% Case 2: single dense Lorentz block, one dense norm-bound column.
dense2.l = 0;
dense2.q = 1;
dense2.cols = [50; 6];
d2.det = 4.0;
qblkstart2 = [5; 10];
S2 = run_case(dense2, d2, qblkstart2);
save('-v7', fullfile(out_dir, 'case2.mat'), '-struct', 'S2');
fprintf('case2 smult: %s\n', mat2str(S2.smult'));

%% Case 3: no dense norm-bound columns at all (nden=0), only LP-dense cols
%% (dense.l>0) plus a Lorentz block with zero norm-bound entries selected.
dense3.l = 2;
dense3.q = 1;
dense3.cols = [1; 2; 77];  % nl=2 placeholders, nq=1 placeholder, nden=0
d3.det = 9.0;
qblkstart3 = [5; 10];
S3 = run_case(dense3, d3, qblkstart3);
save('-v7', fullfile(out_dir, 'case3.mat'), '-struct', 'S3');
fprintf('case3 smult: %s\n', mat2str(S3.smult'));

fprintf('adenscale oracle written to %s\n', out_dir);

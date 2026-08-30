% generate_incorder_oracle.m
%
% Runs the real Octave/MEX incorder.c (greedy incremental column
% ordering used by getada3 for the PSD dense-column-elimination
% pivoting order) across a few sparsity patterns.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_incorder_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'incorder');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

rand('seed', 77);

%% Case 1: default args (whole matrix, ifirst=1), moderate sparsity.
At1 = sprandn(8, 6, 0.4);
[perm1, dz1] = incorder(At1);
save('-v7', fullfile(out_dir, 'case1.mat'), 'At1', 'perm1', 'dz1');

%% Case 2: explicit Ajc1/ifirst restricting to a PSD sub-range, mirroring
%% sedumi.m's real call `incorder(A, Ablkjc(:,3), K.mainblks(3))`.
At2 = sprandn(10, 5, 0.5);
blkstart2 = [1; 4; 11];             % 2 row-blocks: [1,3], [4,10]
Ablkjc2 = partitA(At2, blkstart2);
ifirst2 = blkstart2(2);             % first "PSD" row = 4
[perm2, dz2] = incorder(At2, Ablkjc2(:,2), ifirst2);
save('-v7', fullfile(out_dir, 'case2.mat'), 'At2', 'Ablkjc2', 'ifirst2', 'perm2', 'dz2');

%% Case 3: several tied (equal nnz-count) columns, to check tie-break
%% order matches the greedy "first encountered" rule.
At3 = sparse(6, 4);
At3(1,1) = 1; At3(2,1) = 1;
At3(1,2) = 1; At3(3,2) = 1;
At3(4,3) = 1; At3(5,3) = 1;
At3(2,4) = 1; At3(4,4) = 1;
[perm3, dz3] = incorder(At3);
save('-v7', fullfile(out_dir, 'case3.mat'), 'At3', 'perm3', 'dz3');

%% Case 4: a column with zero PSD nonzeros (degenerate "deg" column).
At4 = sprandn(7, 5, 0.4);
At4(:,3) = 0; %#ok<SPRIX>
[perm4, dz4] = incorder(At4);
save('-v7', fullfile(out_dir, 'case4.mat'), 'At4', 'perm4', 'dz4');

fprintf('incorder oracle written to %s\n', out_dir);

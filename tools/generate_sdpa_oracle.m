% generate_sdpa_oracle.m
%
% Oracle for sdpa.py's read_sdpa() (Phase 4): runs the real
% conversion/fromsdpa.m against a small, hand-crafted sparse SDPA file
% (tests/fixtures/sdpa/test_problem.dat-s -- 2 constraints, one -2
% diagonal/LP block plus one 2x2 PSD block) and saves the resulting
% (At,b,c,K) as the oracle test_sdpa.py compares sedumipy.sdpa.read_sdpa's
% own parse against.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_sdpa_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
addpath(fullfile(vendor_root, 'conversion'));

fixture_dir = fullfile(repo_root, 'tests', 'fixtures', 'sdpa');
dat_file = fullfile(fixture_dir, 'test_problem.dat-s');

[At, b, c, K] = fromsdpa(dat_file);
save('-mat', fullfile(fixture_dir, 'test_problem_oracle.mat'), 'At', 'b', 'c', 'K');
printf('wrote %s\n', fullfile(fixture_dir, 'test_problem_oracle.mat'));

% search_blkchol_diagadd.m
%
% Same idea as search_blkchol_edge_cases.m, over the more targeted batch
% in tests/fixtures/blkchol_diagadd_candidates/*.mat -- looking
% specifically for nadd>0.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; search_blkchol_diagadd"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
cand_dir = fullfile(repo_root, 'tests', 'fixtures', 'blkchol_diagadd_candidates');

addpath(vendor_root);
global ADA_sedumi_

files = dir(fullfile(cand_dir, '*.mat'));
nfound = 0;
for i = 1:numel(files)
  name = files(i).name;
  data = load(fullfile(cand_dir, name), 'X');
  ADA_sedumi_ = data.X;
  try
    L = symbchol();
    [LL, Ld, Lskip, Ladd] = blkchol(L, data.X);
    nskip = nnz(Lskip);
    nadd = nnz(Ladd);
    if nadd > 0
      fprintf('%-40s nskip=%3d nadd=%3d  <-- FOUND\n', name, nskip, nadd);
      nfound = nfound + 1;
    end
  catch err
    fprintf('%-40s ERROR: %s\n', name, err.message);
  end
end
fprintf('\n%d/%d candidates with nadd>0\n', nfound, numel(files));

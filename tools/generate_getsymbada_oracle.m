% generate_getsymbada_oracle.m
%
% Runs the real Octave/MEX getsymbada.m (+ its partitA/findblks/extractA
% dependencies) exactly as sedumi.m's own pre-main-loop setup does
% (see sedumi.m lines ~356-382), for a handful of small LP+Lorentz+PSD
% problems, to produce oracle fixtures for this port's getsymbada.py.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_getsymbada_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'getsymbada');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

function save_case(out_dir, name, A, K, Ablkjc, q_pattern, SYMBADA)
    mainblks = K.mainblks; qblkstart = K.qblkstart; sblkstart = K.sblkstart; %#ok<NASGU>
    save('-v7', fullfile(out_dir, [name '.mat']), ...
         'A', 'mainblks', 'qblkstart', 'sblkstart', 'Ablkjc', 'q_pattern', 'SYMBADA');
end

function [Ablkjc, q_pattern, SYMBADA] = run_case(A, K)
    Ablkjc = partitA(A, K.mainblks);
    q_pattern = findblks(A, Ablkjc, 2, 3, K.qblkstart);
    if ~isempty(q_pattern)
        q_pattern = q_pattern + spones(extractA(A, Ablkjc, 1, 2, K.mainblks(1), K.mainblks(2)));
    end
    SYMBADA = getsymbada(A, Ablkjc, struct('q', q_pattern), K.sblkstart);
end

rand('seed', 88);

%% Case 1: sparse LP + one Lorentz cone + one small PSD block.
K1.l = 3; K1.q = 4; K1.s = 3;
n1 = K1.l + K1.q + K1.s^2;
m1 = 5;
A1 = sprandn(m1, n1, 0.3);
[x1, y1, K1out, prep1] = deal([]); %#ok<ASGLU>
[Anorm1, ~, ~, K1n] = pretransfo(A1, ones(m1,1), ones(n1,1), K1, struct());
[Ablkjc1, q_pattern1, SYMBADA1] = run_case(Anorm1, K1n);
save_case(out_dir, 'case1', Anorm1, K1n, Ablkjc1, q_pattern1, SYMBADA1);

%% Case 2: sparse LP + two Lorentz cones + two PSD blocks, denser A to
%% approach (but not necessarily hit) the 0.9 density fallback.
K2.l = 2; K2.q = [3;3]; K2.s = [2;3];
n2 = K2.l + sum(K2.q) + sum(K2.s.^2);
m2 = 6;
A2 = sprandn(m2, n2, 0.5);
[Anorm2, ~, ~, K2n] = pretransfo(A2, ones(m2,1), ones(n2,1), K2, struct());
[Ablkjc2, q_pattern2, SYMBADA2] = run_case(Anorm2, K2n);
save_case(out_dir, 'case2', Anorm2, K2n, Ablkjc2, q_pattern2, SYMBADA2);

%% Case 3: fully dense A, to trigger the spars(...)==1 dense fallback.
K3.l = 2; K3.q = 3; K3.s = 2;
n3 = K3.l + K3.q + K3.s^2;
m3 = 4;
A3 = rand(m3, n3);
[Anorm3, ~, ~, K3n] = pretransfo(A3, ones(m3,1), ones(n3,1), K3, struct());
[Ablkjc3, q_pattern3, SYMBADA3] = run_case(Anorm3, K3n);
save_case(out_dir, 'case3', Anorm3, K3n, Ablkjc3, q_pattern3, SYMBADA3);

%% Case 4: no Lorentz cones at all (q_pattern has 0 rows), PSD-only +
%% LP, to exercise the isempty(DAt.q) branch.
K4.l = 3; K4.s = [2;2];
n4 = K4.l + sum(K4.s.^2);
m4 = 5;
A4 = sprandn(m4, n4, 0.35);
[Anorm4, ~, ~, K4n] = pretransfo(A4, ones(m4,1), ones(n4,1), K4, struct());
[Ablkjc4, q_pattern4, SYMBADA4] = run_case(Anorm4, K4n);
save_case(out_dir, 'case4', Anorm4, K4n, Ablkjc4, q_pattern4, SYMBADA4);

fprintf('getsymbada oracle written to %s\n', out_dir);

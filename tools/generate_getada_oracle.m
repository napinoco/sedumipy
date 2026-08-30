% generate_getada_oracle.m
%
% Runs the real Octave/MEX getada1/getada2/getada3 (getada{1,2,3}.c)
% through the exact pre-loop Aord/DAt setup and main-loop call sequence
% sedumi.m uses (see sedumi.m lines ~356-378 and ~450-452), for a
% handful of small LP+Lorentz+PSD problems, to produce oracle fixtures
% for this port's _native.getada1/getada2/getada3 ctypes bindings.
%
% `d.l`, `d.det`, DAt.q's numeric values, and `udsqr` are NOT taken from
% a real sdinit()/invcholfac() solve state -- getada1/2/3 are agnostic
% to where their scaling inputs come from (they just consume arrays of
% the documented sizes), so this fabricates plausible positive random
% values of the right sizes directly. This keeps the oracle focused on
% the kernels' own sparsity/indexing logic, which is exactly what this
% port's bindings need to match; the full real-`d` integration is
% covered end-to-end by the Step 7 solve-level fixtures instead.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_getada_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'getada');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end
addpath(vendor_root);

function save_case(out_dir, name, S)
    save('-v7', fullfile(out_dir, [name '.mat']), '-struct', 'S');
end

function S = run_case(A, K)
    m = size(A, 2);
    Ablkjc = partitA(A, K.mainblks);
    Aord.lqperm = sortnnz(A, [], Ablkjc(:,3));
    DAt.q = findblks(A, Ablkjc, 2, 3, K.qblkstart);
    if ~isempty(DAt.q)
        DAt.q = DAt.q + spones(extractA(A, Ablkjc, 1, 2, K.mainblks(1), K.mainblks(2)));
        [qi, qj, qv] = find(DAt.q); %#ok<ASGLU>
        DAt.q = sparse(qi, qj, rand(size(qv)), size(DAt.q,1), size(DAt.q,2));
        Aord.qperm = sortnnz(DAt.q, [], []);
    else
        DAt.q = sparse(0, m);
        Aord.qperm = (1:m)';
    end
    [Aord.sperm, Aord.dz] = incorder(A, Ablkjc(:,3), K.mainblks(3));

    ADA0 = getsymbada(A, Ablkjc, DAt, K.sblkstart);

    d.l = rand(K.l, 1) + 0.1;
    d.det = rand(length(K.q), 1) + 0.1;

    lenud = K.rLen + 2*K.hLen; %#ok<NASGU>
    % rLen/hLen not always present on K from pretransfo unless SDP present;
    % fall back to summing K.s^2 directly (real case only in these fixtures).
    lenud = sum(K.s .^ 2);
    udsqr = rand(lenud, 1) + 0.1;

    ADA1 = getada1(ADA0, A, Ablkjc(:,3), Aord.lqperm, d, K.qblkstart);
    ADA2 = getada2(ADA1, DAt, Aord, K);
    [ADA3, absd] = getada3(ADA2, A, Ablkjc(:,3), Aord, udsqr, K);

    S.A = A;
    S.mainblks = K.mainblks;
    S.qblkstart = K.qblkstart;
    S.sblkstart = K.sblkstart;
    S.blkstart = K.blkstart;
    S.Kl = K.l;
    S.Kq = K.q;
    S.Ks = K.s;
    S.Ablkjc = Ablkjc;
    S.lqperm = Aord.lqperm;
    S.qperm = Aord.qperm;
    S.sperm = Aord.sperm;
    S.dz = Aord.dz;
    S.DAt_q = DAt.q;
    S.d_l = d.l;
    S.d_det = d.det;
    S.udsqr = udsqr;
    S.ADA0 = ADA0;
    S.ADA1 = ADA1;
    S.ADA2 = ADA2;
    S.ADA3 = ADA3;
    S.absd = absd;
end

rand('seed', 77);

%% Case 1: LP + one Lorentz cone + one small PSD block.
K1.l = 3; K1.q = 4; K1.s = 3;
n1 = K1.l + K1.q + K1.s^2;
m1 = 5;
A1 = sprandn(m1, n1, 0.3);
[Anorm1, ~, ~, K1n] = pretransfo(A1, ones(m1,1), ones(n1,1), K1, struct());
S1 = run_case(Anorm1, K1n);
save_case(out_dir, 'case1', S1);

%% Case 2: LP + two Lorentz cones + two PSD blocks, denser A.
K2.l = 2; K2.q = [3;3]; K2.s = [2;3];
n2 = K2.l + sum(K2.q) + sum(K2.s.^2);
m2 = 6;
A2 = sprandn(m2, n2, 0.5);
[Anorm2, ~, ~, K2n] = pretransfo(A2, ones(m2,1), ones(n2,1), K2, struct());
S2 = run_case(Anorm2, K2n);
save_case(out_dir, 'case2', S2);

%% Case 3: no Lorentz cones at all, PSD-only + LP (exercises getada2's
%% "no Lorentz blocks" short-circuit).
K3.l = 3; K3.s = [2;2];
n3 = K3.l + sum(K3.s.^2);
m3 = 5;
A3 = sprandn(m3, n3, 0.35);
[Anorm3, ~, ~, K3n] = pretransfo(A3, ones(m3,1), ones(n3,1), K3, struct());
S3 = run_case(Anorm3, K3n);
save_case(out_dir, 'case3', S3);

%% Case 4: no PSD blocks at all, LP + Lorentz only (exercises getada3's
%% "no PSD blocks" cpspdiag() fallback for absd).
K4.l = 3; K4.q = [4;3];
n4 = K4.l + sum(K4.q);
m4 = 6;
A4 = sprandn(m4, n4, 0.4);
[Anorm4, ~, ~, K4n] = pretransfo(A4, ones(m4,1), ones(n4,1), K4, struct());
S4 = run_case(Anorm4, K4n);
save_case(out_dir, 'case4', S4);

fprintf('getada oracle written to %s\n', out_dir);

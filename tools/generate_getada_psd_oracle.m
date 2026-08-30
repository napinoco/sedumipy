% generate_getada_psd_oracle.m
%
% Oracle for this port's build_aord()/getada_psd() (getada_psd.py),
% covering the full one-time pre-loop Aord/getsymbada setup AND one
% main-loop-style getada1->getada2->getada3 iteration, exactly as
% sedumi.m's own pre-loop code (lines ~330-382) and main loop (lines
% ~450-452) run them -- including real invcholfac() (not a fabricated
% udsqr, unlike tools/generate_getada_oracle.m's narrower per-kernel
% fixtures), to validate the full Step 4/5 orchestration end-to-end.
%
% Run from the repository root:
%   octave-cli --no-gui --eval "cd tools; generate_getada_psd_oracle"

repo_root = fileparts(fileparts(mfilename('fullpath')));
vendor_root = fullfile(repo_root, 'vendor', 'sedumi-upstream');
out_dir = fullfile(repo_root, 'tests', 'fixtures', 'getada_psd');
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
        q_pattern = DAt.q + spones(extractA(A, Ablkjc, 1, 2, K.mainblks(1), K.mainblks(2)));
        Aord.qperm = sortnnz(q_pattern, [], []);
        [qi, qj, qv] = find(DAt.q); %#ok<ASGLU>
        DAt.q = sparse(qi, qj, rand(size(qv)), size(DAt.q,1), size(DAt.q,2));
    else
        q_pattern = sparse(0, m);
        Aord.qperm = (1:m)';
    end
    [Aord.sperm, Aord.dz] = incorder(A, Ablkjc(:,3), K.mainblks(3));

    ADA0 = getsymbada(A, Ablkjc, struct('q', q_pattern), K.sblkstart);

    d.l = rand(K.l, 1) + 0.1;
    d.det = rand(length(K.q), 1) + 0.1;
    lenud = sum(K.s .^ 2);
    d.u = rand(lenud, 1) + 0.1;
    d.perm = [];

    udsqr = invcholfac(d.u, K, d.perm);

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
    S.q_pattern = q_pattern;
    S.SYMBADA = ADA0;
    S.DAt_q = DAt.q;
    S.d_l = d.l;
    S.d_det = d.det;
    S.d_u = d.u;
    S.udsqr = udsqr;
    S.ADA1 = ADA1;
    S.ADA2 = ADA2;
    S.ADA3 = ADA3;
    S.absd = absd;
end

rand('seed', 99);

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

%% Case 3: no Lorentz cones at all, PSD-only + LP.
K3.l = 3; K3.s = [2;2];
n3 = K3.l + sum(K3.s.^2);
m3 = 5;
A3 = sprandn(m3, n3, 0.35);
[Anorm3, ~, ~, K3n] = pretransfo(A3, ones(m3,1), ones(n3,1), K3, struct());
S3 = run_case(Anorm3, K3n);
save_case(out_dir, 'case3', S3);

fprintf('getada_psd oracle written to %s\n', out_dir);

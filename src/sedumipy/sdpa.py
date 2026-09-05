"""Port of conversion/fromsdpa.m: reads an SDP problem from a sparse
SDPA-formatted file (the format SDPLIB and most other public SDP
benchmark collections use) into SeDuMi's own (At, b, c, K) form.

::

    At, b, c, K = read_sdpa("arch0.dat-s")
    x, y, info = sedumipy.sedumi(At, b, c, K)

SDPA sparse format, in order: mDIM (m, the number of constraints),
nBLOCK (the number of blocks), a line of nBLOCK block sizes (negative
size n means n diagonal/LP variables, not a PSD block), the b vector (m
entries), then one line per nonzero coefficient as
"<matno> <blkno> <i> <j> <value>" -- matno=0 is the cost matrix C,
matno=1..m is constraint matrix A_matno; only the upper or lower triangle
of each symmetric block needs to be given (both fromsdpa.m and this port
mirror whichever triangle is present across the diagonal). SDPA's own
convention is to MAXIMIZE trace(C,X); SeDuMi minimizes, so C's sign is
flipped on the way in (and back on the way out, in write_sdpa).

write_sdpa (the inverse direction) is new functionality, not a port --
original SeDuMi's own conversion/writesdp.m instead writes the unrelated
SDPpack format, not SDPA sparse. It's the direct inverse of the reading
convention above, verified by round-tripping through read_sdpa itself
(tests/test_sdpa.py) as well as against a real Octave fromsdpa.m oracle.
"""

from __future__ import annotations

import re

import numpy as np
import scipy.sparse as sp


def read_sdpa(fname):
    """(At, b, c, K) = read_sdpa(fname): see module docstring. `fname` is
    a path to an uncompressed SDPA sparse (.dat-s) file."""
    with open(fname) as f:
        lines = f.readlines()

    idx = 0

    def next_line():
        nonlocal idx
        line = lines[idx]
        idx += 1
        return line

    # ---- mDIM: skip comment lines (fromsdpa.m just skips any line whose
    # leading token isn't parseable as the integer m; SDPA comment lines
    # conventionally start with '"' or '*', but real files are lenient) ----
    m = None
    while m is None:
        line = next_line()
        match = re.match(r"\s*(\d+)", line)
        if match:
            m = int(match.group(1))
    nblocks = int(re.match(r"\s*(\d+)", next_line()).group(1))

    # ---- block sizes: '.,(){}' are treated as whitespace, like fromsdpa.m ----
    dims_line = re.sub(r"[.,(){}]", " ", next_line())
    dims = np.array([int(tok) for tok in dims_line.split()][:nblocks], dtype=np.int64)
    if dims.size != nblocks or np.any(dims == 0):
        raise ValueError(f"{fname}: invalid SDPA block dimensions")

    N = int(-dims[dims < 0].sum() + (dims[dims > 0] ** 2).sum())

    # ---- offset[k] = 0-indexed flat position where block k+1 starts;
    # diagonal/1x1 blocks (dims<=1) come first, contiguously, then PSD
    # blocks (dims>1) -- matches fromsdpa.m's loffset/sdpoffset exactly ----
    offset = np.zeros(nblocks, dtype=np.int64)
    loffset = 0
    sdpoffset = int(np.abs(dims[dims <= 1]).sum())
    for k in range(nblocks):
        if dims[k] <= 1:
            offset[k] = loffset
            loffset += abs(int(dims[k]))
        else:
            offset[k] = sdpoffset
            sdpoffset += int(dims[k]) ** 2
    stride = np.where(dims < 0, 0, dims)

    # ---- b vector ----
    b_tokens = []
    while len(b_tokens) < m:
        b_tokens.extend(re.sub(r"[,(){}]", " ", next_line()).split())
    b = np.array([float(tok) for tok in b_tokens[:m]], dtype=np.float64)

    # ---- coefficients: matno blkno i j value, one nonzero per line ----
    matno, blkno, ii, jj, val = [], [], [], [], []
    for line in lines[idx:]:
        tokens = line.split()
        if not tokens:
            continue
        matno.append(int(tokens[0]))
        blkno.append(int(tokens[1]))
        ii.append(int(tokens[2]))
        jj.append(int(tokens[3]))
        val.append(float(tokens[4]))
    matno = np.array(matno, dtype=np.int64)
    blkno = np.array(blkno, dtype=np.int64) - 1  # 0-index into dims/offset/stride
    ii = np.array(ii, dtype=np.int64)
    jj = np.array(jj, dtype=np.int64)
    val = np.array(val, dtype=np.float64)

    def flat_positions(sel):
        """0-indexed flat row position of (i,j) and its symmetric (j,i)
        counterpart, for the selected entries -- fromsdpa.m's
        offset(blk)+(i-1)*stride(blk)+j, 0-indexed."""
        b_ = blkno[sel]
        i_ = ii[sel]
        j_ = jj[sel]
        pos_ij = offset[b_] + (i_ - 1) * stride[b_] + (j_ - 1)
        pos_ji = offset[b_] + (j_ - 1) * stride[b_] + (i_ - 1)
        return pos_ij, pos_ji

    c_sel = matno == 0
    pos_ij, pos_ji = flat_positions(c_sel)
    v = val[c_sel]
    on_diag = ii[c_sel] == jj[c_sel]
    rows = np.concatenate([pos_ij, pos_ji])
    vals = np.concatenate([v, np.where(on_diag, 0.0, v)])
    c = sp.coo_matrix((-vals, (rows, np.zeros_like(rows))), shape=(N, 1)).tocsc()
    c.eliminate_zeros()

    a_sel = matno != 0
    pos_ij, pos_ji = flat_positions(a_sel)
    v = val[a_sel]
    on_diag = ii[a_sel] == jj[a_sel]
    cols = matno[a_sel] - 1
    rows = np.concatenate([pos_ij, pos_ji])
    vals = np.concatenate([v, np.where(on_diag, 0.0, v)])
    cols = np.concatenate([cols, cols])
    At = sp.coo_matrix((vals, (rows, cols)), shape=(N, m)).tocsc()
    At.eliminate_zeros()

    K = {
        "l": int(-dims[dims < 0].sum() + (dims == 1).sum()),
        "s": dims[dims > 1],
    }
    return At, b, c, K


def write_sdpa(fname, At, b, c, K, comment="sedumipy") -> None:
    """write_sdpa(fname, At, b, c, K): writes (At, b, c, K) out as a
    sparse SDPA (.dat-s) file -- the inverse of read_sdpa. LP-only (K.l)
    problems are written as a single diagonal block; K.q/K.r (SDPA has no
    second-order-cone block type) are rejected. Only the upper triangle
    of each PSD block is written, matching common SDPA sparse-format
    practice (see module docstring)."""
    Kl = int(K.get("l", 0) or 0)
    Ks = np.atleast_1d(np.asarray(K.get("s", []), dtype=np.int64)).ravel()
    if len(np.atleast_1d(np.asarray(K.get("q", [])))) or len(
        np.atleast_1d(np.asarray(K.get("r", [])))
    ):
        raise ValueError("write_sdpa: K.q/K.r (second-order cones) are not representable in SDPA format")

    At = At.tocsc() if sp.issparse(At) else sp.csc_matrix(At)
    c = np.asarray(c.todense() if sp.issparse(c) else c).ravel()
    b = np.asarray(b.todense() if sp.issparse(b) else b).ravel()
    m = At.shape[1]

    dims = ([-Kl] if Kl else []) + [int(s) for s in Ks]
    nblocks = len(dims)
    offset = []
    pos = 0
    for d in dims:
        offset.append(pos)
        pos += abs(d) if d < 0 else d * d
    stride = [0 if d < 0 else d for d in dims]

    def block_entries(vec):
        """yield (blkno, i, j, value) for vec's upper triangle, one
        entry per structurally-nonzero position, 1-indexed i<=j."""
        for k, d in enumerate(dims):
            lo, hi = offset[k], offset[k] + (abs(d) if d < 0 else d * d)
            block = vec[lo:hi]
            if d < 0:
                for i, value in enumerate(block):
                    if value != 0.0:
                        yield k + 1, i + 1, i + 1, value
            else:
                mat = block.reshape(d, d)
                for i in range(d):
                    for j in range(i, d):
                        value = mat[i, j]
                        if value != 0.0:
                            yield k + 1, i + 1, j + 1, value

    with open(fname, "w") as f:
        f.write(f'"{comment}"\n')
        f.write(f"{m}\n{nblocks}\n")
        f.write(" ".join(str(d) for d in dims) + "\n")
        f.write(" ".join(f"{v:.18e}" for v in b) + "\n")
        for blkno, i, j, value in block_entries(-c):
            f.write(f"0 {blkno} {i} {j} {value:.18e}\n")
        for col in range(m):
            for blkno, i, j, value in block_entries(np.asarray(At[:, col].todense()).ravel()):
                f.write(f"{col + 1} {blkno} {i} {j} {value:.18e}\n")

"""Phase 4: .mat I/O for SeDuMi-style problem/solution files.

Not a port of any particular upstream .m file -- real SeDuMi problems and
solutions are just plain MATLAB structs/arrays saved with MATLAB's own
`save`, so there is no "reference .m implementation" to port here, unlike
the rest of this package. This module exists so callers don't have to
hand-write scipy.io.loadmat/savemat plus the K-struct field unwrapping
this port's own tests already do ad hoc (see e.g. test_sedumi.py's
`_load_K`); it follows the same conventions those tests settled on.

Field-name conventions matched here (both are common in the wild, e.g.
vendor/sedumi-upstream/examples/*.mat all use "At"):
  - the constraint matrix may be stored as "A" or "At" (SeDuMi's own
    sedumi() accepts either orientation and disambiguates via shape --
    see pretransfo.py's own (At,K) size check -- so this module does too,
    simply forwarding whichever one it finds along with its name);
  - K's fields arrive as a MATLAB struct (a (1,1) structured array once
    loaded via scipy.io.loadmat); scalar fields ("f", "l") are unwrapped
    to a plain Python int, everything else (K.q, K.r, K.s, K.xcomplex,
    ...) to a flat NumPy array, exactly like checkpars.py/pretransfo.py
    already expect;
  - pars, if present, is unwrapped the same way test_sedumi.py's
    `_load_pars` does (scalar-only fields, one MATLAB scalar per pars
    option -- true of every pars field SeDuMi itself defines).
"""

from __future__ import annotations

import numpy as np
import scipy.io
import scipy.sparse as sp


def _unwrap_K(Kmat) -> dict:
    """Kmat: a scipy.io.loadmat (1,1) structured-array K struct. "f"/"l"
    (if present) come back as a plain Python int; every other field
    (K.q, K.r, K.s, K.xcomplex, ...) stays a flat NumPy array even when
    it has just one element (a single Lorentz/PSD block) -- matching
    pretransfo.py's/checkpars.py's own expectations."""
    K = {}
    for fld in Kmat.dtype.names:
        val = Kmat[fld]
        K[fld] = int(val.item()) if (val.size == 1 and fld in ("f", "l")) else val.ravel()
    return K


def _unwrap_pars(parsmat) -> dict:
    """parsmat: a scipy.io.loadmat (1,1) structured-array pars struct.
    Every pars field SeDuMi defines is itself a scalar option."""
    return {fld: parsmat[fld].item() for fld in parsmat.dtype.names}


def read_mat(path):
    """(A, b, c, K, pars) = read_mat(path): loads a SeDuMi-style problem
    .mat file. `A` is returned exactly as stored (sparse or dense,
    whichever orientation the file used) -- sedumi() accepts either.
    `pars` is `{}` if the file has no "pars" field (SeDuMi problem files
    saved from real examples, like vendor/sedumi-upstream/examples/*.mat,
    normally don't)."""
    data = scipy.io.loadmat(path)

    if "At" in data:
        A = data["At"]
    elif "A" in data:
        A = data["A"]
    else:
        raise ValueError(f"{path}: no 'A' or 'At' field found")

    b = np.asarray(data["b"].todense() if sp.issparse(data["b"]) else data["b"]).ravel()
    c = np.asarray(data["c"].todense() if sp.issparse(data["c"]) else data["c"]).ravel()
    K = _unwrap_K(data["K"][0, 0])
    pars = _unwrap_pars(data["pars"][0, 0]) if "pars" in data else {}

    return A, b, c, K, pars


def write_solution_mat(path, x, y, info) -> None:
    """write_solution_mat(path, x, y, info): saves sedumi()'s own return
    values (x, y, info) as a MATLAB struct-compatible .mat file, the
    layout real SeDuMi's own examples/*.mat solution dumps use (a plain
    "x"/"y" vector plus an "info" struct)."""
    scipy.io.savemat(path, {"x": np.asarray(x), "y": np.asarray(y), "info": info})

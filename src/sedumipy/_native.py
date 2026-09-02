"""ctypes bindings over libsedumi.so (the Phase 1 standalone C library).

Per this project's own coding convention (see CONTRIBUTING.md section 9),
every C kernel call in this port goes through this module -- other
modules reach the C library only via `_native.xxx`, never directly. This
is the final binding surface (Phase 2's cluster 1-5 work, completed):
every MEX target that real SeDuMi's `install_sedumi.m` actually builds
has a binding here, wired into the higher-level module that needs it
(getdense.py, getdatm.py, pcg.py, cone.py, updtransfo.py, wregion.py,
sdinit.py, getada_psd.py, symbchol.py, symbcholden.py, ...) -- see each
binding's own docstring for its `.c`/`.m` source and its callers.

Two of `install_sedumi.m`'s MEX targets are deliberately NOT bound here
via ctypes, and live as pure-Python reimplementations elsewhere instead
(`incorder.py`, `neighborhood.py`): `incorder.c` and `iswnbr.c` both feed
a `qsort()` comparator cast to the wrong function-pointer type (undefined
behavior in C), which was confirmed to actually produce build-dependent
orderings -- see `neighborhood.py`'s and `incorder.py`'s own docstrings.
`sortnnz` below has the same issue but stays in this module as a
pure-Python function (no `_lib.sortnnz` call) rather than moving out,
since it's a simple sort with no separate state to manage.

A handful of bindings in this module are intentionally unused by the
rest of this port -- this was confirmed by checking each one against
`install_sedumi.m`'s build list and, where relevant, its actual mex
argument list, not just by grepping this port's own callers:

- `realdot`/`realssqr`/`scalarmul`/`addscalarmul`: generic BLAS-like
  helpers (sdmauxRdot.c/sdmauxScalarmul.c) with no `.m` wrapper or MEX
  target of their own in real SeDuMi -- they're linked into many other
  kernels' object code, not independently callable there either. Bound
  here from the Phase 1 kernel smoke test (kernel_smoke/smoke_test.c);
  this port's own algorithm code just uses NumPy directly for these
  trivial vectorized ops instead of round-tripping through ctypes.
- `blkmul` and `mJdetd`: their `.c` files exist in csrc/ but neither
  `blkmul.c` nor `mJdetd.c` appears anywhere in `install_sedumi.m`'s MEX
  build list -- they are dead code in real upstream SeDuMi itself (no
  `.m` file ever calls them), not a gap in this port.
- `cholsplit`: real `symbchol.m` does call `cholsplit()` and store the
  result as `L.split`, but no other real `.m`/`.c` file ever reads
  `L.split` back (confirmed via `blkchol.c`'s own mex argument list,
  which reads `L.perm`/`L.L`/`L.xsuper`/`L.tmpsiz` but not `L.split`) --
  so it's vestigial in real SeDuMi too, and omitting its call from
  `symbchol.py` doesn't change any result.

Struct layouts (SedumiKRaw, ConeK) are hand-mirrored from blksdp.h and
must be kept in sync with it field-for-field, in declaration order --
ctypes.Structure uses the platform's normal C struct layout rules (the
System V x86_64 ABI on Linux/macOS), so as long as the field list here
matches blksdp.h's, the memory layout matches what the C side expects.
"""

from __future__ import annotations

import ctypes
import os
import shutil
import subprocess
import sys
from pathlib import Path

_PKG_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _PKG_DIR.parent.parent
if sys.platform == "darwin":
    _LIB_PATH = _PKG_DIR / "libsedumi.dylib"
elif sys.platform == "win32":
    _LIB_PATH = _PKG_DIR / "libsedumi.dll"
else:
    _LIB_PATH = _PKG_DIR / "libsedumi.so"


def _ensure_built() -> Path:
    if _LIB_PATH.exists():
        return _LIB_PATH
    build_script = _REPO_ROOT / "tools" / "build_libsedumi.sh"
    if not build_script.exists():
        raise RuntimeError(
            f"{_LIB_PATH} not found and {build_script} is missing; "
            "cannot build libsedumi automatically."
        )
    # tools/build_libsedumi.sh is a bash script; Windows has no shebang
    # support, so it needs `bash` (an MSYS2 MinGW64 shell -- see
    # CONTRIBUTING.md's Windows note) invoked explicitly -- and by its
    # full path, not the bare name "bash": subprocess.run() on Windows
    # goes through CreateProcess(), which searches the System32
    # directory *before* PATH, and every Windows install since 10
    # 1607 ships a `bash.exe` stub there that just prints "Windows
    # Subsystem for Linux has no installed distributions" and exits
    # nonzero -- it would silently shadow MSYS2's real bash.exe even
    # with MSYS2 correctly first on PATH.
    if sys.platform == "win32":
        bash = shutil.which("bash")
        if bash is None:
            raise RuntimeError(
                "libsedumi needs to be built but no `bash` was found on "
                "PATH (expected an MSYS2 MinGW64 install -- see "
                "CONTRIBUTING.md's Windows note)."
            )
        command = [bash, str(build_script), str(_LIB_PATH)]
    else:
        command = [str(build_script), str(_LIB_PATH)]
    subprocess.run(command, check=True, cwd=_REPO_ROOT)
    if not _LIB_PATH.exists():
        raise RuntimeError(f"build_libsedumi.sh ran but did not produce {_LIB_PATH}")
    return _LIB_PATH


_built_lib_path = _ensure_built()

if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
    # libsedumi.dll dynamically depends on libopenblas.dll (and a couple
    # of mingw runtime DLLs). Python 3.8+ changed ctypes.CDLL()'s default
    # search behavior on Windows to LOAD_LIBRARY_SEARCH_DEFAULT_DIRS,
    # which does NOT include the PATH environment variable -- only the
    # loaded DLL's own directory, the application directory, System32,
    # and directories explicitly registered via os.add_dll_directory().
    # Confirmed the hard way: this used to just add _PKG_DIR (where a
    # wheel's `delvewheel repair` bundles those DLLs alongside
    # libsedumi.dll itself), on the assumption a dev/MSYS2 environment's
    # own PATH would cover the rest -- it does not, and ctypes.CDLL()
    # failed with "Could not find module ... or one of its dependencies"
    # in a real dev-install CI run despite gcc/openblas being correctly
    # on PATH.
    #
    # A first attempt registered *every* PATH directory instead -- that
    # backfired: a failed module import is evicted from sys.modules, so
    # pytest's importorskip() re-executes this file's whole module body
    # (including this loop) on every one of the ~47 test files that
    # import sedumipy, without ever calling os.remove_dll_directory().
    # Windows' internal DLL search path table isn't unbounded, and dozens
    # of PATH entries times dozens of retries overflowed it -- confirmed
    # in CI: even the single, first, always-valid _PKG_DIR call itself
    # started failing with "WinError 206: filename or extension is too
    # long" (a real limit on the cumulative registered search path, not
    # a literal complaint about that one short path).
    #
    # Target the toolchain directory directly instead, the same way
    # tools/build_libsedumi.sh resolves it (MSYS2_ROOT, defaulting to the
    # standard C:\msys64 install location) -- bounded to two directories
    # regardless of how many times this module gets re-executed.
    os.add_dll_directory(str(_PKG_DIR))
    _msys2_root = os.environ.get("MSYS2_ROOT", r"C:\msys64")
    for _candidate in (
        os.path.join(_msys2_root, "mingw64", "bin"),
        os.path.join(_msys2_root, "usr", "bin"),
    ):
        if os.path.isdir(_candidate):
            try:
                os.add_dll_directory(_candidate)
            except OSError:
                pass

_lib = ctypes.CDLL(str(_built_lib_path))

c_size_t_p = ctypes.POINTER(ctypes.c_size_t)
c_double_p = ctypes.POINTER(ctypes.c_double)


class SedumiKRaw(ctypes.Structure):
    """Mirrors `sedumiKRaw` in blksdp.h -- field order matters."""

    _fields_ = [
        ("f", ctypes.c_double),
        ("l", ctypes.c_double),
        ("q", c_double_p),
        ("qN", ctypes.c_size_t),
        ("r", c_double_p),
        ("rN", ctypes.c_size_t),
        ("s", c_double_p),
        ("sN", ctypes.c_size_t),
        ("rsdpNgiven", ctypes.c_char),
        ("rsdpN", ctypes.c_double),
        ("statsGiven", ctypes.c_char),
        ("rLen", ctypes.c_double),
        ("hLen", ctypes.c_double),
        ("qMaxn", ctypes.c_double),
        ("rMaxn", ctypes.c_double),
        ("hMaxn", ctypes.c_double),
        ("blkstart", c_double_p),
        ("blkstartN", ctypes.c_size_t),
    ]


class ConeK(ctypes.Structure):
    """Mirrors `coneK` in blksdp.h -- field order matters."""

    _fields_ = [
        ("frN", ctypes.c_size_t),
        ("lpN", ctypes.c_size_t),
        ("lorN", ctypes.c_size_t),
        ("rconeN", ctypes.c_size_t),
        ("sdpN", ctypes.c_size_t),
        ("rsdpN", ctypes.c_size_t),
        ("qMaxn", ctypes.c_size_t),
        ("rMaxn", ctypes.c_size_t),
        ("hMaxn", ctypes.c_size_t),
        ("rLen", ctypes.c_size_t),
        ("hLen", ctypes.c_size_t),
        ("qDim", ctypes.c_size_t),
        ("rDim", ctypes.c_size_t),
        ("hDim", ctypes.c_size_t),
        ("lorNL", c_double_p),
        ("rconeNL", c_double_p),
        ("sdpNL", c_double_p),
    ]


_lib.conepars_raw.argtypes = [ctypes.POINTER(SedumiKRaw), ctypes.POINTER(ConeK)]
_lib.conepars_raw.restype = None

_lib.realdot.argtypes = [c_double_p, c_double_p, ctypes.c_size_t]
_lib.realdot.restype = ctypes.c_double

_lib.realssqr.argtypes = [c_double_p, ctypes.c_size_t]
_lib.realssqr.restype = ctypes.c_double

_lib.scalarmul.argtypes = [c_double_p, ctypes.c_double, c_double_p, ctypes.c_size_t]
_lib.scalarmul.restype = None

_lib.addscalarmul.argtypes = [c_double_p, ctypes.c_double, c_double_p, ctypes.c_size_t]
_lib.addscalarmul.restype = None

_lib.fwsolve.argtypes = [
    c_double_p,
    c_size_t_p,
    c_size_t_p,
    c_double_p,
    c_size_t_p,
    ctypes.c_size_t,
    c_double_p,
]
_lib.fwsolve.restype = None

_lib.bwsolve.argtypes = [
    c_double_p,
    c_size_t_p,
    c_size_t_p,
    c_double_p,
    c_size_t_p,
    ctypes.c_size_t,
    c_double_p,
]
_lib.bwsolve.restype = None

# ordmmd_'s parameters are all `integer *` in ordmmd.h, and `integer` is
# `mwSignedIndex` (a *signed*, pointer-width int -- ptrdiff_t under
# SEDUMI_STANDALONE), unlike the `mwIndex` (size_t, unsigned) used
# elsewhere in blksdp.h. Using the wrong signedness/width here would be a
# real ABI mismatch, not just a style choice, so this is bound with its
# own pointer type (c_ssize_t) rather than reusing c_size_t_p.
c_ssize_t_p = ctypes.POINTER(ctypes.c_ssize_t)

_lib.ordmmd_.argtypes = [c_ssize_t_p] * 9
_lib.ordmmd_.restype = ctypes.c_int

# sfinit_/symfct_ (symfct.h) use the same `integer` = mwSignedIndex
# convention as ordmmd_.
_lib.sfinit_.argtypes = [c_ssize_t_p] * 15
_lib.sfinit_.restype = ctypes.c_int

_lib.symfct_.argtypes = [c_ssize_t_p] * 17
_lib.symfct_.restype = ctypes.c_int

# expandsub (symfctmex.c) is a plain (non-Fortran-derived) C helper: n and
# nsuper are passed *by value* as mwSize (size_t, unsigned), unlike the
# by-reference mwSignedIndex convention above -- and its pointer params
# are mwIndex* (size_t*, unsigned) even though the buffers they point to
# were just filled by sfinit_/symfct_ through mwSignedIndex* (signed)
# parameters. That mixed signedness is exactly what the original
# mexFunction does too (same buffers, reinterpreted): safe because both
# are pointer-width integers and every value involved is non-negative.
_lib.expandsub.argtypes = [
    ctypes.c_size_t,
    ctypes.c_size_t,
    c_size_t_p,
    c_size_t_p,
    c_size_t_p,
    c_size_t_p,
]
_lib.expandsub.restype = None

_lib.gettmpsiz.argtypes = [c_size_t_p, c_size_t_p, c_size_t_p, ctypes.c_size_t, c_size_t_p]
_lib.gettmpsiz.restype = ctypes.c_size_t

_lib.permuteP.argtypes = [
    c_size_t_p, c_size_t_p, c_double_p,
    c_size_t_p, c_size_t_p, c_double_p,
    c_size_t_p, c_double_p, ctypes.c_size_t,
]
_lib.permuteP.restype = None

_lib.spchol.argtypes = [
    ctypes.c_size_t, ctypes.c_size_t, c_size_t_p,
    c_size_t_p, c_size_t_p, c_size_t_p, c_double_p,
    c_size_t_p, c_double_p, c_double_p, c_size_t_p,
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    c_size_t_p, c_size_t_p,
    ctypes.c_size_t, c_size_t_p, ctypes.c_size_t, c_double_p,
]
_lib.spchol.restype = ctypes.c_size_t

_SIZE_T_ERROR = (1 << 64) - 1  # (mwIndex)-1 sentinel: spchol/blkLDL's
# "insufficient workspace" return value, wrapped around through mwIndex
# being unsigned (size_t) -- ctypes.c_size_t surfaces it as this value
# rather than -1.


def _compute_snode(xsuper, m):
    """Mirrors the small "map each column to its supernode" loop that
    appears identically in choltmpsiz.c and blkchol.c's spchol():
        j = xsuper[0]
        for jsup in range(nsuper): while j < xsuper[jsup+1]: snode[j++] = jsup
    """
    import numpy as np

    nsuper = len(xsuper) - 1
    counts = np.diff(np.asarray(xsuper, dtype=np.int64))
    return np.repeat(np.arange(nsuper, dtype=np.int64), counts).astype(np.uintp)


def _as_double_array(x):
    import numpy as np

    arr = np.ascontiguousarray(x, dtype=np.float64)
    ptr = arr.ctypes.data_as(c_double_p)
    return arr, ptr


def realdot(x, y) -> float:
    """r = sum(x_i * y_i), via the BLAS-backed C kernel (sdmauxRdot.c)."""
    xa, xp = _as_double_array(x)
    ya, yp = _as_double_array(y)
    if xa.shape != ya.shape:
        raise ValueError(f"shape mismatch: {xa.shape} vs {ya.shape}")
    return _lib.realdot(xp, yp, xa.size)


def realssqr(x) -> float:
    """r = sum(x_i^2), via the BLAS-backed C kernel (sdmauxRdot.c)."""
    xa, xp = _as_double_array(x)
    return _lib.realssqr(xp, xa.size)


def scalarmul(alpha: float, x):
    """r = alpha * x, via the BLAS-backed C kernel (sdmauxScalarmul.c)."""
    import numpy as np

    xa, xp = _as_double_array(x)
    out = np.empty_like(xa)
    outp = out.ctypes.data_as(c_double_p)
    _lib.scalarmul(outp, float(alpha), xp, xa.size)
    return out


def addscalarmul(r, alpha: float, x):
    """r += alpha * x, via the BLAS-backed C kernel (sdmauxScalarmul.c)."""
    ra, rp = _as_double_array(r)
    xa, xp = _as_double_array(x)
    _lib.addscalarmul(rp, float(alpha), xp, xa.size)
    return ra


def _as_index_array(x):
    import numpy as np

    arr = np.ascontiguousarray(x, dtype=np.uintp)
    ptr = arr.ctypes.data_as(c_size_t_p)
    return arr, ptr


def _cached_csc_solve_arrays(L_csc):
    """(Ljc, Ljc_p, Lir, Lir_p, Lpr, Lpr_p) for fwsolve()/bwsolve(), memoized
    on the L_csc object itself.

    L_csc.indptr/.indices/.data never change across the many fwsolve()/
    bwsolve() calls PCG makes against the same factorization within one
    outer IPM iteration (only the right-hand side `y` does) -- but scipy
    always normalizes constructed csc_matrix indices to int32/int64
    (verified empirically: passing uintp indices into
    scipy.sparse.csc_matrix() does not stick), so _as_index_array()'s
    uintp cast was previously redone from scratch on every single call.
    For a problem with a large nnz(L) (e.g. DIMACS qssp180old, nnz(L) ~
    8.3e7) that reconversion cost dominates runtime -- cProfile on 5
    outer iterations showed numpy.ascontiguousarray alone at 130s/462s
    (28%) of total time, mostly attributable to exactly this. Caching it
    per L_csc object (a fresh object every outer iteration, so this can
    never go stale) removes the redundant work entirely without changing
    any computed value."""
    cache = getattr(L_csc, "_sedumipy_solve_cache", None)
    if cache is None:
        Ljc, _ = _as_index_array(L_csc.indptr)
        Lir, _ = _as_index_array(L_csc.indices)
        Lpr, _ = _as_double_array(L_csc.data)
        cache = (Ljc, Ljc.ctypes.data_as(c_size_t_p), Lir, Lir.ctypes.data_as(c_size_t_p), Lpr, Lpr.ctypes.data_as(c_double_p))
        L_csc._sedumipy_solve_cache = cache
    return cache


def fwsolve(L_csc, xsuper, y):
    """Forward-solve y := L \\ y in place, where L is the unit-lower-
    triangular factor from SeDuMi's supernodal Cholesky (blkchol), stored
    exactly like `scipy.sparse.csc_matrix` (indptr/indices/data), and
    `xsuper` marks each supernode's first column (length nsuper+1, as
    produced by symbchol). L's stored diagonal entries are never read --
    only their presence in the sparsity pattern matters, per fwblkslv.c.

    This wraps fwsolve() in fwblkslv.c (sedumi's forward substitution
    kernel) with no MATLAB/Octave/MEX layer at all.
    """
    Ljc, Ljc_p, Lir, Lir_p, Lpr, Lpr_p = _cached_csc_solve_arrays(L_csc)
    xs, xs_p = _as_index_array(xsuper)
    ya, y_p = _as_double_array(y)

    nsuper = xs.size - 1
    m = L_csc.shape[0]
    fwork = ctypes.create_string_buffer(8 * m)  # generous upper bound
    fwork_p = ctypes.cast(fwork, c_double_p)

    _lib.fwsolve(y_p, Ljc_p, Lir_p, Lpr_p, xs_p, nsuper, fwork_p)
    return ya


def bwsolve(L_csc, xsuper, y):
    """Backward-solve y := L' \\ y in place -- see fwsolve()'s docstring;
    wraps bwsolve() in bwblkslv.c."""
    Ljc, Ljc_p, Lir, Lir_p, Lpr, Lpr_p = _cached_csc_solve_arrays(L_csc)
    xs, xs_p = _as_index_array(xsuper)
    ya, y_p = _as_double_array(y)

    nsuper = xs.size - 1
    m = L_csc.shape[0]
    fwork = ctypes.create_string_buffer(8 * m)  # generous upper bound
    fwork_p = ctypes.cast(fwork, c_double_p)

    _lib.bwsolve(y_p, Ljc_p, Lir_p, Lpr_p, xs_p, nsuper, fwork_p)
    return ya


def _to_fortran_adjacency(A_csc):
    """Mirrors getadj() in ordmmdmex.c/symfctmex.c exactly: converts a
    0-indexed CSC sparsity pattern (diagonal entries included or not, both
    fine) into the 1-indexed (xadj, adjncy) adjacency-list form Liu's
    Fortran ordmmd_/symfct_ expect, dropping diagonal entries. `A_csc`
    must be symmetric (only one triangle need be nonzero for the sparsity
    pattern, but symbchol.m always passes a symmetric ADA_sedumi_) and is
    sorted by row index within each column first, to match MATLAB/
    Octave's own sparse storage -- this determines the exact adjacency
    order ordmmd_/symfct_ see, which affects tie-breaking, so must match
    bit-for-bit to reproduce the same ordering as the MEX build.
    """
    import numpy as np

    A = A_csc.copy()
    A.sort_indices()
    n = A.shape[0]
    cjc = A.indptr
    cir = A.indices

    xadj = np.empty(n + 1, dtype=np.intp)
    adjncy = np.empty(cjc[n], dtype=np.intp)  # upper bound; diag entries excluded
    inz = 0
    for j in range(n):
        xadj[j] = inz + 1
        for ix in range(cjc[j], cjc[j + 1]):
            i = cir[ix]
            if i != j:
                adjncy[inz] = i + 1
                inz += 1
    xadj[n] = inz + 1
    return xadj, adjncy[:inz].copy()


def ordmmd(A_csc):
    """Multiple minimum degree ordering (Liu's genmmd, via ordmmd.c) --
    the exact same algorithm and C source SeDuMi's ordmmdmex() MEX
    function calls, just without any MATLAB/Octave/MEX in the calling
    path. `A_csc` must be a square, symmetric scipy.sparse.csc_matrix (as
    symbchol.m always passes -- only the sparsity pattern matters).

    Returns a 0-indexed permutation array `perm` (Python/NumPy
    convention); MATLAB/Octave's ordmmdmex returns the 1-indexed version
    of the exact same array.
    """
    import numpy as np

    n = A_csc.shape[0]
    if A_csc.shape[0] != A_csc.shape[1]:
        raise ValueError("A_csc must be square")

    xadj, adjncy = _to_fortran_adjacency(A_csc)
    perm = np.zeros(n, dtype=np.intp)
    invp = np.zeros(n, dtype=np.intp)
    iwsiz = 4 * n
    iwork = np.zeros(max(iwsiz, 1), dtype=np.intp)

    neqns = ctypes.c_ssize_t(n)
    iwsiz_c = ctypes.c_ssize_t(iwsiz)
    nofsub = ctypes.c_ssize_t(0)
    iflag = ctypes.c_ssize_t(0)

    def p(arr):
        return arr.ctypes.data_as(c_ssize_t_p)

    _lib.ordmmd_(
        ctypes.byref(neqns),
        p(xadj),
        p(adjncy),
        p(invp),
        p(perm),
        ctypes.byref(iwsiz_c),
        p(iwork),
        ctypes.byref(nofsub),
        ctypes.byref(iflag),
    )
    if iflag.value == -1:
        raise RuntimeError("ordmmd: insufficient working storage (iwsiz too small)")
    return perm - 1


def symbolic_cholesky(A_csc, perm0):
    """Symbolic block-sparse Cholesky factorization: mirrors
    `L = symfctmex(X, perm)` (symfctmex.c, Liu's SPARSPAK-A sfinit_/
    symfct_) exactly -- same C/Fortran source, no MATLAB/Octave/MEX in
    the calling path. This is the real (non-dense-fallback) branch of
    symbchol.m.

    Parameters
    ----------
    A_csc : scipy.sparse.csc_matrix, symmetric m x m sparsity pattern
        (as symbchol.m's ADA_sedumi_).
    perm0 : 0-indexed initial permutation (e.g. from ordmmd()).

    Returns
    -------
    dict with:
      "L"      -- scipy.sparse.csc_matrix, the symbolic factor's
                  sparsity pattern (all data entries are 1.0 -- this is
                  a symbolic factorization, not a numeric one), matching
                  the MATLAB/Octave L.L field.
      "perm"   -- 0-indexed final permutation (matching L.perm; sfinit_
                  can refine the input ordmmd() permutation to merge
                  same-pattern columns into supernodes, so this is not
                  always identical to perm0).
      "xsuper" -- 0-indexed supernode boundaries, length nsuper+1
                  (matching L.xsuper).
    """
    import numpy as np

    m = A_csc.shape[0]
    if A_csc.shape[0] != A_csc.shape[1]:
        raise ValueError("A_csc must be square")

    xadj, adjncy = _to_fortran_adjacency(A_csc)
    nnza = A_csc.nnz  # matches Xjc[m] in symfctmex.c: X's own total stored
    # nnz (diagonal included), used only as an upper bound for workspace
    # sizing, not literally len(adjncy) (which excludes the diagonal).

    perm = (np.asarray(perm0, dtype=np.intp) + 1).copy()  # 1-indexed, mutable:
    # sfinit_ updates (perm, invp) in place to an equivalent ordering, so
    # this buffer's *final* contents (after both Fortran calls below) is
    # the true output permutation -- exactly what symfctmex.c reads back.
    invp = np.zeros(m, dtype=np.intp)
    for i in range(m):
        invp[perm[i] - 1] = i + 1

    colcnt = np.zeros(m, dtype=np.intp)
    snode = np.zeros(m, dtype=np.intp)
    xsuper = np.zeros(m + 1, dtype=np.intp)
    iwsiz = 7 * m + 3
    iwork = np.zeros(max(iwsiz, 1), dtype=np.intp)

    def p(arr):
        return arr.ctypes.data_as(c_ssize_t_p)

    m_c = ctypes.c_ssize_t(m)
    nnza_c = ctypes.c_ssize_t(nnza)
    nnzl = ctypes.c_ssize_t(0)
    nsub = ctypes.c_ssize_t(0)
    nsuper_c = ctypes.c_ssize_t(0)
    iwsiz_c = ctypes.c_ssize_t(iwsiz)
    flag = ctypes.c_ssize_t(0)

    _lib.sfinit_(
        ctypes.byref(m_c), ctypes.byref(nnza_c), p(xadj), p(adjncy),
        p(perm), p(invp), p(colcnt),
        ctypes.byref(nnzl), ctypes.byref(nsub), ctypes.byref(nsuper_c),
        p(snode), p(xsuper), ctypes.byref(iwsiz_c), p(iwork), ctypes.byref(flag),
    )
    if flag.value == -1:
        raise RuntimeError("sfinit: insufficient working storage (iwsiz too small)")

    xlindx = np.zeros(m + 1, dtype=np.intp)
    # Lir/Ljc are allocated once and reused for two different purposes in
    # sequence -- exactly as symfctmex.c does: symfct_ first fills them as
    # (lindx, xlnz) -- a compact, per-supernode, 1-indexed representation
    # -- then expandsub() overwrites the same buffers in place with the
    # standard, per-column, 0-indexed CSC representation (Ljc capacity
    # m+1 is already enough for both; Lir's nnzl capacity from sfinit_ is
    # an upper bound valid for both the compact and the expanded form).
    Lir = np.zeros(max(nnzl.value, 1), dtype=np.intp)
    Ljc = np.zeros(m + 1, dtype=np.intp)

    flag2 = ctypes.c_ssize_t(0)
    _lib.symfct_(
        ctypes.byref(m_c), ctypes.byref(nnza_c), p(xadj), p(adjncy),
        p(perm), p(invp), p(colcnt),
        ctypes.byref(nsuper_c), p(xsuper), p(snode),
        ctypes.byref(nsub), p(xlindx), p(Lir), p(Ljc),
        ctypes.byref(iwsiz_c), p(iwork), ctypes.byref(flag2),
    )
    if flag2.value == -1:
        raise RuntimeError("symfct: insufficient working storage (iwsiz too small)")
    if flag2.value == -2:
        raise RuntimeError("symfct: input error")

    def p_size_t(arr):
        return arr.ctypes.data_as(c_size_t_p)

    _lib.expandsub(m, nsuper_c.value, p_size_t(xsuper), p_size_t(xlindx),
                    p_size_t(Ljc), p_size_t(Lir))

    nnz_L = int(Ljc[m])
    import scipy.sparse

    L_csc = scipy.sparse.csc_matrix(
        (np.ones(nnz_L, dtype=np.float64), Lir[:nnz_L].astype(np.int64),
         Ljc.astype(np.int64)),
        shape=(m, m),
    )

    final_xsuper = xsuper[: nsuper_c.value + 1] - 1
    snode_for_tmpsiz = _compute_snode(final_xsuper, m).copy()
    tmpsiz = _lib.gettmpsiz(
        Ljc.astype(np.uintp).ctypes.data_as(c_size_t_p),
        Lir[:nnz_L].astype(np.uintp).ctypes.data_as(c_size_t_p),
        final_xsuper.astype(np.uintp).ctypes.data_as(c_size_t_p),
        nsuper_c.value,
        snode_for_tmpsiz.ctypes.data_as(c_size_t_p),
    )

    return {
        "L": L_csc,
        "perm": perm - 1,
        "xsuper": final_xsuper,
        "tmpsiz": int(tmpsiz),
    }


def symbolic_cholesky_dense(m: int) -> dict:
    """Mirrors symbchol.m's `else` branch (`spars(ADA)==1`, a fully
    dense matrix): skips minimum-degree ordering entirely and returns
    the trivial one-big-supernode symbolic structure directly --
    perm=identity, L = tril(ones(m,m)), xsuper=[0,m] (0-indexed).
    tmpsiz is still computed via the real gettmpsiz() kernel (same one
    symbolic_cholesky() uses) so numeric_cholesky()'s workspace sizing
    matches blkchol.c exactly for this branch too.
    """
    import numpy as np
    import scipy.sparse

    perm = np.arange(m, dtype=np.int64)
    xsuper = np.array([0, m], dtype=np.int64)
    L_csc = scipy.sparse.csc_matrix(np.tril(np.ones((m, m), dtype=np.float64)))
    Ljc = np.ascontiguousarray(L_csc.indptr, dtype=np.uintp)
    Lir = np.ascontiguousarray(L_csc.indices, dtype=np.uintp)
    snode = _compute_snode(xsuper, m)
    tmpsiz = _lib.gettmpsiz(
        Ljc.ctypes.data_as(c_size_t_p),
        Lir.ctypes.data_as(c_size_t_p),
        xsuper.astype(np.uintp).ctypes.data_as(c_size_t_p),
        1,
        snode.ctypes.data_as(c_size_t_p),
    )
    return {"L": L_csc, "perm": perm, "xsuper": xsuper, "tmpsiz": int(tmpsiz)}


def numeric_cholesky(sym: dict, X_csc, pars: dict | None = None, absd=None) -> dict:
    """Numeric block sparse LDL' factorization: mirrors
    `[L.L, L.d, skip, diagadd] = blkchol(L, X, pars, absd)` (blkchol.c's
    permuteP + spchol/blkLDL) exactly, no MATLAB/Octave/MEX in the
    calling path.

    Parameters
    ----------
    sym : dict from symbolic_cholesky() -- needs "L" (pattern), "perm"
        (0-indexed), "xsuper" (0-indexed), "tmpsiz".
    X_csc : scipy.sparse.csc_matrix, the numeric matrix to factor (only
        its lower triangle is read, matching P(perm,perm) -- SeDuMi
        always builds X to already be symmetric with only tril stored).
    pars : optional dict overriding canceltol (1e-12), maxu (5e2),
        abstol (1e-20), delay (False) -- same names/defaults as blkchol.c.
    absd : optional length-m array of "before cancellation" diagonal
        magnitudes (pars.absd in the .m API); defaults to X's own
        diagonal.

    Returns dict with "L" (scipy.sparse.csc_matrix, same pattern as
    sym["L"], numeric values), "d" (length-m diagonal of D, with
    d[skip]==0 as blkchol always reports it), "skip" (0-indexed columns
    where the pivot was too unstable and got replaced with a unit
    column -- matching `L.d(find(L.skip)) = inf` in blkchol.m's own
    solve recipe), "diagadd" (values added to the diagonal at the
    OTHER unstable pivots that were stabilized instead of skipped).
    """
    import numpy as np
    import scipy.sparse

    m = X_csc.shape[0]
    L_pattern = sym["L"].tocsc()
    perm = np.ascontiguousarray(sym["perm"], dtype=np.uintp)
    xsuper = np.ascontiguousarray(sym["xsuper"], dtype=np.uintp)
    nsuper = xsuper.size - 1
    tmpsiz = sym["tmpsiz"]

    pars = pars or {}
    canceltol = float(pars.get("canceltol", 1e-12))
    maxu = float(pars.get("maxu", 5e2))
    abstol = max(float(pars.get("abstol", 1e-20)), 0.0)
    use_delay = bool(pars.get("delay", False))

    def p_size_t(arr):
        return arr.ctypes.data_as(c_size_t_p)

    def p_double(arr):
        return arr.ctypes.data_as(c_double_p)

    Ljc = np.ascontiguousarray(L_pattern.indptr, dtype=np.uintp)
    Lir_original = np.ascontiguousarray(L_pattern.indices, dtype=np.uintp)
    Lir = Lir_original.copy()  # spchol uses this as scratch, restored after
    Lpr = np.zeros(int(Ljc[-1]), dtype=np.float64)

    Pj = np.zeros(max(m, 1), dtype=np.float64)  # permuteP's own scratch space
    X = X_csc.tocsc()
    Pjc = np.ascontiguousarray(X.indptr, dtype=np.uintp)
    Pir = np.ascontiguousarray(X.indices, dtype=np.uintp)
    Ppr = np.ascontiguousarray(X.data, dtype=np.float64)

    _lib.permuteP(
        p_size_t(Ljc), p_size_t(Lir), p_double(Lpr),
        p_size_t(Pjc), p_size_t(Pir), p_double(Ppr),
        p_size_t(perm), p_double(Pj), m,
    )

    if absd is not None:
        absd_arr = np.ascontiguousarray(absd, dtype=np.float64)
        orgd = absd_arr[perm.astype(np.int64)].copy()
    else:
        orgd = np.array([Lpr[Ljc[j]] for j in range(m)], dtype=np.float64)

    snode = np.zeros(max(m, 1), dtype=np.uintp)
    xlindx = np.zeros(m + 1, dtype=np.uintp)
    d = np.zeros(m, dtype=np.float64)
    skip = np.zeros(max(m, 1), dtype=np.uintp)
    nadd_c = ctypes.c_size_t(0)

    iwsiz = max(2 * (m + nsuper), 1)
    fwsiz = max(tmpsiz, 1)
    iwork = np.zeros(iwsiz, dtype=np.uintp)
    fwork = np.zeros(fwsiz, dtype=np.float64)

    nskip = _lib.spchol(
        m, nsuper, p_size_t(xsuper),
        p_size_t(snode), p_size_t(xlindx), p_size_t(Lir), p_double(orgd),
        p_size_t(Ljc), p_double(Lpr), p_double(d), p_size_t(perm),
        abstol, canceltol, maxu,
        p_size_t(skip), ctypes.byref(nadd_c),
        iwsiz, p_size_t(iwork), fwsiz, p_double(fwork),
    )
    if nskip == _SIZE_T_ERROR:
        raise RuntimeError("spchol: insufficient working storage (iwsiz/fwsiz too small)")
    nadd = nadd_c.value

    # spchol used Lir as scratch (the "compress subscripts" step turns it
    # into a per-supernode compact array); the *sparsity pattern* of the
    # output L is unchanged from the input, so restore it -- exactly the
    # memcpy(L.ir, LINir, ...) in blkchol.c's mexFunction.
    Lir[:] = Lir_original

    skip = skip[:nskip]
    diagadd_idx = np.zeros(nadd, dtype=np.uintp)
    diagadd_val = np.zeros(nadd, dtype=np.float64)

    skip_out = []
    skip_val = []
    for j in range(nskip):
        i = int(skip[j])
        if use_delay:
            skip_val.append(1.0)
        else:
            skip_val.append(Lpr[Ljc[i]])
            Lpr[Ljc[i]] = 1.0
            Lpr[Ljc[i] + 1 : Ljc[i + 1]] = 0.0
        skip_out.append(i)

    # iwork[:nadd] holds the diagadd indices, written by spchol.
    for j in range(nadd):
        i = int(iwork[j])
        diagadd_idx[j] = i
        diagadd_val[j] = orgd[i]

    L_csc = scipy.sparse.csc_matrix(
        (Lpr, Lir.astype(np.int64), Ljc.astype(np.int64)), shape=(m, m)
    )
    # Ljc/Lir are already uintp and Lpr already float64 right here -- feed
    # them straight into fwsolve()/bwsolve()'s cache (see
    # _cached_csc_solve_arrays()) instead of letting it re-derive uintp
    # copies later from L_csc.indptr/.indices, which are the int64 copies
    # scipy just normalized them into two lines up. Skips one of the two
    # redundant dtype round-trips a large nnz(L) would otherwise pay for.
    L_csc._sedumipy_solve_cache = (
        Ljc, Ljc.ctypes.data_as(c_size_t_p),
        Lir, Lir.ctypes.data_as(c_size_t_p),
        Lpr, Lpr.ctypes.data_as(c_double_p),
    )

    return {
        "L": L_csc,
        "d": d,
        "skip": np.array(skip_out, dtype=np.int64),
        "skip_values": np.array(skip_val, dtype=np.float64),
        "diagadd_index": diagadd_idx.astype(np.int64),
        "diagadd": diagadd_val,
    }


_lib.ddotxj.argtypes = [c_double_p, c_double_p, c_double_p, c_size_t_p, ctypes.c_size_t]
_lib.ddotxj.restype = None

_lib.blkmul.argtypes = [c_double_p, c_double_p, c_double_p, c_size_t_p, ctypes.c_size_t, ctypes.c_size_t]
_lib.blkmul.restype = ctypes.c_int

_lib.vecsymPSD.argtypes = [c_double_p, c_double_p, ctypes.c_size_t, ctypes.c_size_t, c_double_p]
_lib.vecsymPSD.restype = None

_lib.rquaddadd.argtypes = [c_double_p, ctypes.c_double, ctypes.c_double, ctypes.c_double]
_lib.rquaddadd.restype = ctypes.c_double


def ddot(d, X, blkstart):
    """ddot(d, X, blkstart) -- dense-X path of ddot.c/ddotxj: for each
    column of X and each Lorentz block k (spanning blkstart[k]:blkstart[k+1]
    in 0-indexed, half-open convention), computes d[k]'*X[block,column].
    Wraps ddotxj() directly (no MATLAB/Octave/MEX); the sparse-X path
    (spddotxj) is not yet wrapped -- Phase 3 will add it if/when the .m
    port actually needs ddot on sparse X.

    `blkstart` is the real, absolute (possibly >1) 1-indexed .m/MEX-style
    array; `X` (like `d`) may be given either as the exact qDim-row block
    span, or as a full, absolute-indexed array -- ddot.c's own
    mexFunction (not just its ddotxj() core) applies a row offset in the
    latter case, so this wrapper replicates that same case analysis
    (nrows==qDim / nrows==nblk+qDim / general nrows>=blkstart[-1]) rather
    than just handling the qDim-exact case ddotxj() itself assumes.
    """
    import numpy as np

    d = np.ascontiguousarray(d, dtype=np.float64)
    X = np.ascontiguousarray(X, dtype=np.float64)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    blkstart = np.ascontiguousarray(blkstart, dtype=np.int64)
    nblk = blkstart.size - 1
    nrows, ncols = X.shape

    qDim = int(blkstart[-1] - blkstart[0])
    if d.size != qDim:
        d = d[int(blkstart[0]) :]

    # 0-indexed absolute positions -- mirrors the mexFunction's own
    # `blkstart[i] = (mwIndex)blkstartPr[i] - 1` conversion; the row-offset
    # comparisons below are against these, not the raw 1-indexed values.
    blkstart0 = blkstart - 1
    row_off = 0
    if nrows != qDim:
        if nrows < int(blkstart0[-1]):
            if nrows != nblk + qDim:
                raise ValueError("ddot: X size mismatch")
            row_off = nblk  # Lorentz trace + norm-bound layout
        else:
            row_off = int(blkstart0[0])  # X is the full, absolute-indexed array

    out = np.empty((nblk, ncols), dtype=np.float64, order="F")
    bs = (blkstart - blkstart[0]).astype(np.uintp)  # ddotxj asserts blkstart[0]==0
    dptr = d.ctypes.data_as(c_double_p)
    bsptr = bs.ctypes.data_as(c_size_t_p)
    for j in range(ncols):
        col = np.ascontiguousarray(X[row_off:, j])
        _lib.ddotxj(
            out[:, j].ctypes.data_as(c_double_p), dptr,
            col.ctypes.data_as(c_double_p), bsptr, nblk,
        )
    return out.squeeze(axis=1) if ncols == 1 else out


def blkmul(mu, d, nL):
    """y[block k] = mu[k] * d[block k], blocks given by nL (block LENGTHS,
    not offsets -- see blkmul.m). Wraps blkmul() (blkmul.c) directly."""
    import numpy as np

    mu = np.ascontiguousarray(mu, dtype=np.float64).ravel()
    d = np.ascontiguousarray(d, dtype=np.float64).ravel()
    nL_arr = np.ascontiguousarray(nL, dtype=np.uintp).ravel()
    kappa = mu.size
    n = d.size
    y = np.zeros(n, dtype=np.float64)

    remaining = _lib.blkmul(
        y.ctypes.data_as(c_double_p), mu.ctypes.data_as(c_double_p),
        d.ctypes.data_as(c_double_p), nL_arr.ctypes.data_as(c_size_t_p),
        kappa, n,
    )
    if remaining != 0:
        raise ValueError("blkmul: nL size mismatch (sum(nL) != len(d))")
    return y


def qblkmul(mu, d, blkstart):
    """y[block k] = mu[k] * d[block k], blocks given by blkstart
    (1-indexed, as in the .m/MEX convention -- see qblkmul.m). qblkmul.c
    has no separable core function (the whole computation lives in its
    mexFunction), so this ports that logic directly to NumPy rather than
    binding a C function that doesn't exist as such.

    Vectorized via np.repeat rather than a per-block Python for loop: on
    a problem with many Lorentz blocks (e.g. DIMACS qssp180old, 65341 of
    them), a `for k in range(nblk)` loop pays Python-interpreter
    per-iteration overhead tens of thousands of times per call (cProfile
    on 5 outer iterations showed this function's own loop body at 49.7s
    of self time across 334 calls). np.repeat(mu, block_sizes) expands
    each mu[k] to its block's width in one vectorized call, computing the
    exact same per-element products mu[k]*d[...] in the same order --
    bug-for-bug identical output, no Python-level loop."""
    import numpy as np

    mu = np.ascontiguousarray(mu, dtype=np.float64).ravel()
    d = np.ascontiguousarray(d, dtype=np.float64).ravel()
    blkstart = np.ascontiguousarray(blkstart, dtype=np.int64).ravel() - 1
    nblk = mu.size
    span = int(blkstart[-1] - blkstart[0])

    if d.size != span:
        if d.size == nblk + span:
            d = d[nblk:]
        else:
            d = d[int(blkstart[0]) :]
    # The loop this replaces only ever reads positions [0, span) of the
    # (already offset-adjusted) d above -- match that truncation exactly,
    # since d can still be longer than span here (e.g. when the caller
    # passes the whole state vector and blkstart[0] > 0 trims only the
    # leading LP/arrow part, not any unrelated trailing data).
    d = d[:span]

    block_sizes = np.diff(blkstart)
    return np.repeat(mu, block_sizes) * d


def vecsym(x, K: dict):
    """y = vecsym(x, K): copies the LP+SOCP part of x unchanged, then
    symmetrizes each real PSD block ((Xk+Xk')/2) and Hermitianizes each
    complex one, via vecsymPSD() (vecsym.c) directly."""
    import numpy as np

    cK = cone_from_dict(K)
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    lqDim = cK.lpN + cK.qDim
    lenfull = lqDim + cK.rDim + cK.hDim
    if x.size != lenfull:
        raise ValueError(f"x must have length {lenfull}, got {x.size}")

    y = x.copy()
    sdpNL = cK._keepalive[2]  # the "s" array cone_from_dict built cK from
    _lib.vecsymPSD(
        y[lqDim:].ctypes.data_as(c_double_p),
        x[lqDim:].ctypes.data_as(c_double_p),
        cK.rsdpN, cK.sdpN,
        sdpNL.ctypes.data_as(c_double_p) if cK.sdpN else None,
    )
    return y


def quadadd(xhi, xlo, y):
    """(zhi, zlo) = quadadd(xhi, xlo, y): extended-precision (double-
    double style) addition xhi+xlo+y, elementwise, via rquaddadd()
    (quadadd.c) directly -- not reimplemented in Python, since the whole
    point of this kernel is the specific extended-precision arithmetic
    sequence, which is easy to get subtly wrong by "simplifying"."""
    import numpy as np

    xhi = np.ascontiguousarray(xhi, dtype=np.float64).ravel()
    xlo = np.ascontiguousarray(xlo, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y, dtype=np.float64).ravel()
    m = xhi.size
    zhi = np.empty(m, dtype=np.float64)
    zlo = np.empty(m, dtype=np.float64)
    for i in range(m):
        lo = ctypes.c_double(0.0)
        zhi[i] = _lib.rquaddadd(ctypes.byref(lo), xhi[i], xlo[i], y[i])
        zlo[i] = lo.value
    return zhi, zlo


class KeyDouble(ctypes.Structure):
    """Mirrors blksdp.h's `keydouble` (double r; mwIndex k;)."""

    _fields_ = [("r", ctypes.c_double), ("k", ctypes.c_size_t)]


c_ubyte_p = ctypes.POINTER(ctypes.c_ubyte)  # for `char*`/`bool*` buffers --
# deliberately not ctypes.c_char_p, which has Python-string marshaling
# semantics that are the wrong fit for a plain output byte buffer.

_lib.fwprodform.argtypes = [
    c_double_p, c_size_t_p, c_size_t_p, c_double_p, c_double_p, c_size_t_p,
    c_ubyte_p, ctypes.c_size_t,
]
_lib.fwprodform.restype = None

_lib.bwprodform.argtypes = [
    c_double_p, c_size_t_p, c_size_t_p, c_double_p, c_double_p, c_size_t_p,
    c_ubyte_p, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
]
_lib.bwprodform.restype = None

_lib.prodformfact.argtypes = [
    c_double_p, c_size_t_p, c_double_p, c_size_t_p,
    c_double_p, c_ubyte_p, c_size_t_p,
    c_size_t_p, c_size_t_p,
    c_double_p, ctypes.c_size_t, c_size_t_p, c_size_t_p,
    ctypes.c_double, c_double_p, ctypes.POINTER(KeyDouble),
]
_lib.prodformfact.restype = None


def _dpr1_apply(direction, Lden: dict, b):
    """Shared implementation for fwdpr1()/bwdpr1(): PROD_k L(pk,betak) *
    ynew = yold (forward) or PROD_k L(pk,betak)' * ynew = yold (backward),
    where L(p,beta) = eye(m) + tril(p*beta',-1). Wraps fwprodform/
    bwprodform (fwdpr1.c/bwdpr1.c) directly.

    Lden needs "betajc" (1-indexed, length nden+1; nden==0 means no dense
    columns at all -- y is returned unchanged, exactly like the MEX
    build's early return), "p", "beta", "pivperm", "dopiv", and "dz" (a
    sparse matrix whose .indptr is used as the per-column cumulative
    "reach" xsuper and whose .indices give the row-compaction mapping --
    see dpr1fact()'s docstring for what this "dz" scheme actually means).
    """
    import numpy as np

    betajc_1indexed = np.ascontiguousarray(Lden["betajc"], dtype=np.int64).ravel()
    nden = betajc_1indexed.size - 1
    b = np.ascontiguousarray(b, dtype=np.float64)
    if nden == 0:
        return b.copy()

    single_col = b.ndim == 1
    B = b.reshape(-1, 1) if single_col else b
    m = B.shape[0]

    betajc = (betajc_1indexed - 1).astype(np.uintp)
    p = np.ascontiguousarray(Lden["p"], dtype=np.float64)
    beta = np.ascontiguousarray(Lden["beta"], dtype=np.float64)
    ordered = np.ascontiguousarray(Lden["dopiv"], dtype=np.uint8).ravel()
    # pivperm is an opaque, internal 0-indexed array private to the
    # dpr1fact()<->fwdpr1()/bwdpr1() pairing -- never interpreted as a
    # MATLAB-facing 1-indexed permutation anywhere, including by the
    # original C mexFunctions (they round-trip it unchanged), so no
    # +-1 conversion happens here either.
    pivperm = np.ascontiguousarray(Lden["pivperm"], dtype=np.uintp).ravel()
    dz_jc = np.ascontiguousarray(Lden["dz"].indptr, dtype=np.uintp)
    dz_ir = np.ascontiguousarray(Lden["dz"].indices, dtype=np.uintp)
    dznnz = int(dz_jc[nden])

    Y = B.copy()
    fwork = np.empty(max(dznnz, 1), dtype=np.float64)
    for j in range(Y.shape[1]):
        col = Y[:, j]
        for i in range(dznnz):
            fwork[i] = col[dz_ir[i]]
        if direction == "fw":
            _lib.fwprodform(
                fwork.ctypes.data_as(c_double_p), dz_jc.ctypes.data_as(c_size_t_p),
                pivperm.ctypes.data_as(c_size_t_p), p.ctypes.data_as(c_double_p),
                beta.ctypes.data_as(c_double_p), betajc.ctypes.data_as(c_size_t_p),
                ordered.ctypes.data_as(c_ubyte_p), nden,
            )
        else:
            # bwprodform additionally needs the *total* lengths of p and
            # pivperm up front (it walks backward, decrementing into
            # them), unlike fwprodform which only needs cumulative
            # offsets it can derive from dz_jc as it goes forward.
            _lib.bwprodform(
                fwork.ctypes.data_as(c_double_p), dz_jc.ctypes.data_as(c_size_t_p),
                pivperm.ctypes.data_as(c_size_t_p), p.ctypes.data_as(c_double_p),
                beta.ctypes.data_as(c_double_p), betajc.ctypes.data_as(c_size_t_p),
                ordered.ctypes.data_as(c_ubyte_p), nden, p.size, pivperm.size,
            )
        for i in range(dznnz):
            col[dz_ir[i]] = fwork[i]
    return Y[:, 0] if single_col else Y


def fwdpr1(Lden: dict, b):
    """y = fwdpr1(Lden, b): solve PROD_k L(pk,betak) * y = b. See
    _dpr1_apply()'s docstring for Lden's fields."""
    return _dpr1_apply("fw", Lden, b)


def bwdpr1(Lden: dict, b):
    """y = bwdpr1(Lden, b): solve PROD_k L(pk,betak)' * y = b. See
    _dpr1_apply()'s docstring for Lden's fields."""
    return _dpr1_apply("bw", Lden, b)


def dpr1fact(x, d, Lsym: dict, smult, maxu: float):
    """[Lden, d_out] = dpr1fact(x, d, Lsym, smult, maxu): factors
    diag(d) + x*diag(smult)*x' = (PROD_k L(pk,betak)) * diag(d_out) *
    (PROD_k L(pk,betak))', wrapping prodformfact() (dpr1fact.c) directly
    -- the same C computation the MEX build uses for SeDuMi's
    "dense column" handling in the PCG preconditioner.

    x : scipy.sparse.csc_matrix, m x n (n = number of dense columns).
    d : length-m array, diagonal to update.
    Lsym : dict with "dz" (scipy.sparse.csc_matrix -- see below), "perm"
        (1-indexed length-n column order, as symfctmex-style outputs
        use), "first" (1-indexed length-n array).
    smult : length-n array of per-column multipliers (x*diag(smult)*x').
    maxu : stability threshold -- a column gets pivoted (reordered) if a
        pivot magnitude ratio would otherwise exceed this.

    Returns a dict shaped like Lden (see _dpr1_apply()'s docstring) plus
    the updated diagonal, ready to feed to fwdpr1()/bwdpr1() after also
    copying over Lsym's dz/perm/first fields (as deninfac.m does):
        Lden["dz"], Lden["first"], Lden["perm"] = Lsym["dz"], Lsym["first"], Lsym["perm"]
    """
    import numpy as np
    import scipy.sparse

    X = x.tocsc()
    m, n = X.shape
    dz = Lsym["dz"].tocsc()
    dz_jc = np.ascontiguousarray(dz.indptr, dtype=np.uintp)
    dz_ir = np.ascontiguousarray(dz.indices, dtype=np.uintp)
    dznnz = int(dz_jc[n])
    if dznnz > m:
        raise ValueError("Lsym.dz size mismatch: more compact rows than m")

    colperm = (np.ascontiguousarray(Lsym["perm"], dtype=np.int64).ravel() - 1).astype(np.uintp)
    firstpiv = (np.ascontiguousarray(Lsym["first"], dtype=np.int64).ravel() - 1).astype(np.uintp)

    pnnz = int(sum(int(dz_jc[j + 1]) for j in range(n)))
    d_compact = np.empty(max(dznnz, 1), dtype=np.float64)
    lab = np.ascontiguousarray(d, dtype=np.float64).copy()
    for i in range(dznnz):
        d_compact[i] = lab[dz_ir[i]]

    dep = np.zeros(dznnz + 1, dtype=np.uintp)
    ndep = 0
    for i in range(dznnz):
        if d_compact[i] <= 0.0:
            dep[ndep] = i
            ndep += 1
    dep[ndep] = m

    invrowperm = np.zeros(max(m, 1), dtype=np.uintp)
    for i in range(dznnz):
        invrowperm[dz_ir[i]] = i

    p = np.zeros(max(pnnz + m, 1), dtype=np.float64)
    pos = 0
    for j in range(n):
        pos += int(dz_jc[j])
        permj = int(colperm[j])
        for i in range(X.indptr[permj], X.indptr[permj + 1]):
            p[pos + int(invrowperm[X.indices[i]])] = X.data[i]

    p = p[:pnnz].copy() if pnnz > 0 else np.zeros(0, dtype=np.float64)
    beta = np.zeros(max(pnnz, 1), dtype=np.float64)
    betajc = np.zeros(n + 1, dtype=np.uintp)
    ordered = np.zeros(max(n, 1), dtype=np.uint8)
    pivperm = np.zeros(max(pnnz, 1), dtype=np.uintp)
    fwork = np.zeros(max(dznnz, 1), dtype=np.float64)
    kdwork = (KeyDouble * max(dznnz, 1))()
    ndep_c = ctypes.c_size_t(ndep)
    smult_arr = np.ascontiguousarray(smult, dtype=np.float64)

    _lib.prodformfact(
        p.ctypes.data_as(c_double_p), pivperm.ctypes.data_as(c_size_t_p),
        beta.ctypes.data_as(c_double_p), betajc.ctypes.data_as(c_size_t_p),
        d_compact.ctypes.data_as(c_double_p), ordered.ctypes.data_as(c_ubyte_p),
        dz_jc.ctypes.data_as(c_size_t_p),
        colperm.ctypes.data_as(c_size_t_p), firstpiv.ctypes.data_as(c_size_t_p),
        smult_arr.ctypes.data_as(c_double_p), n, dep.ctypes.data_as(c_size_t_p),
        ctypes.byref(ndep_c),
        maxu, fwork.ctypes.data_as(c_double_p), kdwork,
    )

    for i in range(dznnz):
        lab[dz_ir[i]] = d_compact[i]

    # permnnz = sum{dz.jc[j+1] | ordered[j]==1} -- exactly dpr1fact.c's
    # mexFunction; pivperm[:permnnz] is meaningful, the rest is scratch.
    permnnz = 0
    for i in range(n):
        if ordered[i]:
            permnnz += int(dz_jc[i + 1])

    Lden = {
        "betajc": (betajc[: n + 1].astype(np.int64) + 1),  # 1-indexed, .m-facing
        "beta": beta[: int(betajc[n])].copy(),
        "p": p,
        "dopiv": ordered[:n].copy(),
        "pivperm": pivperm[:permnnz].copy(),  # opaque, see _dpr1_apply()
    }
    return Lden, lab


for _name, _argtypes, _restype in [
    ("matgivens", [c_double_p, c_double_p, c_size_t_p, ctypes.c_size_t], None),
    ("rotorder", [c_size_t_p, c_double_p, c_size_t_p, c_double_p, c_double_p,
                  ctypes.c_double, ctypes.c_size_t], None),
    ("qdivv", [c_double_p, c_double_p, c_double_p, ctypes.c_size_t], None),
    ("psdframeit", [c_double_p, c_double_p, c_double_p, c_size_t_p,
                     ctypes.c_size_t, ctypes.c_size_t, c_double_p], None),
    ("psdinvjmul", [c_double_p, c_double_p, c_double_p, c_double_p, c_size_t_p,
                     ctypes.c_size_t, ctypes.c_size_t, c_double_p], None),
    ("qrfac", [c_double_p, c_double_p, c_double_p, ctypes.c_size_t], None),
    ("utmulx", [c_double_p, c_double_p, c_double_p, ctypes.c_size_t], None),
    ("triu2sym", [c_double_p, ctypes.c_size_t], None),
    ("uperm", [c_double_p, c_double_p, c_size_t_p, ctypes.c_size_t], None),
    ("invmatperm", [c_double_p, c_double_p, c_size_t_p, ctypes.c_size_t], None),
]:
    _fn = getattr(_lib, _name)
    _fn.argtypes = _argtypes
    _fn.restype = _restype


def _real_sdp_blocks(cK):
    """Yields (offset_into_x_or_frms_style_array, nk) for each REAL
    (non-complex-Hermitian) PSD block, i.e. k in range(cK.rsdpN). Complex
    Hermitian PSD blocks (k in rsdpN:sdpN) are not covered by any
    function in this cluster yet -- a documented gap, not an oversight;
    see each function's docstring."""
    sdpNL = cK._keepalive[2]
    return [int(sdpNL[k]) for k in range(cK.rsdpN)]


def givensrot(gjc, g, x, K: dict):
    """y = givensrot(gjc, g, x, K): apply a sequence of precomputed
    Givens rotations to each real PSD block of x. Real-symmetric blocks
    only (K.s[:K.rsdpN]) -- complex Hermitian is not covered. Wraps
    matgivens() (givensrot.c) per block."""
    import numpy as np

    cK = cone_from_dict(K)
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    y = x.copy()
    gjc_full = np.ascontiguousarray(gjc, dtype=np.int64).ravel()  # already
    # 0-indexed/C-style per givensrot.c's comment ("don't subtract 1")
    g_full = np.ascontiguousarray(g, dtype=np.float64).ravel()

    xoff = goff = joff = 0
    for nk in _real_sdp_blocks(cK):
        nksqr = nk * nk
        gjc_blk = np.ascontiguousarray(gjc_full[joff : joff + nk], dtype=np.uintp)
        _lib.matgivens(
            y[xoff : xoff + nksqr].ctypes.data_as(c_double_p),
            g_full[goff:].ctypes.data_as(c_double_p),
            gjc_blk.ctypes.data_as(c_size_t_p), nk,
        )
        goff += 2 * int(gjc_blk[nk - 1]) if nk > 0 else 0
        xoff += nksqr
        joff += nk
    return y


def urotorder(u, K: dict, maxu: float, perm_in=None):
    """[u_out, perm, gjc, g] = urotorder(u, K, maxu, perm_in=None):
    stably reorders each real PSD block's upper-triangular factor via
    Givens rotations. Real-symmetric blocks only. Wraps rotorder()
    (urotorder.c) per block, then uperm()+triu2sym() to physically
    permute and symmetrize.

    `perm_in`, if given (matching urotorder.m's optional 4th argument),
    is composed with the freshly computed per-block permutation instead
    of just converting it to 1-indexed
    (`perm_out[i] = perm_in[perm[i]]` per urotorder.c's own
    `permPr[i] = permOld[perm[i]]`) -- done here in Python rather than
    by extending the C call, since the underlying rotation math is
    identical either way; only how the output `perm` is labeled changes.
    """
    import numpy as np

    cK = cone_from_dict(K)
    u_in = np.ascontiguousarray(u, dtype=np.float64).ravel()
    u_out = np.empty_like(u_in)
    perm_out = np.zeros(cK.rLen + cK.hLen, dtype=np.float64)
    gjc_out = np.zeros(cK.rLen + cK.hLen, dtype=np.float64)
    g_chunks = []
    perm_in_arr = None
    if perm_in is not None and np.asarray(perm_in).size:
        perm_in_arr = np.ascontiguousarray(perm_in, dtype=np.int64).ravel()

    xoff = poff = 0
    for nk in _real_sdp_blocks(cK):
        nksqr = nk * nk
        fwork = u_in[xoff : xoff + nksqr].copy()
        perm = np.zeros(nk, dtype=np.uintp)
        gjc = np.zeros(nk, dtype=np.uintp)
        d = np.zeros(nk, dtype=np.float64)
        g = np.zeros(max(nk * (nk - 1), 1), dtype=np.float64)  # upper bound

        _lib.rotorder(
            perm.ctypes.data_as(c_size_t_p), fwork.ctypes.data_as(c_double_p),
            gjc.ctypes.data_as(c_size_t_p), g.ctypes.data_as(c_double_p),
            d.ctypes.data_as(c_double_p), float(maxu) ** 2, nk,
        )
        block_out = np.empty(nksqr, dtype=np.float64)
        _lib.uperm(
            block_out.ctypes.data_as(c_double_p), fwork.ctypes.data_as(c_double_p),
            perm.ctypes.data_as(c_size_t_p), nk,
        )
        _lib.triu2sym(block_out.ctypes.data_as(c_double_p), nk)
        u_out[xoff : xoff + nksqr] = block_out

        ginz = int(gjc[nk - 1]) if nk > 0 else 0
        if perm_in_arr is not None:
            perm_out[poff : poff + nk] = perm_in_arr[poff : poff + nk][perm.astype(np.int64)]
        else:
            perm_out[poff : poff + nk] = perm + 1  # 1-indexed, .m-facing
        gjc_out[poff : poff + nk] = gjc
        g_chunks.append(g[: 2 * ginz])
        xoff += nksqr
        poff += nk
    g_out = np.concatenate(g_chunks) if g_chunks else np.zeros(0)
    return u_out, perm_out, gjc_out, g_out


def sqrtinv(q, vlab, K: dict):
    """y = sqrtinv(q, vlab, K): y = (Q / diag(sqrt(vlab)))' per real PSD
    block, so Y'*Y = inv(Q*diag(vlab)*Q'). Real-symmetric blocks only.
    Wraps qdivv() (sqrtinv.c) per block."""
    import numpy as np

    cK = cone_from_dict(K)
    q = np.ascontiguousarray(q, dtype=np.float64).ravel()
    vlab = np.ascontiguousarray(vlab, dtype=np.float64).ravel()
    diagskip = cK.lpN + 2 * cK.lorN
    v = vlab[diagskip:]
    y = np.empty_like(q)

    qoff = voff = 0
    for nk in _real_sdp_blocks(cK):
        nksqr = nk * nk
        _lib.qdivv(
            y[qoff : qoff + nksqr].ctypes.data_as(c_double_p),
            q[qoff : qoff + nksqr].ctypes.data_as(c_double_p),
            v[voff : voff + nk].ctypes.data_as(c_double_p), nk,
        )
        qoff += nksqr
        voff += nk
    return y


def psdframeit(lab, frms, K: dict):
    """x = psdframeit(lab, frms, K): x = FRM*lab, FRM a product-form
    Householder reflection (as produced by qrK()). Real-symmetric blocks
    only. Wraps psdframeit() (psdframeit.c) directly -- it already loops
    over every block itself, given sdpNL/rsdpN/sdpN."""
    import numpy as np

    cK = cone_from_dict(K)
    lab = np.ascontiguousarray(lab, dtype=np.float64).ravel()
    frms = np.ascontiguousarray(frms, dtype=np.float64).ravel()
    lenud = cK.rDim + cK.hDim
    x = np.zeros(lenud, dtype=np.float64)
    fwsiz = max(cK.rMaxn**2, 2 * cK.hMaxn**2, 1)
    fwork = np.zeros(fwsiz, dtype=np.float64)
    sdpNL = np.ascontiguousarray(cK._keepalive[2], dtype=np.uintp)

    _lib.psdframeit(
        x.ctypes.data_as(c_double_p), frms.ctypes.data_as(c_double_p),
        lab.ctypes.data_as(c_double_p), sdpNL.ctypes.data_as(c_size_t_p),
        cK.rsdpN, cK.sdpN, fwork.ctypes.data_as(c_double_p),
    )
    return x


def psdinvjmul(x, frms, y, K: dict):
    """z = psdinvjmul(x, frms, y, K): solves X*Z+Z*X = 2*Y in the PSD
    cone, given X's eigenvalues `x` and eigenbasis `frms` (product-form
    Householder, as from qrK()). Real-symmetric blocks only. Wraps
    psdinvjmul() (psdinvjmul.c) directly (loops over blocks itself).

    `x` and `y` may each be given either as exactly the PSD-only length
    (cK.rLen+cK.hLen for x, cK.rDim+cK.hDim for y) or as a full internal
    vector (L+Q+S) -- matching psdinvjmul.c's own mexFunction, which
    auto-detects a full-length input and skips its L+Q prefix
    (`x += cK.lpN + 2*cK.lorN`, `y += cK.lpN + cK.qDim`) rather than
    silently misreading a full vector as if it started at the PSD block
    (the bug this replicated fix corrects: found while porting
    wregion.m, which always calls this with full-length vTAR/dxmdz)."""
    import numpy as np

    cK = cone_from_dict(K)
    frms = np.ascontiguousarray(frms, dtype=np.float64).ravel()
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y, dtype=np.float64).ravel()
    lenud = cK.rDim + cK.hDim
    if x.size != cK.rLen + cK.hLen:
        x = x[cK.lpN + 2 * cK.lorN :]
    if y.size != lenud:
        y = y[cK.lpN + cK.qDim :]
    z = np.zeros(lenud, dtype=np.float64)
    fwsiz = max(cK.rMaxn, 2 * cK.hMaxn, 1)
    fwork = np.zeros(fwsiz, dtype=np.float64)
    sdpNL = np.ascontiguousarray(cK._keepalive[2], dtype=np.uintp)

    _lib.psdinvjmul(
        z.ctypes.data_as(c_double_p), frms.ctypes.data_as(c_double_p),
        x.ctypes.data_as(c_double_p), y.ctypes.data_as(c_double_p),
        sdpNL.ctypes.data_as(c_size_t_p), cK.rsdpN, cK.sdpN,
        fwork.ctypes.data_as(c_double_p),
    )
    return z


def qrK(x, K: dict):
    """[frms, r] = qrK(x, K): QR-factorizes each real PSD block of x via
    Householder reflections. Real-symmetric blocks only. Wraps qrfac()
    (qrK.c) per block.

    Per qrK.c's mexFunction, a real nxn block of `frms` is qrfac()'s `q`
    and `beta` outputs packed into ONE nksqr(=n^2)-length buffer: the
    first n*(n-1) values are the Householder vectors (qrfac's `q`,
    exactly as psdframeit()/psdinvjmul() expect them), and the LAST n
    values are `beta` (qrfac is called with beta pointing at
    frms_block[nksqr-n:]. `r`'s block is qrfac's mutated `u` verbatim --
    triu(u) is the real upper-triangular R factor, but tril(u,-1) is
    left as whatever qrfac's arithmetic happened to leave there
    ("undefined" per qrfac's own doc comment) rather than cleaned to
    zero, matching the real MEX build bit-for-bit.
    """
    import numpy as np

    cK = cone_from_dict(K)
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    lenud = cK.rDim + cK.hDim
    if x.size != lenud:
        raise ValueError(f"x must have length {lenud}, got {x.size}")

    frms = np.zeros(lenud, dtype=np.float64)
    r = x.copy()

    off = 0
    for nk in _real_sdp_blocks(cK):
        nksqr = nk * nk
        u = r[off : off + nksqr]  # qrfac mutates this block of r in place
        frms_blk = frms[off : off + nksqr]
        beta_ptr = frms_blk[nksqr - nk :].ctypes.data_as(c_double_p)
        _lib.qrfac(
            beta_ptr, frms_blk.ctypes.data_as(c_double_p),
            u.ctypes.data_as(c_double_p), nk,
        )
        off += nksqr
    return frms, r


def invcholfac(u, K: dict, perm=None):
    """y = invcholfac(u, K, perm=None): y = U'*U per real PSD block (or
    invperm(U'*U) if `perm` given). Real-symmetric blocks only. Wraps
    utmulx()+triu2sym() (+invmatperm() if perm given), from
    invcholfac.c/triuaux.c."""
    import numpy as np

    cK = cone_from_dict(K)
    u = np.ascontiguousarray(u, dtype=np.float64).ravel()
    lenud = cK.rDim + cK.hDim
    y = np.zeros(lenud, dtype=np.float64)

    perm_arr = None
    if perm is not None:
        perm_arr = (np.ascontiguousarray(perm, dtype=np.int64).ravel() - 1).astype(np.uintp)

    off = poff = 0
    for nk in _real_sdp_blocks(cK):
        nksqr = nk * nk
        block = np.zeros(nksqr, dtype=np.float64)
        _lib.utmulx(
            block.ctypes.data_as(c_double_p), u[off : off + nksqr].ctypes.data_as(c_double_p),
            u[off : off + nksqr].ctypes.data_as(c_double_p), nk,
        )
        _lib.triu2sym(block.ctypes.data_as(c_double_p), nk)
        if perm_arr is not None:
            _lib.invmatperm(
                y[off : off + nksqr].ctypes.data_as(c_double_p),
                block.ctypes.data_as(c_double_p),
                perm_arr[poff : poff + nk].ctypes.data_as(c_size_t_p), nk,
            )
        else:
            y[off : off + nksqr] = block
        off += nksqr
        poff += nk
    return y


# extractA()/adendotd() take a `jcir` struct ({double *pr; mwIndex *jc,
# *ir;}) BY VALUE (not by pointer) -- ctypes marshals a ctypes.Structure
# argument per the platform ABI automatically, so this just needs the
# matching field layout, not any special handling on the call site.
class JcIr(ctypes.Structure):
    _fields_ = [("pr", c_double_p), ("jc", c_size_t_p), ("ir", c_size_t_p)]


_lib.extractA.argtypes = [JcIr, c_size_t_p, c_size_t_p, c_size_t_p, c_double_p,
                            ctypes.c_size_t, ctypes.c_size_t]
_lib.extractA.restype = None

_lib.findblks.argtypes = [
    c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
    c_size_t_p, c_size_t_p, ctypes.c_size_t, ctypes.c_size_t,
    ctypes.c_size_t, c_ubyte_p, c_size_t_p,
]
_lib.findblks.restype = None

_lib.partitA.argtypes = [
    c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
    ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, c_ubyte_p, c_size_t_p,
]
_lib.partitA.restype = None

_lib.adendotd.argtypes = [
    JcIr, JcIr, JcIr, c_double_p, c_double_p, c_size_t_p, c_size_t_p,
    c_size_t_p, ctypes.c_size_t, ctypes.c_size_t, c_double_p,
]
_lib.adendotd.restype = None

_lib.adenscale.argtypes = [
    c_double_p, c_double_p, c_size_t_p, c_size_t_p, c_size_t_p,
    ctypes.c_size_t, ctypes.c_size_t,
]
_lib.adenscale.restype = None


def adenscale(dense: dict, d: dict, qblkstart):
    """smult = adenscale(dense, d, qblkstart): length-nden vector with
    smult[j] = det(dk) for the Lorentz block k that dense column j (a
    dense Lorentz norm-bound column) belongs to -- the Woodbury-update
    scale factor s.t. AP(d)A' = ADA + Ad*diag(smult)*Ad'. `dense` needs
    "l"/"q"/"cols" (see getdense.py; `dense["cols"][dense["l"]+len(dense["q"]):]`
    is the dense Lorentz norm-bound column subset this operates on).
    `d` needs "det" (length len(K["q"]), see sdinit.py/updtransfo.py).
    `qblkstart` is K["qblkstart"] (1-indexed cumulative boundaries).
    Wraps adenscale() (adenscale.c) exactly as deninfac.m calls it.
    """
    import numpy as np

    qblkstart_arr = np.ascontiguousarray(qblkstart, dtype=np.int64).ravel()

    nl = int(dense["l"])
    q_field = np.ascontiguousarray(dense["q"], dtype=np.int64).ravel()
    nq = q_field.size
    cols = np.ascontiguousarray(dense["cols"], dtype=np.int64).ravel()
    nden = cols.size - nl - nq

    dencols0 = (cols[nl + nq :] - 1).astype(np.uintp)
    q0 = (q_field - 1).astype(np.int64)
    blkend0 = (qblkstart_arr[q0 + 1] - 1).astype(np.uintp)
    q0 = q0.astype(np.uintp)

    detd = np.ascontiguousarray(d["det"], dtype=np.float64).ravel()
    smult = np.zeros(max(nden, 1), dtype=np.float64)

    _lib.adenscale(
        smult.ctypes.data_as(c_double_p), detd.ctypes.data_as(c_double_p),
        dencols0.ctypes.data_as(c_size_t_p), q0.ctypes.data_as(c_size_t_p),
        blkend0.ctypes.data_as(c_size_t_p), nq, nden,
    )
    return smult[:nden]


def adendotd(dense: dict, d: dict, adotd, ablk, qblkstart):
    """ad = adendotd(dense, d, adotd, ablk, qblkstart): (ai[k]+Adeni[k])'*
    d[k] for each dense Lorentz block k in dense["q"] -- adds the
    dense-column contribution (dense["A"]'s trace + norm-bound columns)
    on top of `adotd`'s sparse-part contribution, filled into `ablk`'s
    sparsity pattern (a fresh m x nq matrix with that pattern; `ablk`'s
    own data is ignored/overwritten, only its pattern is reused, matching
    the mex wrapper's `mxDuplicateArray` semantics). `dense` additionally
    needs "A" (sparse m x (nl+nq+nden), see getdense.py: first nl columns
    LP-dense, next nq Lorentz-trace, remaining nden Lorentz-norm-bound).
    `d` needs "q1"/"q2" (see sdinit.py/updtransfo.py; "q2" is already
    0-indexed relative to `qblkstart[0]-1`, the same convention used
    throughout cone.py/updtransfo.py). `adotd` is a sparse m x nq matrix
    (DAt["q"] restricted to the dense["q"] columns). Wraps adendotd()
    (adendotd.c) exactly as getDAtm.m calls it.
    """
    import numpy as np
    import scipy.sparse as sp

    qblkstart_arr = np.ascontiguousarray(qblkstart, dtype=np.int64).ravel()
    firstQ = int(qblkstart_arr[0]) - 1

    nl = int(dense["l"])
    q_field = np.ascontiguousarray(dense["q"], dtype=np.int64).ravel()
    nq = q_field.size
    cols = np.ascontiguousarray(dense["cols"], dtype=np.int64).ravel()
    nden = cols.size - nl - nq

    dencols0 = (cols[nl + nq :] - 1).astype(np.uintp)
    q0 = (q_field - 1).astype(np.int64)
    blkend0 = (qblkstart_arr[q0 + 1] - 1).astype(np.uintp)
    q0 = q0.astype(np.uintp)

    A = dense["A"].tocsc()
    m = A.shape[0]
    aden = A[:, nl : nl + nq + nden].tocsc()
    aden_pr = np.ascontiguousarray(aden.data, dtype=np.float64)
    aden_jc = np.ascontiguousarray(aden.indptr, dtype=np.uintp)
    aden_ir = np.ascontiguousarray(aden.indices, dtype=np.uintp)
    aden_struct = JcIr(pr=aden_pr.ctypes.data_as(c_double_p), jc=aden_jc.ctypes.data_as(c_size_t_p),
                        ir=aden_ir.ctypes.data_as(c_size_t_p))

    ADOTD = adotd.tocsc() if sp.issparse(adotd) else sp.csc_matrix(adotd)
    adotd_pr = np.ascontiguousarray(ADOTD.data, dtype=np.float64)
    adotd_jc = np.ascontiguousarray(ADOTD.indptr, dtype=np.uintp)
    adotd_ir = np.ascontiguousarray(ADOTD.indices, dtype=np.uintp)
    adotd_struct = JcIr(pr=adotd_pr.ctypes.data_as(c_double_p), jc=adotd_jc.ctypes.data_as(c_size_t_p),
                         ir=adotd_ir.ctypes.data_as(c_size_t_p))

    ABLK = ablk.tocsc() if sp.issparse(ablk) else sp.csc_matrix(ablk)
    ad_pr = np.zeros(max(int(ABLK.indptr[-1]), 1), dtype=np.float64)
    ad_jc = np.ascontiguousarray(ABLK.indptr, dtype=np.uintp)
    ad_ir = np.ascontiguousarray(ABLK.indices, dtype=np.uintp)
    ad_struct = JcIr(pr=ad_pr.ctypes.data_as(c_double_p), jc=ad_jc.ctypes.data_as(c_size_t_p),
                      ir=ad_ir.ctypes.data_as(c_size_t_p))

    d1 = np.ascontiguousarray(d["q1"], dtype=np.float64).ravel()
    # d.c's mexFunction passes `d2 - firstQ` (pointer arithmetic) so that
    # C's `d2[i]` for a GLOBAL 0-indexed row i reads d["q2"][i-firstQ];
    # ctypes has no clean way to express a negative-offset pointer, so
    # instead zero-pad the front of the buffer by firstQ entries and keep
    # `dencols0` as global (unshifted) subscripts -- d2_padded[i] then
    # equals d["q2"][i-firstQ] for every i>=firstQ, exactly reproducing
    # the C pointer's addressing without needing pointer arithmetic.
    d2 = np.ascontiguousarray(d["q2"], dtype=np.float64).ravel()
    d2_padded = np.concatenate([np.zeros(firstQ, dtype=np.float64), d2])

    nnz = int(ABLK.indptr[-1])
    fwork = np.zeros(max(m, 1), dtype=np.float64)

    _lib.adendotd(
        ad_struct, adotd_struct, aden_struct,
        d1.ctypes.data_as(c_double_p), d2_padded.ctypes.data_as(c_double_p),
        q0.ctypes.data_as(c_size_t_p), dencols0.ctypes.data_as(c_size_t_p),
        blkend0.ctypes.data_as(c_size_t_p), nq, nden, fwork.ctypes.data_as(c_double_p),
    )
    return sp.csc_matrix(
        (ad_pr[:nnz], ad_ir[:nnz].astype(np.int64), ad_jc.astype(np.int64)), shape=(m, nq)
    )


def extractA(At, Ajc_table, blk0: int, blk1, blkstart):
    """Apart = extractA(At, Ajc, blk0, blk1, blkstart): fast row-range
    slice of a sparse matrix At (m columns, transposed convention as
    SeDuMi uses throughout), restricted per-column to the nonzero range
    Ajc_table[:, blk0-1] .. Ajc_table[:, blk1-1] (1-indexed block
    columns; blk0<=0 means "from the start of the column", blk1 empty/
    None means "to the end"). blkstart = (ifirst, n) 1-indexed row range
    (MATLAB's 2-element blkstart form, not the 6-arg blkstart2 form).
    Wraps extractA() (extractA.c) directly."""
    import numpy as np
    import scipy.sparse

    A = At.tocsc()
    m = A.shape[1]
    ifirst, n_end = int(blkstart[0]), int(blkstart[1])
    ifirst -= 1
    n = n_end - (ifirst + 1)

    Ajc_table = None if Ajc_table is None else np.ascontiguousarray(Ajc_table, dtype=np.int64)
    njc = 0 if Ajc_table is None else Ajc_table.shape[1]

    Ajc = np.zeros(2 * m, dtype=np.uintp)
    if blk0 <= 0 or Ajc_table is None:
        Ajc[:m] = A.indptr[:m]
    else:
        Ajc[:m] = Ajc_table[:, blk0 - 1]
    if blk1 is None or (njc and blk1 - 1 >= njc):
        Ajc[m:] = A.indptr[1:]
    else:
        Ajc[m:] = Ajc_table[:, blk1 - 1]

    ynnz = int(np.sum(Ajc[m:].astype(np.int64) - Ajc[:m].astype(np.int64)))
    Y_jc = np.zeros(m + 1, dtype=np.uintp)
    Y_ir = np.zeros(max(ynnz, 1), dtype=np.uintp)
    Y_pr = np.zeros(max(ynnz, 1), dtype=np.float64)

    Air = np.ascontiguousarray(A.indices, dtype=np.uintp)
    Apr = np.ascontiguousarray(A.data, dtype=np.float64)
    Y = JcIr(pr=Y_pr.ctypes.data_as(c_double_p), jc=Y_jc.ctypes.data_as(c_size_t_p),
             ir=Y_ir.ctypes.data_as(c_size_t_p))

    _lib.extractA(
        Y, Ajc[:m].ctypes.data_as(c_size_t_p), Ajc[m:].ctypes.data_as(c_size_t_p),
        Air.ctypes.data_as(c_size_t_p), Apr.ctypes.data_as(c_double_p), ifirst, m,
    )
    return scipy.sparse.csc_matrix((Y_pr, Y_ir.astype(np.int64), Y_jc.astype(np.int64)), shape=(n, m))


def findblks(At, Ablkjc_table, blk0: int, blk1, blkstart):
    """Ablk = findblks(At, Ablkjc, blk0, blk1, blkstart): sparse nblk x m
    0/1 indicator of which of the nblk row-ranges (given by 1-indexed
    blkstart, length nblk+1) each column of At has a nonzero in (within
    the blk0..blk1 restricted nnz range, same convention as extractA()).
    Wraps findblks() (findblks.c) directly."""
    import numpy as np
    import scipy.sparse

    A = At.tocsc()
    m = A.shape[1]
    # findblks.c's mexFunction decrements each blkstart(i) TWICE (asserting
    # positivity after each): blkstart[i] = blkstartPr[i]-1, and separately
    # blkstart[nblk+i] = blkstartPr[i]-2 -- the SAME source value both
    # times, not the next boundary. (blkstart here therefore needs every
    # entry >= 2 in 1-indexed terms; a leading boundary of 1 fails the
    # real MEX build's own assertion, not just this port.)
    blkstart_arr = np.ascontiguousarray(blkstart, dtype=np.int64)
    nblk = blkstart_arr.size - 1
    blkstart0 = np.zeros(2 * nblk, dtype=np.uintp)
    blkstart0[:nblk] = blkstart_arr[:nblk] - 1
    blkstart0[nblk:] = blkstart_arr[:nblk] - 2

    Ablkjc_table = None if Ablkjc_table is None else np.ascontiguousarray(Ablkjc_table, dtype=np.int64)
    njc = 0 if Ablkjc_table is None else Ablkjc_table.shape[1]

    Ajc = np.zeros(2 * m, dtype=np.uintp)
    if blk0 <= 0 or Ablkjc_table is None:
        Ajc[:m] = A.indptr[:m]
    else:
        Ajc[:m] = Ablkjc_table[:, blk0 - 1]
    if blk1 is None or (njc and blk1 - 1 >= njc):
        Ajc[m:] = A.indptr[1:]
    else:
        Ajc[m:] = Ablkjc_table[:, blk1 - 1]

    blknnz = max(int(np.sum(Ajc[m:].astype(np.int64) - Ajc[:m].astype(np.int64))), 1)
    Ablk_jc = np.zeros(m + 1, dtype=np.uintp)
    Ablk_ir = np.zeros(blknnz, dtype=np.uintp)
    iwsize = nblk + 2 + int(np.floor(np.log(1.0 + nblk) / np.log(2.0)))
    iwork = np.zeros(max(iwsize, 1), dtype=np.uintp)
    cwork = np.zeros(max(nblk, 1), dtype=np.uint8)
    Air = np.ascontiguousarray(A.indices, dtype=np.uintp)

    _lib.findblks(
        Ablk_ir.ctypes.data_as(c_size_t_p), Ablk_jc.ctypes.data_as(c_size_t_p),
        Ajc[:m].ctypes.data_as(c_size_t_p), Ajc[m:].ctypes.data_as(c_size_t_p),
        Air.ctypes.data_as(c_size_t_p),
        blkstart0[:nblk].ctypes.data_as(c_size_t_p), blkstart0[nblk:].ctypes.data_as(c_size_t_p),
        m, nblk, iwsize, cwork.ctypes.data_as(c_ubyte_p), iwork.ctypes.data_as(c_size_t_p),
    )
    nnz_final = int(Ablk_jc[m])
    data = np.ones(max(nnz_final, 1), dtype=np.float64)[:nnz_final]
    return scipy.sparse.csc_matrix(
        (data, Ablk_ir[:nnz_final].astype(np.int64), Ablk_jc.astype(np.int64)), shape=(nblk, m)
    )


def partitA(At, blkstart):
    """Ablkjc = partitA(At, blkstart): m x nblk table where column k
    gives, for each column j of At, the first nonzero-row subscript at
    or beyond blkstart[k] (1-indexed blkstart, length nblk). This is
    exactly the "Ajc table" extractA()/findblks() take as their
    Ajc_table/Ablkjc_table argument. Wraps partitA() (partitA.c)
    directly."""
    import numpy as np

    A = At.tocsc()
    m = A.shape[1]
    blkstart_arr = (np.ascontiguousarray(blkstart, dtype=np.int64) - 1).astype(np.uintp)
    nblk = blkstart_arr.size
    L = nblk + 2

    iwsize = max(int(np.floor(np.log(1.0 + nblk) / np.log(2.0))), 0)
    iwork = np.zeros(max(iwsize, 1), dtype=np.uintp)
    Ablkjc_work = np.zeros(L * m, dtype=np.uintp)
    cwork = np.zeros(max(nblk, 1), dtype=np.uint8)
    Ajc = np.ascontiguousarray(A.indptr, dtype=np.uintp)
    Air = np.ascontiguousarray(A.indices, dtype=np.uintp)

    _lib.partitA(
        Ablkjc_work.ctypes.data_as(c_size_t_p), Ajc.ctypes.data_as(c_size_t_p),
        Air.ctypes.data_as(c_size_t_p), blkstart_arr.ctypes.data_as(c_size_t_p),
        m, nblk, iwsize, cwork.ctypes.data_as(c_ubyte_p), iwork.ctypes.data_as(c_size_t_p),
    )
    out = np.empty((m, nblk), dtype=np.int64)
    Ablkjc_work_r = Ablkjc_work.reshape(m, L)
    out[:, :] = Ablkjc_work_r[:, 1 : nblk + 1].astype(np.int64)
    return out


_lib.getada1.argtypes = [
    JcIr, JcIr, c_double_p, c_double_p, c_size_t_p, c_size_t_p,
    c_size_t_p, c_size_t_p, ctypes.c_size_t, ctypes.c_size_t, c_double_p,
]
_lib.getada1.restype = None

_lib.getada2.argtypes = [
    JcIr, JcIr, c_size_t_p, c_size_t_p, ctypes.c_size_t, ctypes.c_size_t, c_double_p,
]
_lib.getada2.restype = None

_lib.getada3.argtypes = [
    JcIr, c_double_p, JcIr, c_double_p,
    c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
    c_size_t_p, c_size_t_p, c_size_t_p,
    c_size_t_p, c_size_t_p,
    ctypes.c_size_t, ctypes.c_size_t, ctypes.POINTER(ConeK),
    c_double_p, ctypes.c_size_t, c_size_t_p, ctypes.c_size_t,
    c_ubyte_p,
]
_lib.getada3.restype = None

_lib.dzblkpartit.argtypes = [
    c_size_t_p, c_size_t_p, c_size_t_p, ctypes.c_size_t, ctypes.c_size_t,
]
_lib.dzblkpartit.restype = None

_lib.spmakesym.argtypes = [JcIr, ctypes.c_size_t, c_size_t_p]
_lib.spmakesym.restype = None

_lib.cpspdiag.argtypes = [c_double_p, JcIr, ctypes.c_size_t]
_lib.cpspdiag.restype = None


def getada1(ada, At, Ajc2, perm, d: dict, qblkstart):
    """ADA = getada1(ADA, At, Ajc2, perm, d, qblkstart): fresh
    ADA(i,j) = (D(d^2; LP,Lorentz)*At(:,i))'*At(:,j) on
    triu(ADA(perm,perm)) only (matching the sparsity pattern of `ada`;
    entries outside that triangle come back as 0 regardless of what
    `ada` held there). Wraps getada1() (getada1.c) exactly as
    sedumi.m's main loop calls it: `Ajc2` is `Ablkjc[:, 2]` (0-indexed,
    partitA()'s 3rd column), `perm` is `Aord["lqperm"]` (1-indexed,
    sortnnz()'s output convention), `qblkstart` is `K["qblkstart"]`
    (1-indexed cumulative boundaries, length 1+len(K["q"])). `d` needs
    "l" and "det" arrays (see sdinit.py/updtransfo.py).
    """
    import numpy as np
    import scipy.sparse as sp

    A = At.tocsc()
    m = A.shape[1]

    qblkstart_arr = np.ascontiguousarray(qblkstart, dtype=np.int64).ravel()
    nblk = qblkstart_arr.size
    d_l = np.ascontiguousarray(d["l"], dtype=np.float64).ravel()
    d_det = np.ascontiguousarray(d["det"], dtype=np.float64).ravel()

    blkstart = np.zeros(nblk + 1, dtype=np.uintp)
    blkstart[0] = d_l.size
    blkstart[1:] = (qblkstart_arr - 1).astype(np.uintp)

    perm_arr = (np.ascontiguousarray(perm, dtype=np.int64).ravel() - 1).astype(np.uintp)
    invperm = np.zeros(m, dtype=np.uintp)
    invperm[perm_arr] = np.arange(m, dtype=np.uintp)

    Ajc2_arr = np.ascontiguousarray(Ajc2, dtype=np.uintp).ravel()

    Apr = np.ascontiguousarray(A.data, dtype=np.float64)
    Ajc = np.ascontiguousarray(A.indptr, dtype=np.uintp)
    Air = np.ascontiguousarray(A.indices, dtype=np.uintp)
    at_struct = JcIr(pr=Apr.ctypes.data_as(c_double_p), jc=Ajc.ctypes.data_as(c_size_t_p),
                      ir=Air.ctypes.data_as(c_size_t_p))

    pattern = ada.tocsc()
    ada_jc = np.ascontiguousarray(pattern.indptr, dtype=np.uintp)
    ada_ir = np.ascontiguousarray(pattern.indices, dtype=np.uintp)
    nnz = int(pattern.indptr[m])
    ada_pr = np.zeros(max(nnz, 1), dtype=np.float64)
    ada_struct = JcIr(pr=ada_pr.ctypes.data_as(c_double_p), jc=ada_jc.ctypes.data_as(c_size_t_p),
                       ir=ada_ir.ctypes.data_as(c_size_t_p))

    fwsiz = max(2 * int(blkstart[nblk]), 1)
    fwork = np.zeros(fwsiz, dtype=np.float64)

    _lib.getada1(
        ada_struct, at_struct, d_l.ctypes.data_as(c_double_p), d_det.ctypes.data_as(c_double_p),
        Ajc2_arr.ctypes.data_as(c_size_t_p), blkstart.ctypes.data_as(c_size_t_p),
        perm_arr.ctypes.data_as(c_size_t_p), invperm.ctypes.data_as(c_size_t_p),
        m, nblk, fwork.ctypes.data_as(c_double_p),
    )
    return sp.csc_matrix(
        (ada_pr[:nnz], ada_ir[:nnz].astype(np.int64), ada_jc.astype(np.int64)), shape=(m, m)
    )


def getada2(ada, DAt: dict, Aord: dict, K: dict):
    """ADA = getada2(ADA, DAt, Aord, K): ADA += DAt["q"]'*DAt["q"] on
    triu(ADA(qperm,qperm)) only; a plain copy of `ada` (unchanged) when
    K has no Lorentz blocks or DAt["q"] is completely empty -- both
    short-circuits present in the real mex wrapper. Wraps getada2()
    (getada2.c) exactly as sedumi.m's main loop calls it.
    """
    import numpy as np
    import scipy.sparse as sp

    m = ada.shape[0]
    pattern = ada.tocsc()
    ada_jc = np.ascontiguousarray(pattern.indptr, dtype=np.uintp)
    ada_ir = np.ascontiguousarray(pattern.indices, dtype=np.uintp)
    nnz = int(pattern.indptr[m])
    ada_pr = np.array(pattern.data, dtype=np.float64, copy=True)

    def _wrap():
        return sp.csc_matrix(
            (ada_pr[:nnz], ada_ir[:nnz].astype(np.int64), ada_jc.astype(np.int64)), shape=(m, m)
        )

    lorN = len(K.get("q", []))
    ddota = DAt.get("q")
    if lorN <= 0 or ddota is None:
        return _wrap()

    ddota = ddota.tocsc() if sp.issparse(ddota) else sp.csc_matrix(ddota)
    if int(ddota.indptr[m]) == 0:
        return _wrap()

    perm_arr = (np.ascontiguousarray(Aord["qperm"], dtype=np.int64).ravel() - 1).astype(np.uintp)
    invperm = np.zeros(m, dtype=np.uintp)
    invperm[perm_arr] = np.arange(m, dtype=np.uintp)

    dd_pr = np.ascontiguousarray(ddota.data, dtype=np.float64)
    dd_jc = np.ascontiguousarray(ddota.indptr, dtype=np.uintp)
    dd_ir = np.ascontiguousarray(ddota.indices, dtype=np.uintp)
    ddota_struct = JcIr(pr=dd_pr.ctypes.data_as(c_double_p), jc=dd_jc.ctypes.data_as(c_size_t_p),
                         ir=dd_ir.ctypes.data_as(c_size_t_p))
    ada_struct = JcIr(pr=ada_pr.ctypes.data_as(c_double_p), jc=ada_jc.ctypes.data_as(c_size_t_p),
                       ir=ada_ir.ctypes.data_as(c_size_t_p))

    fwork = np.zeros(max(lorN, 1), dtype=np.float64)

    _lib.getada2(
        ada_struct, ddota_struct,
        perm_arr.ctypes.data_as(c_size_t_p), invperm.ctypes.data_as(c_size_t_p),
        m, lorN, fwork.ctypes.data_as(c_double_p),
    )
    return _wrap()


def getada3(ada, At, Ajc1, Aord: dict, udsqr, K: dict):
    """(ADA, absd) = getada3(ADA, At, Ajc1, Aord, udsqr, K): ADA(i,j) +=
    ai'*D(d^2;PSD)*aj on triu(ADA(sperm,sperm)) (a no-op when K has no
    PSD blocks -- `absd` is then `diag(ADA)` via cpspdiag() instead of
    the PSD cancellation measure), followed unconditionally by
    ADA := (ADA+ADA')/2 via spmakesym(). `Ajc1` is `Ablkjc[:, 2]`
    (0-indexed, partitA()'s 3rd column -- the same array getada1() calls
    `Ajc2`). `Aord` needs "dz" (sparse lenfull x m, incorder()'s dz
    output) and "sperm" (1-indexed, incorder()'s perm output). `udsqr`
    is `_native.invcholfac(d["u"], K, d["perm"])`. Wraps getada3()
    (getada3.c), dzblkpartit(), cpspdiag(), and spmakesym() exactly as
    sedumi.m's main loop / getada3.m's mex wrapper does.
    """
    import numpy as np
    import scipy.sparse as sp

    cK = cone_from_dict(K)
    m = int(At.shape[1])
    lenud = int(cK.rDim) + int(cK.hDim)

    A = At.tocsc()
    Apr = np.ascontiguousarray(A.data, dtype=np.float64)
    Ajc = np.ascontiguousarray(A.indptr, dtype=np.uintp)
    Air = np.ascontiguousarray(A.indices, dtype=np.uintp)
    at_struct = JcIr(pr=Apr.ctypes.data_as(c_double_p), jc=Ajc.ctypes.data_as(c_size_t_p),
                      ir=Air.ctypes.data_as(c_size_t_p))

    udsqr_arr = np.ascontiguousarray(udsqr, dtype=np.float64).ravel()
    Ajc1_arr = np.ascontiguousarray(Ajc1, dtype=np.uintp).ravel()

    pattern = ada.tocsc()
    ada_jc = np.ascontiguousarray(pattern.indptr, dtype=np.uintp)
    ada_ir = np.ascontiguousarray(pattern.indices, dtype=np.uintp)
    nnz = int(pattern.indptr[m])
    ada_pr = np.array(pattern.data, dtype=np.float64, copy=True)
    ada_struct = JcIr(pr=ada_pr.ctypes.data_as(c_double_p), jc=ada_jc.ctypes.data_as(c_size_t_p),
                       ir=ada_ir.ctypes.data_as(c_size_t_p))

    absd = np.zeros(m, dtype=np.float64)

    perm_arr = (np.ascontiguousarray(Aord["sperm"], dtype=np.int64).ravel() - 1).astype(np.uintp)
    invperm = np.zeros(m, dtype=np.uintp)
    invperm[perm_arr] = np.arange(m, dtype=np.uintp)

    sdpN = int(cK.sdpN)
    if sdpN > 0:
        blkstart_full = np.ascontiguousarray(K["blkstart"], dtype=np.float64).ravel()
        offset = int(cK.lorN) + 1
        blkstart = (blkstart_full[offset : offset + sdpN + 1] - 1).astype(np.uintp)

        dz = Aord["dz"].tocsc()
        dzstructjc = np.ascontiguousarray(dz.indptr, dtype=np.uintp)
        dzstructir = np.ascontiguousarray(dz.indices, dtype=np.uintp)

        col_widths = np.diff(dzstructjc.astype(np.int64))
        maxadd = int(col_widths.max()) if col_widths.size else 0
        dznnz = int(dzstructjc[m])

        iwsiz = int(np.floor(np.log(1.0 + maxadd) / np.log(2.0))) if maxadd > 0 else 0
        iwsiz += maxadd + 2
        iwsiz = 2 * sdpN + dznnz + max(iwsiz, max(int(cK.rMaxn), int(cK.hMaxn)))
        iwork = np.zeros(max(iwsiz, m), dtype=np.uintp)

        psdNL_arr = np.ascontiguousarray(cK._keepalive[2], dtype=np.float64).astype(np.uintp)

        # xblk[j] = k iff blkstart[k] <= j < blkstart[k+1]; sized to cover
        # the whole 0:lenfull range (not just blkstart[0]:blkstart[-1])
        # rather than replicating the real C code's `xblk -= blkstart[0]`
        # pointer-offset trick -- dzstructir values are always within
        # [blkstart[0], blkstart[sdpN]) (incorder() only ever emits PSD
        # row subscripts there), so the unused low range is simply never
        # read.
        xblk = np.zeros(int(blkstart[sdpN]), dtype=np.uintp)
        for k in range(sdpN):
            xblk[int(blkstart[k]) : int(blkstart[k + 1])] = k

        dzjc = np.zeros(sdpN + 1, dtype=np.uintp)
        _lib.dzblkpartit(
            dzjc.ctypes.data_as(c_size_t_p), dzstructir.ctypes.data_as(c_size_t_p),
            xblk.ctypes.data_as(c_size_t_p), dznnz, sdpN,
        )

        cwork = np.zeros(max(maxadd, 1), dtype=np.uint8)
        fwsiz = lenud + 2 * max(int(cK.rMaxn) ** 2, 2 * int(cK.hMaxn) ** 2)
        fwork = np.zeros(max(fwsiz, 1), dtype=np.float64)

        _lib.getada3(
            ada_struct, absd.ctypes.data_as(c_double_p), at_struct, udsqr_arr.ctypes.data_as(c_double_p),
            Ajc1_arr.ctypes.data_as(c_size_t_p), dzjc.ctypes.data_as(c_size_t_p),
            dzstructjc.ctypes.data_as(c_size_t_p), dzstructir.ctypes.data_as(c_size_t_p),
            blkstart.ctypes.data_as(c_size_t_p), xblk.ctypes.data_as(c_size_t_p),
            psdNL_arr.ctypes.data_as(c_size_t_p),
            perm_arr.ctypes.data_as(c_size_t_p), invperm.ctypes.data_as(c_size_t_p),
            m, lenud, ctypes.byref(cK),
            fwork.ctypes.data_as(c_double_p), fwsiz,
            iwork.ctypes.data_as(c_size_t_p), iwsiz,
            cwork.ctypes.data_as(c_ubyte_p),
        )
    else:
        iwork = np.zeros(max(m, 1), dtype=np.uintp)
        # NOT wrapped via ctypes, deliberately: cpspdiag.c's diagonal
        # lookup goes through blksdp.h's `ibsearch` macro, which casts
        # icmp() -- declared to return `char` -- to a `COMPFUN`
        # (`int(*)(const void*, const void*)`) for bsearch()'s
        # comparator: undefined behavior in C, the same pattern already
        # documented (and worked around the same way, by not binding the
        # C kernel at all) for sortnnz.c/iswnbr.c's qsort comparators --
        # see neighborhood.py's docstring. Confirmed to actually bite
        # here too: on this port's build, bsearch() never finds the
        # diagonal entry at all (every `absd` entry silently comes back
        # 0.0, even though the diagonal is genuinely present and
        # correctly sorted in the input). pattern.diagonal() is exactly
        # cpspdiag's own documented intent (`d := diag(X)`), computed
        # directly by scipy without going through the broken bsearch.
        absd[:] = np.asarray(pattern.diagonal(), dtype=np.float64)

    _lib.spmakesym(ada_struct, m, iwork.ctypes.data_as(c_size_t_p))

    ADA = sp.csc_matrix(
        (ada_pr[:nnz], ada_ir[:nnz].astype(np.int64), ada_jc.astype(np.int64)), shape=(m, m)
    )
    return ADA, absd


_lib.getsplit.argtypes = [c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
                            ctypes.c_size_t, ctypes.c_size_t]
_lib.getsplit.restype = None

_lib.getfirstpiv.argtypes = [c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, ctypes.c_size_t]
_lib.getfirstpiv.restype = None


def mJdetd(detd, K: dict):
    """y = mJdetd(detd, K): y[k-th Lorentz block] = [-detd(k); detd(k)*
    ones(nk-1,1)] i.e. -detd(k)*J with J=diag([1,-I]). mJdetd.c has no
    separable core function (the whole computation lives directly in its
    mexFunction), so this ports that logic straight to NumPy."""
    import numpy as np

    cK = cone_from_dict(K)
    detd = np.ascontiguousarray(detd, dtype=np.float64).ravel()
    y = np.empty(cK.qDim, dtype=np.float64)
    lorNL = cK._keepalive[0]  # "q" array
    i = 0
    for k in range(cK.lorN):
        nk = int(lorNL[k])
        y[i] = -detd[k]
        y[i + 1 : i + nk] = detd[k]
        i += nk
    return y


def sortnnz(At, Ajc1=None, Ajc2=None):
    """perm = sortnnz(At, Ajc1, Ajc2): 1-indexed column permutation of At
    sorting columns by ascending (Ajc2-Ajc1) nnz count (defaults to each
    column's own full nnz range).

    NOT wrapped via ctypes, deliberately: sortnnz.c's kicmp() returns
    `char`, but is called through qsort() via a cast to a `COMPFUN`
    (`int(*)(const void*, const void*)`) function pointer -- undefined
    behavior in C (a qsort comparator must genuinely return int), which
    was confirmed to actually bite here: this port's libsedumi.so build
    and the Octave/MEX build, both compiled from the identical
    sortnnz.c/sdmauxCmp.c, produced DIFFERENT (both internally
    consistent-looking, but different) orderings for tied nnz counts on
    the same input, i.e. this specific C code path is not reliably
    reproducible across builds/compilers at all -- so matching "whatever
    the MEX build happens to do" isn't a well-defined target to bind
    against. This instead implements sortnnz's clearly-stated intent
    (ascending nnz, stable order for ties) directly, which is exactly
    what happened to match the Octave build's output on every fixture
    tried (its qsort() also preserved original order among ties on this
    data) without depending on the same undefined behavior to keep
    doing so.
    """
    import numpy as np

    A = At.tocsc()
    m = A.shape[1]
    Ajc1_arr = A.indptr[:m] if Ajc1 is None else np.ascontiguousarray(Ajc1)
    Ajc2_arr = A.indptr[1:] if Ajc2 is None else np.ascontiguousarray(Ajc2)
    nnz_per_col = np.asarray(Ajc2_arr, dtype=np.int64) - np.asarray(Ajc1_arr, dtype=np.int64)
    order = sorted(range(m), key=lambda k: nnz_per_col[k])
    return (np.array(order, dtype=np.int64) + 1)


def cholsplit(L: dict, cachesize_kb: float = 512):
    """split = cholsplit(L, cachesize_kb): recommends splitting each
    supernode into cache-sized column groups for blkchol's dense
    sub-block updates. Wraps getsplit() (cholsplit.c) directly."""
    import numpy as np

    L_pattern = L["L"].tocsc()
    m = L_pattern.shape[0]
    xsuper = (np.ascontiguousarray(L["xsuper"], dtype=np.int64) - 1).astype(np.uintp)
    nsuper = xsuper.size - 1
    ljc = np.ascontiguousarray(L_pattern.indptr, dtype=np.uintp)
    lir = np.ascontiguousarray(L_pattern.indices, dtype=np.uintp)
    cachesiz = int((0.9 * (1024 / 8)) * cachesize_kb)  # 90% of floats-per-KB
    split = np.zeros(m, dtype=np.uintp)

    _lib.getsplit(
        split.ctypes.data_as(c_size_t_p), ljc.ctypes.data_as(c_size_t_p),
        lir.ctypes.data_as(c_size_t_p), xsuper.ctypes.data_as(c_size_t_p),
        nsuper, cachesiz,
    )
    return split.astype(np.int64)


def finsymbden(LAD, perm, dz, firstq: int):
    """Lden = finsymbden(LAD, perm, dz, firstq): inserts Lorentz-trace
    columns into (perm, dz), producing the "symbolic dense-column"
    structure dpr1fact()/fwdpr1()/bwdpr1() consume (see
    sedumi_port._native.dpr1fact's docstring for the "dz" cumulative-
    compact-row-set convention this also produces). Wraps
    getfirstpiv() (finsymbden.c) directly; the perm/dz remapping itself
    (inserting trace columns) is small enough to port straight to NumPy,
    mirroring finsymbden.c's mexFunction line for line.
    """
    import numpy as np
    import scipy.sparse

    LAD = LAD.tocsc()
    m, n = LAD.shape
    dz = dz.tocsc()
    nperm = dz.shape[1]
    firstQ = firstq - 1
    lastQ = firstQ + n - nperm

    dz_jc = np.ascontiguousarray(dz.indptr, dtype=np.int64)
    dz_ir = np.ascontiguousarray(dz.indices, dtype=np.int64)
    nnzdz = int(dz_jc[nperm])
    invdz = np.zeros(max(m, 1), dtype=np.uintp)
    for i in range(dz_jc[0], nnzdz):
        invdz[dz_ir[i]] = i

    perm_in = np.ascontiguousarray(perm, dtype=np.int64).ravel()
    new_perm = np.zeros(n, dtype=np.uintp)
    dznewJc = np.zeros(n + 1, dtype=np.uintp)
    inz = 0
    for i in range(nperm):
        j = int(perm_in[i]) - 1
        new_perm[inz] = j
        dznewJc[inz] = dz_jc[i]
        inz += 1
        if firstQ <= j < lastQ:
            new_perm[inz] = nperm + j - firstQ
            dznewJc[inz] = dz_jc[i + 1]
            inz += 1
    assert inz == n
    dznewJc[n] = dz_jc[nperm]

    firstpiv = np.zeros(max(n, 1), dtype=np.uintp)
    LADjc = np.ascontiguousarray(LAD.indptr, dtype=np.uintp)
    LADir = np.ascontiguousarray(LAD.indices, dtype=np.uintp)
    _lib.getfirstpiv(
        firstpiv.ctypes.data_as(c_size_t_p), invdz.ctypes.data_as(c_size_t_p),
        dznewJc.ctypes.data_as(c_size_t_p), LADjc.ctypes.data_as(c_size_t_p),
        LADir.ctypes.data_as(c_size_t_p), n,
    )

    dz_new_ir = dz_ir[: int(dznewJc[n])].astype(np.int64)
    dz_new = scipy.sparse.csc_matrix(
        (np.ones(max(len(dz_new_ir), 1), dtype=np.float64)[: len(dz_new_ir)],
         dz_new_ir, dznewJc.astype(np.int64)),
        shape=(m, n),
    )

    return {
        "LAD": LAD,
        "perm": (new_perm + 1).astype(np.int64),
        "dz": dz_new,
        "first": (firstpiv + 1).astype(np.int64),
    }


def cone_from_dict(K: dict) -> ConeK:
    """Build a ConeK from a plain dict shaped like SeDuMi's K struct, e.g.
    {"f": 2, "l": 3, "q": [4], "s": [2, 3]}. Mirrors what conepars() does
    for the MATLAB/Octave K struct, with no mxArray/MATLAB/Octave layer
    anywhere in the path.
    """
    import numpy as np

    raw = SedumiKRaw()
    raw.f = float(K.get("f", 0.0))
    raw.l = float(K.get("l", 0.0))

    # Keep the backing numpy arrays alive for the duration of the call by
    # returning them alongside; conepars_raw() only reads them, so this
    # scope is enough (it doesn't retain the pointers past the call).
    q = np.ascontiguousarray(K.get("q", []), dtype=np.float64)
    r = np.ascontiguousarray(K.get("r", []), dtype=np.float64)
    s = np.ascontiguousarray(K.get("s", []), dtype=np.float64)
    raw.q = q.ctypes.data_as(c_double_p) if q.size else None
    raw.qN = q.size
    raw.r = r.ctypes.data_as(c_double_p) if r.size else None
    raw.rN = r.size
    raw.s = s.ctypes.data_as(c_double_p) if s.size else None
    raw.sN = s.size

    if "rsdpN" in K:
        raw.rsdpNgiven = b"\x01"
        raw.rsdpN = float(K["rsdpN"])
    else:
        raw.rsdpNgiven = b"\x00"

    raw.statsGiven = b"\x00"  # always recompute; the precomputed-stats
    # path exists to mirror the MEX adapter, not needed from Python

    cone = ConeK()
    _lib.conepars_raw(ctypes.byref(raw), ctypes.byref(cone))
    # cone.{lorNL,rconeNL,sdpNL} are just copies of raw.{q,r,s}, i.e.
    # pointers into q/r/s above -- keep them alive as long as `cone` is,
    # or they'd dangle the moment this function returns and q/r/s get
    # garbage collected.
    cone._keepalive = (q, r, s)
    return cone


# symbfwmat() grows its output row-index buffer with mxRealloc() (=
# realloc(), see sedumi_platform.h) as it fills it in, then shrinks it
# down to the exact final size at the end -- it must therefore own a
# buffer allocated via the C allocator family, NOT a numpy-owned buffer
# (realloc()'ing memory numpy itself allocated would corrupt numpy's own
# bookkeeping), and freed with a calloc/free that shares the same heap
# as whatever realloc() libsedumi's C code itself resolves to, or
# free()'ing memory realloc() may have moved would corrupt the CRT's own
# bookkeeping. On macOS/Linux, CDLL(None) resolves symbols already
# linked into the process, which always includes the C runtime -- there
# is no separate libc to load. Windows has no such "the process itself"
# handle (ctypes.CDLL(None) raises TypeError there), and no single
# universal CRT either; libsedumi.dll here is always built by MSYS2's
# MINGW64 gcc (see tools/build_libsedumi.sh), which -- unlike UCRT64/
# CLANG64 -- links the legacy msvcrt.dll, so loading that by name gets
# the exact same heap its realloc() uses.
_libc = ctypes.CDLL("msvcrt") if sys.platform == "win32" else ctypes.CDLL(None)
_libc.calloc.argtypes = [ctypes.c_size_t, ctypes.c_size_t]
_libc.calloc.restype = ctypes.c_void_p
_libc.free.argtypes = [ctypes.c_void_p]
_libc.free.restype = None

c_size_t_pp = ctypes.POINTER(c_size_t_p)

_lib.snodeCompress.argtypes = [
    c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
    ctypes.c_size_t,
]
_lib.snodeCompress.restype = None

_lib.symbfwmat.argtypes = [
    c_size_t_p, c_size_t_pp, c_size_t_p, c_size_t_p, c_size_t_p,
    c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p, c_size_t_p,
    ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
    c_size_t_p, c_ubyte_p,
]
_lib.symbfwmat.restype = None


def symbfwblk(L: dict, B):
    """X = symbfwblk(L, B): symbolic sparsity pattern of `L\\B(L["perm"],:)`
    -- i.e. which entries of the forward solve with the sparse Cholesky
    factor `L["L"]` (restricted to `B`'s columns, permuted by `L["perm"]`)
    could be nonzero, WITHOUT doing any arithmetic (matching
    `symbfwmat()`'s pure symbolic nature -- every returned entry's value
    is 1.0). `L` needs "perm" (1-indexed, length m), "L" (sparse m x m
    Cholesky-factor pattern), "xsuper" (1-indexed supernode boundaries,
    length nsuper+1) -- exactly the `L` struct `symbchol()`/`sfinit`
    already produce. `B` is a sparse m x n matrix (only its pattern is
    read). Wraps symbfwmat() (symbfwblk.c) plus its own snodeCompress()
    helper, exactly as symbcholden.m's mex wrapper does.
    """
    import numpy as np
    import scipy.sparse as sp

    L_pattern = L["L"].tocsc()
    m = L_pattern.shape[0]
    ljc = np.ascontiguousarray(L_pattern.indptr, dtype=np.uintp)
    lir = np.ascontiguousarray(L_pattern.indices, dtype=np.uintp)

    xsuper0 = (np.ascontiguousarray(L["xsuper"], dtype=np.int64).ravel() - 1).astype(np.uintp)
    nsuper = xsuper0.size - 1

    perm0 = (np.ascontiguousarray(L["perm"], dtype=np.int64).ravel() - 1).astype(np.uintp)
    invperm = np.zeros(m, dtype=np.uintp)
    invperm[perm0] = np.arange(m, dtype=np.uintp)

    B_csc = B.tocsc() if sp.issparse(B) else sp.csc_matrix(B)
    n = B_csc.shape[1]
    bjc = np.ascontiguousarray(B_csc.indptr, dtype=np.uintp)
    bir = np.ascontiguousarray(B_csc.indices, dtype=np.uintp)

    nnzL = int(L_pattern.indptr[m])
    xlindx = np.zeros(nsuper + 1, dtype=np.uintp)
    lindx = np.zeros(max(nnzL, 1), dtype=np.uintp)
    snode = np.zeros(m, dtype=np.uintp)

    _lib.snodeCompress(
        xlindx.ctypes.data_as(c_size_t_p), lindx.ctypes.data_as(c_size_t_p),
        snode.ctypes.data_as(c_size_t_p), ljc.ctypes.data_as(c_size_t_p),
        lir.ctypes.data_as(c_size_t_p), xsuper0.ctypes.data_as(c_size_t_p), nsuper,
    )

    xjc = np.zeros(n + 1, dtype=np.uintp)
    snodefrom = np.zeros(max(nsuper, 1), dtype=np.uintp)
    processed = np.zeros(max(nsuper, 1), dtype=np.uint8)

    # Every column of x=L\b(perm,:) has at most m nonzeros (rows partition
    # disjointly across supernodes), so m*n is a hard upper bound on the
    # total nnz -- allocating it up front guarantees symbfwmat() never
    # needs to grow this buffer, only shrink it at the very end.
    maxnnz0 = max(m * n, 1)
    raw = _libc.calloc(maxnnz0, ctypes.sizeof(ctypes.c_size_t))
    if not raw:
        raise MemoryError("symbfwblk: calloc failed")
    xir_ptr = ctypes.cast(raw, c_size_t_p)
    maxnnz_c = ctypes.c_size_t(maxnnz0)

    try:
        _lib.symbfwmat(
            xjc.ctypes.data_as(c_size_t_p), ctypes.byref(xir_ptr), ctypes.byref(maxnnz_c),
            bjc.ctypes.data_as(c_size_t_p), bir.ctypes.data_as(c_size_t_p),
            invperm.ctypes.data_as(c_size_t_p), snode.ctypes.data_as(c_size_t_p),
            xsuper0.ctypes.data_as(c_size_t_p), xlindx.ctypes.data_as(c_size_t_p),
            lindx.ctypes.data_as(c_size_t_p), nsuper, m, n,
            snodefrom.ctypes.data_as(c_size_t_p), processed.ctypes.data_as(c_ubyte_p),
        )
        final_nnz = int(xjc[n])
        if final_nnz > 0:
            xir_out = np.ctypeslib.as_array(xir_ptr, shape=(final_nnz,)).astype(np.int64).copy()
        else:
            xir_out = np.zeros(0, dtype=np.int64)
    finally:
        _libc.free(ctypes.cast(xir_ptr, ctypes.c_void_p))

    data = np.ones(final_nnz, dtype=np.float64)
    return sp.csc_matrix((data, xir_out, xjc.astype(np.int64)), shape=(m, n))

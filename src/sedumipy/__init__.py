"""sedumipy: a MATLAB/Octave-free port of SeDuMi (SDP/SOCP interior-point
solver), in pure Python (NumPy/SciPy) plus a small standalone C kernel
library (libsedumi.so, via ctypes -- see _native.py).

    import sedumipy
    x, y, info = sedumipy.sedumi(A, b, c, K)

is equivalent to real SeDuMi's own `[x,y,info] = sedumi(A,b,c,K)`. See
sedumi.py's own docstring for the exact scope (LP + second-order-cone +
PSD cones) and CONTRIBUTING.md for the full project status.
"""

from . import _native  # noqa: F401
from .matio import read_mat, write_solution_mat  # noqa: F401
from .sedumi import sedumi  # noqa: F401

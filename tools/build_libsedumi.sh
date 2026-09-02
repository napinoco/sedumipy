#!/usr/bin/env bash
# build_libsedumi.sh
#
# Phase 1: builds every SeDuMi C kernel into a single standalone shared
# library with NO MATLAB/Octave/MEX dependency at all (no mex.h, no
# mxArray, nothing). This is the artifact Phase 2 (Python bindings) will
# link against.
#
# bwblkslv2.c is intentionally excluded: it is dead code (not referenced
# by install_sedumi.m or any .m file) that already failed to compile
# before this port -- it #includes a "blkchol.h" that does not exist in
# this repository. Left untouched rather than guessed at.
#
# Usage (from the repository root):
#   tools/build_libsedumi.sh [output_path]

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
csrc_dir="$repo_root/csrc"
out="${1:-/tmp/libsedumi.so}"

cd "$csrc_dir"

sources=()
for f in *.c; do
  if [ "$f" = "bwblkslv2.c" ]; then
    continue
  fi
  sources+=("$f")
done

echo "Building $out from ${#sources[@]} source files (SEDUMI_STANDALONE, no mex.h)..."

# symfct.c (a forked f2c-style translation) declares several functions
# the old K&R way: an `extern int foo();` forward declaration with empty
# parens, meaning "unspecified arguments" in C17 and earlier, followed by
# a real definition that takes arguments. C23 redefines empty parens to
# mean "zero arguments" (matching "(void)"), which turns that into a
# hard prototype mismatch -- and GCC 14+ defaults to -std=gnu23. This
# hit MSYS2's gcc 16 on Windows first (MSYS2 packages track upstream gcc
# closely), but would eventually hit any platform's gcc/clang once its
# default catches up. -std=gnu17 pins the pre-C23 semantics this code
# was written against, rather than rewriting decades-old, numerically
# load-bearing C for a language-standard technicality.
std_flag=-std=gnu17

kernel="$(uname -s)"
case "$kernel" in
  Darwin)
    # macOS ships no `libblas`/`-lblas` by default and has no system
    # package manager to install one; every Mac does ship
    # Accelerate.framework (BLAS/LAPACK), so link against that instead of
    # requiring Homebrew. `cc` (not `gcc`, which on macOS is usually just
    # a clang alias anyway) for the same reason -- no assumption of a
    # real GNU toolchain.
    cc "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" \
      -framework Accelerate -lm
    ;;
  MINGW*|MSYS*)
    # Run from an MSYS2 MinGW64 shell, with mingw-w64-x86_64-gcc and
    # mingw-w64-x86_64-openblas installed (`pacman -S`) -- see
    # CONTRIBUTING.md's Windows note. This produces a plain PE DLL
    # ctypes.CDLL() loads like any other shared library (libsedumi.dll is
    # not a PyInit_*-exporting CPython extension, so it doesn't need to
    # be built with the same compiler as Python itself). `-lopenblas`
    # (not `-lblas`): MSYS2's BLAS/LAPACK meta-packages don't provide a
    # plain `libblas`, but OpenBLAS exports the same FORT()-mangled
    # symbol names (see sedumi_platform.h), so it's a drop-in swap.
    # `libopenblas.dll` itself (and the mingw runtime DLLs) are not
    # statically linked here -- see setup.py's Windows note on
    # `delvewheel` bundling them into the wheel instead.
    #
    # Deliberately NOT a bare `gcc`/PATH lookup: GitHub Actions' Windows
    # runners ship their own unrelated MinGW toolchain preinstalled at
    # C:\mingw64 (confirmed: cibuildwheel's build picked up its gcc 15
    # instead of MSYS2's pacman-installed gcc 16, silently, then failed
    # to link -lopenblas since only the MSYS2 install has it) -- the same
    # PATH-shadowing class of problem as bash's own resolution (see
    # setup.py/_native.py's notes), just a second, unrelated toolchain
    # this time instead of WSL's bash.exe stub. `/mingw64` (a plain,
    # constant POSIX path, no computation needed) is always correct
    # inside *any* MSYS2 bash regardless of install location or how bash
    # itself was invoked: MSYS2's runtime resolves `/` relative to
    # msys-2.0.dll's own location, and mingw64/ is always its sibling.
    # (An earlier version of this tried to derive the same path from
    # $BASH's own dirname instead and undercounted a `..` level, landing
    # on /usr/mingw64 -- unnecessary and wrong; this is simpler and
    # actually correct.)
    /mingw64/bin/gcc.exe "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -shared -I. "${sources[@]}" -o "$out" \
      -L/mingw64/lib -lopenblas -lm
    ;;
  *)
    gcc "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" -lblas -lm
    ;;
esac

echo "OK: $out"
case "$kernel" in
  Darwin)
    nm -gU "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions"
    ;;
  MINGW*|MSYS*)
    # PE/COFF has no ELF-style .dynsym for `nm -D` to read; this is an
    # approximate count of defined text symbols instead, diagnostic only.
    /mingw64/bin/nm.exe "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions (approx.)"
    ;;
  *)
    nm -D "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions"
    ;;
esac

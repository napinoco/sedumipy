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
kernel="$(uname -s)"
case "$kernel" in
  Darwin)
    # macOS ships no `libblas`/`-lblas` by default and has no system
    # package manager to install one; every Mac does ship
    # Accelerate.framework (BLAS/LAPACK), so link against that instead of
    # requiring Homebrew. `cc` (not `gcc`, which on macOS is usually just
    # a clang alias anyway) for the same reason -- no assumption of a
    # real GNU toolchain.
    cc -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" \
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
    gcc -DSEDUMI_STANDALONE -O2 -Wall -shared -I. "${sources[@]}" -o "$out" \
      -lopenblas -lm
    ;;
  *)
    gcc -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" -lblas -lm
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
    nm "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions (approx.)"
    ;;
  *)
    nm -D "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions"
    ;;
esac

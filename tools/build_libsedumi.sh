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
    # this time instead of WSL's bash.exe stub.
    #
    # Also deliberately NOT MSYS2's own /mingw64 POSIX mount, and NOT
    # just a single hardcoded Windows-style path either: two prior
    # versions of this each seemed reasonable and each still failed with
    # "No such file or directory" in real CI, despite mingw-w64-x86_64-
    # gcc having just been installed successfully by the very same job.
    # Rather than guess a fourth time, try every plausible form and, if
    # every single one fails, dump enough diagnostics to actually see
    # what's on disk instead of failing silently.
    msys_root="${MSYS2_ROOT:-C:/msys64}"
    gcc_candidates=(
      "$msys_root/mingw64/bin/gcc.exe"
      "/mingw64/bin/gcc.exe"
      "$(cygpath -u "$msys_root" 2>/dev/null || true)/mingw64/bin/gcc.exe"
    )
    gcc_bin=""
    for candidate in "${gcc_candidates[@]}"; do
      if [ -n "$candidate" ] && [ -x "$candidate" ]; then
        gcc_bin="$candidate"
        break
      fi
    done
    if [ -z "$gcc_bin" ]; then
      gcc_bin="$(command -v gcc 2>/dev/null || true)"
    fi
    if [ -z "$gcc_bin" ]; then
      echo "ERROR: could not locate MSYS2's mingw-w64-x86_64-gcc under any" >&2
      echo "of these forms: ${gcc_candidates[*]}, nor via PATH. Diagnostics:" >&2
      echo "  MSYS2_ROOT=${MSYS2_ROOT:-<unset>}  msys_root=$msys_root  PATH=$PATH" >&2
      echo "  ls -la \"$msys_root/mingw64/bin\":" >&2
      ls -la "$msys_root/mingw64/bin" >&2 2>&1 || echo "  (that directory listing itself failed)" >&2
      exit 1
    fi
    echo "Using gcc: $gcc_bin"
    gcc_bin_dir="$(dirname "$gcc_bin")"
    "$gcc_bin" "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -shared -I. "${sources[@]}" -o "$out" \
      -L"$gcc_bin_dir/../lib" -lopenblas -lm
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
    # `$gcc_bin_dir` was resolved above, next to the gcc that built $out.
    "$gcc_bin_dir/nm.exe" "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions (approx.)"
    ;;
  *)
    nm -D "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions"
    ;;
esac

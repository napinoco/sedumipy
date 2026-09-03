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
#   tools/build_libsedumi.sh [output_path] [python_executable]
#
# BLAS: on Linux and Windows, prefers scipy-openblas64 (see below) when
# the given Python can import it, falling back to that OS's own
# system-BLAS story (MSYS2's -lopenblas on Windows, -lblas/-lopenblas on
# Linux) only when it can't -- that fallback is what keeps a plain `pip
# install -e .` working offline, or on a platform scipy-openblas64
# doesn't ship for. macOS always links the system Accelerate framework
# instead, regardless of scipy-openblas64's availability -- see the
# Darwin case below for why.

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
csrc_dir="$repo_root/csrc"
out="${1:-/tmp/libsedumi.so}"
python_bin="${2:-python3}"

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

# scipy-openblas64: a pip-installable, prebuilt ILP64 OpenBLAS -- the
# same package numpy/scipy themselves build against -- with wheels for
# Linux (x86_64/aarch64/...) and Windows (amd64/arm64) (also macOS, but
# see the Darwin case below for why this build doesn't use it there).
# Using it means the BLAS this library links on Linux/Windows no longer
# has to be separately installed or built per OS (no
# libblas-dev/libopenblas-dev, no MSYS2 -lopenblas): `pip install
# scipy-openblas64` covers both with one command. It's a
# *build-time-only* dependency (never installed at runtime for end
# users): a wheel build vendors the resulting shared library into the
# wheel itself (auditwheel/delvewheel, same as this script already
# relied on for -lopenblas), exactly like numpy/scipy's own wheels do.
#
# It ships as an ILP64 build (64-bit Fortran integers, not the 32-bit
# LP64 int reference BLAS/OpenBLAS/Accelerate use -- see
# sedumi_platform.h's SEDUMI_BLAS_ILP64), and its exported symbols carry
# an extra prefix/suffix (e.g. "dscal_" -> "scipy_dscal_64_") so several
# independently-built OpenBLAS copies can coexist in one process without
# clashing. Both are handled by -DSEDUMI_BLAS_ILP64 and the
# -DBLAS_SYMBOL_PREFIX/-DBLAS_SYMBOL_SUFFIX flags below, which come
# straight from the package's own pkg-config metadata rather than being
# hardcoded here -- so a version bump changing either can't silently
# desync this build from what the library actually exports.
scipy_openblas64_flags() {
  "$python_bin" - <<'PYEOF'
import re
import sys

try:
    import scipy_openblas64 as sob
except ImportError:
    sys.exit(1)

includedir = sob.get_include_dir()
libdir = sob.get_lib_dir()
pkg_config = sob.get_pkg_config()


def field(name):
    m = re.search(rf"^{name}:[ \t]*(.*)$", pkg_config, re.MULTILINE)
    if not m:
        sys.exit(f"scipy_openblas64.get_pkg_config() has no {name}: line")
    return m.group(1).strip()


cflags = field("Cflags").replace("${includedir}", includedir)
libs = field("Libs").replace("${libdir}", libdir)
print(cflags)
print(libs)
PYEOF
}

kernel="$(uname -s)"

# Not even attempted on Darwin: see the Accelerate note in the Darwin
# case below for why macOS never prefers scipy-openblas64, regardless of
# whether python_bin can import it -- skip the subprocess (and the
# "Using scipy-openblas64" echo, which would otherwise be actively
# misleading there since that branch ignores these variables entirely).
sob_cflags=""
sob_libs=""
if [ "$kernel" != "Darwin" ] && sob_flags="$(scipy_openblas64_flags 2>/dev/null)"; then
  sob_cflags="$(sed -n '1p' <<<"$sob_flags")"
  sob_libs="$(sed -n '2p' <<<"$sob_flags")"
  echo "Using scipy-openblas64 (ILP64) via $python_bin for BLAS"
fi

case "$kernel" in
  Darwin)
    # Deliberately NOT preferring scipy-openblas64 here even when it's
    # importable, unlike every other branch below: scipy-openblas64
    # exists to avoid the *build/install* pain Linux (system package
    # divergence) and Windows (the whole MSYS2 story) have, and macOS
    # never had that pain -- every Mac already ships Accelerate.framework
    # (BLAS/LAPACK) with zero install step, so there's nothing for a pip
    # package to save here. It would only add wheel bytes (a vendored
    # scipy-openblas64 plus its own libgfortran/libquadmath) and risk
    # being slower than Accelerate's own Apple Silicon tuning, for no
    # offsetting benefit. `cc` (not `gcc`, which on macOS is usually just
    # a clang alias anyway) for the same reason -- no assumption of a
    # real GNU toolchain, no Homebrew dependency.
    cc "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" \
      -framework Accelerate -lm
    ;;
  MINGW*|MSYS*)
    # Still needs an MSYS2 MinGW64 shell with mingw-w64-x86_64-gcc on
    # PATH (`pacman -S` -- see CONTRIBUTING.md's Windows note) even when
    # scipy-openblas64 supplies the BLAS below: MSVC can't build this
    # (see the WINDOWS_BUILD_HELP text in setup.py), so a real gcc is
    # needed regardless of where the BLAS comes from. This produces a
    # plain PE DLL ctypes.CDLL() loads like any other shared library
    # (libsedumi.dll is not a PyInit_*-exporting CPython extension, so it
    # doesn't need to be built with the same compiler as Python itself).
    # Neither the BLAS DLL nor the mingw runtime DLLs are statically
    # linked here -- see setup.py's Windows note on `delvewheel` bundling
    # them into the wheel instead.
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
      # Every candidate location missed -- print what's actually there
      # unconditionally (not just on total failure): a PATH-based `gcc`
      # match below can silently be the *wrong* gcc (GitHub Actions'
      # Windows runners ship an unrelated MinGW toolchain preinstalled
      # at C:\mingw64, with no openblas -- see the comment above), so
      # this is the only chance to see the real state of disk.
      echo "WARNING: none of these MSYS2 gcc candidates matched:" >&2
      echo "  ${gcc_candidates[*]}" >&2
      echo "  MSYS2_ROOT=${MSYS2_ROOT:-<unset>}  msys_root=$msys_root" >&2
      echo "  ls -la \"$msys_root\":" >&2
      ls -la "$msys_root" >&2 2>&1 || echo "  (that listing itself failed -- msys_root doesn't exist)" >&2
      echo "  ls -la \"$msys_root/mingw64/bin\" (if present):" >&2
      ls -la "$msys_root/mingw64/bin" >&2 2>&1 || echo "  (that listing itself failed)" >&2
      echo "  PATH=$PATH" >&2
      # Fall back to PATH, but reject a match under the known-wrong
      # C:\mingw64 (case-insensitively, since bash paths here can come
      # back as either C:\mingw64 or /c/mingw64) rather than silently
      # building with a toolchain that has no -lopenblas.
      path_gcc="$(command -v gcc 2>/dev/null || true)"
      case "${path_gcc,,}" in
        *msys64*) gcc_bin="$path_gcc" ;;
        "") ;;
        *) echo "  PATH's gcc ($path_gcc) looks like the wrong (non-MSYS2) toolchain -- rejecting it." >&2 ;;
      esac
    fi
    if [ -z "$gcc_bin" ]; then
      echo "ERROR: could not locate MSYS2's mingw-w64-x86_64-gcc anywhere usable (see diagnostics above)." >&2
      exit 1
    fi
    echo "Using gcc: $gcc_bin"
    gcc_bin_dir="$(dirname "$gcc_bin")"
    if [ -n "$sob_libs" ]; then
      "$gcc_bin" "$std_flag" -DSEDUMI_STANDALONE -DSEDUMI_BLAS_ILP64 -O2 -Wall -shared -I. \
        $sob_cflags "${sources[@]}" -o "$out" $sob_libs -lm
    else
      # `-lopenblas` (not `-lblas`): MSYS2's BLAS/LAPACK meta-packages
      # don't provide a plain `libblas`, but OpenBLAS exports the same
      # FORT()-mangled symbol names (see sedumi_platform.h), so it's a
      # drop-in swap. Requires mingw-w64-x86_64-openblas installed
      # alongside mingw-w64-x86_64-gcc.
      "$gcc_bin" "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -shared -I. "${sources[@]}" -o "$out" \
        -L"$gcc_bin_dir/../lib" -lopenblas -lm
    fi
    ;;
  *)
    if [ -n "$sob_libs" ]; then
      gcc "$std_flag" -DSEDUMI_STANDALONE -DSEDUMI_BLAS_ILP64 -O2 -Wall -fPIC -shared -I. \
        $sob_cflags "${sources[@]}" -o "$out" $sob_libs -lm
    else
      # Prefer OpenBLAS over the reference Netlib implementation when
      # both are available. They export the same FORT()-mangled symbols
      # (see sedumi_platform.h), so it is a drop-in swap, but not a free
      # one to skip: measured on one box, OpenBLAS runs the Level-1
      # kernels this library actually calls about 3x faster (ddot 0.40s
      # -> 0.12s, daxpy 0.27s -> 0.12s), and all of them are on the hot
      # path -- realdot alone is called from 15 files, including the
      # triangular solves PCG repeats every iteration.
      #
      # Settled by an actual link test rather than by looking for files:
      # distributions disagree about whether -lblas already resolves to
      # OpenBLAS (Debian/Ubuntu route it through update-alternatives, so
      # it often does; RHEL/AlmaLinux keep them separate, so it does
      # not), and the linker is the authority on what it can actually
      # find.
      blas_lib=-lblas
      if echo 'int main(void){return 0;}' | gcc -x c - -lopenblas -o /dev/null 2>/dev/null; then
        blas_lib=-lopenblas
      fi
      echo "Linking BLAS via $blas_lib"
      gcc "$std_flag" -DSEDUMI_STANDALONE -O2 -Wall -fPIC -shared -I. "${sources[@]}" -o "$out" "$blas_lib" -lm
    fi
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
    # $out's actual DLL import table -- libopenblas.dll and the mingw
    # runtime DLLs should appear here, since -lopenblas links against
    # MSYS2's libopenblas.dll.a *import* library (the package ships no
    # static archive). That is what pyproject.toml's delvewheel
    # `--analyze-existing` has to vendor in for a redistributable wheel;
    # without that flag delvewheel silently bundled nothing and claimed
    # "no external dependencies are needed", which this makes checkable
    # rather than assumed. Note `pip` swallows this script's output
    # unless the build fails, so it's visible on a direct invocation
    # (see CONTRIBUTING.md's Windows note), not in a pip install log.
    "$gcc_bin_dir/objdump.exe" -p "$out" 2>/dev/null | grep -i "DLL Name" || true
    ;;
  *)
    nm -D "$out" 2>/dev/null | grep -c ' T ' | xargs -I{} echo "  {} exported functions"
    ;;
esac

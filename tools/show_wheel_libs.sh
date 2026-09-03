#!/usr/bin/env bash
# show_wheel_libs.sh
#
# Prints the shared libraries each given wheel actually carries, and
# fails if a platform that is supposed to vendor its BLAS did not.
#
# This exists because "the repair step ran" is not the same as "the
# repair step bundled anything", and the difference is invisible unless
# you look. Both wheel-repair tools have already been caught claiming
# success while shipping nothing usable:
#
#   - delvewheel logged "no external dependencies are needed" on every
#     Windows build, because it only analyzes .pyd extension modules by
#     default and this wheel has none (libsedumi.dll is a plain
#     ctypes-loaded DLL). Fixed with --analyze-existing; see
#     pyproject.toml.
#   - auditwheel does vendor libblas/libopenblas correctly on Linux, but
#     prints no "Grafting:" line while doing it -- whose absence was
#     briefly and wrongly read as evidence that it had bundled nothing.
#
# So: check the wheel, not the log.
#
# Usage: tools/show_wheel_libs.sh wheelhouse/*.whl

set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "usage: $0 <wheel> [wheel...]" >&2
  exit 2
fi

# macOS is the one platform that legitimately bundles nothing: the
# Darwin build links Accelerate, a system framework present on every
# Mac, so there is no BLAS to carry along.
expect_bundled=1
case "$(uname -s)" in
  Darwin) expect_bundled=0 ;;
esac

# The wheel jobs don't run actions/setup-python (cibuildwheel provisions
# its own interpreters), so don't assume which name the runner's system
# Python answers to: ubuntu-latest and macos-latest may only have
# `python3`, while Windows only has `python`.
py=""
for candidate in python3 python; do
  if command -v "$candidate" >/dev/null 2>&1; then
    py="$candidate"
    break
  fi
done
if [ -z "$py" ]; then
  echo "ERROR: no python3/python on PATH to inspect the wheels with." >&2
  exit 2
fi

status=0
for wheel in "$@"; do
  echo "== $(basename "$wheel")"
  bundled="$("$py" -c '
import sys, zipfile
names = zipfile.ZipFile(sys.argv[1]).namelist()
for n in sorted(names):
    low = n.lower()
    if ".libs/" in low or ".dylibs/" in low or low.endswith((".so", ".dll", ".dylib")):
        if not n.endswith("/"):
            print(n)
' "$wheel")"

  if [ -z "$bundled" ]; then
    echo "  (no shared libraries found at all -- the wheel has no built library?)"
    status=1
    continue
  fi
  echo "$bundled" | sed 's/^/  /'

  if [ "$expect_bundled" -eq 1 ]; then
    # Anything under a .libs/ or .dylibs/ directory was put there by the
    # repair step; the library the build produced itself is not.
    if ! echo "$bundled" | grep -qiE '\.(libs|dylibs)/'; then
      echo "  ERROR: nothing was vendored by the repair step on this platform." >&2
      echo "         The wheel would need a BLAS already installed on the" >&2
      echo "         user's machine -- see this script's header." >&2
      status=1
    fi
  fi
done

exit "$status"

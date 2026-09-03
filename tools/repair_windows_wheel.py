"""Repairs a Windows wheel built by cibuildwheel, via `delvewheel repair`.

libsedumi.dll (see tools/build_libsedumi.sh) dynamically depends on
whichever BLAS backed the build -- normally libscipy_openblas64_.dll,
sitting inside the scipy-openblas64 package's own site-packages
directory, not somewhere delvewheel would find on its own -- plus a
couple of mingw runtime DLLs from the MSYS2 gcc toolchain used to
compile it (see build_libsedumi.sh's Windows note on why MSVC can't).
None of those are on a plain end user's machine, so `delvewheel repair`
has to bundle them into the wheel, the same role auditwheel/delocate
play on Linux/macOS -- just not something cibuildwheel wires up for
Windows by default.

A plain `delvewheel repair --add-path <hardcoded dirs>` command in
pyproject.toml's repair-wheel-command can't express this: the
scipy-openblas64 install directory is a site-packages path that only
exists inside the venv cibuildwheel just built in, so it has to be
resolved at repair time, in Python, rather than hardcoded in TOML.
"""

import subprocess
import sys


def main() -> None:
    wheel, dest_dir = sys.argv[1], sys.argv[2]

    add_paths = [
        # See .github/workflows/wheels.yml's "Verify MSYS2 gcc/openblas,
        # then copy it out of C:\\msys64" step -- the build itself
        # resolves gcc from this copy via MSYS2_ROOT, so the mingw
        # runtime DLLs delvewheel needs to bundle live here too, not in
        # the original MSYS2 install.
        r"C:\toolchain-mingw64\mingw64\bin",
        r"C:\msys64\usr\bin",
    ]

    try:
        import scipy_openblas64 as sob
    except ImportError:
        # Falls back to build_libsedumi.sh's own -lopenblas/MSYS2 path
        # (see its Windows branch) when scipy-openblas64 isn't
        # installed -- nothing extra to add here in that case, since
        # libopenblas.dll already lives under the mingw64 dir above.
        pass
    else:
        add_paths.append(sob.get_lib_dir())

    # --analyze-existing is required here, not optional: delvewheel only
    # analyzes a wheel's *extension modules* (.pyd) by default, and this
    # wheel has none -- libsedumi.dll is a plain ctypes-loaded DLL (see
    # setup.py's docstring), so without this flag delvewheel finds
    # nothing to analyze and reports "no external dependencies are
    # needed", bundling nothing at all. That is not a sign the library
    # is self-contained -- it fails to load on any machine that doesn't
    # happen to have the same DLLs already on PATH.
    command = [
        sys.executable,
        "-m",
        "delvewheel",
        "repair",
        "--analyze-existing",
        "--add-path",
        ";".join(add_paths),
        "-w",
        dest_dir,
        wheel,
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()

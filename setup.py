"""Phase 6: builds libsedumi.so/.dylib as part of `pip install`/`pip
wheel` (via a custom build_ext step that shells out to
tools/build_libsedumi.sh, the same compile command Phase 1 established),
so a built wheel carries a ready-to-load binary instead of relying on
_native.py's on-first-import fallback build (still there, and still what
an editable/development install uses -- see _native.py's own docstring).

An Extension with no sources is registered purely so setuptools marks the
wheel as platform-specific (has_ext_modules() -> True) rather than a
"py3-none-any" pure-Python tag, which would be wrong for a wheel that
bundles a compiled, platform-dependent shared library. The actual
compilation happens in BuildLibsedumi.build_extension below, not through
setuptools' normal C-extension compiler invocation (libsedumi.so is a
plain ctypes-loaded shared library, not a `PyInit_*`-exporting Python
extension module, so it doesn't go through the usual Extension build
path or get a platform/ABI-tagged filename).

Windows: tools/build_libsedumi.sh is a bash script, so it's invoked via
`bash`'s full path explicitly rather than executed directly (Windows has
no shebang support, and passing the bare name "bash" is not enough --
see build_bash_command()'s own docstring for why). This assumes an
MSYS2 MinGW64 environment on PATH (`bash`, `gcc`, and
`mingw-w64-x86_64-openblas` -- see CONTRIBUTING.md's Windows note); the
resulting libsedumi.dll dynamically depends on libopenblas.dll and a
couple of mingw runtime DLLs that are not on a plain end-user's Windows
install, so a wheel built here is not yet redistributable as-is --
.github/workflows/wheels.yml's Windows job runs `delvewheel repair` as a
post-build step to bundle them into the wheel.
"""

import shutil
import subprocess
import sys
from pathlib import Path

# Reaching this means pip is compiling from the source distribution --
# almost always because no prebuilt wheel matched the platform. Say so:
# the bare "no bash on PATH" this replaced left an end user, who never
# asked to compile anything, to work out both why a Python install
# wanted a Unix shell and what to do about it. Keep the wheel coverage
# list in step with pyproject.toml's [tool.cibuildwheel] build/skip.
WINDOWS_BUILD_HELP = """\
sedumipy is being built from source, which on Windows needs an MSYS2
MinGW64 toolchain -- and no `bash` was found on PATH.

You are most likely seeing this because pip fell back to the source
distribution: prebuilt wheels are published for 64-bit x86 Windows,
Linux (x86_64) and macOS (Apple silicon) on CPython 3.10-3.13, so any
other target -- 32-bit or ARM Windows, Alpine/musl, Linux aarch64,
Intel macOS -- is compiled here instead.

To build it, install MSYS2 from https://www.msys2.org/ and run this in
an MSYS2 MinGW64 shell:

    pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-openblas

then add both C:\\msys64\\usr\\bin and C:\\msys64\\mingw64\\bin to PATH
and reinstall. MSVC will not work in its place: libsedumi.dll is a
plain ctypes-loaded DLL rather than a CPython extension, and the build
is a bash script.

Full instructions:
https://github.com/napinoco/sedumipy/blob/main/docs/installation.rst\
"""


def build_bash_command(build_script: Path, out_path: Path) -> list[str]:
    """Non-Windows: run build_script directly (it's executable, with a
    #!/usr/bin/env bash shebang). Windows has no shebang support, so it
    needs `bash` invoked explicitly -- and by its *full path*, not the
    bare name "bash": subprocess.run() on Windows goes through
    CreateProcess(), which searches the System32 directory *before*
    PATH, and every Windows install since 10 1607 ships a `bash.exe`
    stub there that just prints "Windows Subsystem for Linux has no
    installed distributions" and exits nonzero -- it would silently
    shadow MSYS2's real bash.exe even with MSYS2 correctly first on
    PATH."""
    if sys.platform != "win32":
        return [str(build_script), str(out_path)]
    bash = shutil.which("bash")
    if bash is None:
        raise RuntimeError(WINDOWS_BUILD_HELP)
    return [bash, str(build_script), str(out_path)]

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

REPO_ROOT = Path(__file__).resolve().parent
if sys.platform == "darwin":
    LIB_NAME = "libsedumi.dylib"
elif sys.platform == "win32":
    LIB_NAME = "libsedumi.dll"
else:
    LIB_NAME = "libsedumi.so"


class BuildLibsedumi(build_ext):
    def build_extension(self, ext):
        if ext.name != "sedumipy._libsedumi_placeholder":
            super().build_extension(ext)
            return

        out_path = (Path(self.build_lib) / "sedumipy" / LIB_NAME).resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        build_script = REPO_ROOT / "tools" / "build_libsedumi.sh"
        command = build_bash_command(build_script, out_path)
        subprocess.run(command, check=True, cwd=REPO_ROOT)

    def get_ext_filename(self, fullname):
        if fullname == "sedumipy._libsedumi_placeholder":
            return str(Path(*fullname.split(".")[:-1]) / LIB_NAME)
        return super().get_ext_filename(fullname)


setup(
    ext_modules=[Extension("sedumipy._libsedumi_placeholder", sources=[])],
    cmdclass={"build_ext": BuildLibsedumi},
)

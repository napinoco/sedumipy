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
`bash` explicitly rather than executed directly (Windows has no shebang
support). This assumes an MSYS2 MinGW64 environment on PATH (`bash`,
`gcc`, and `mingw-w64-x86_64-openblas` -- see CONTRIBUTING.md's Windows
note); the resulting libsedumi.dll dynamically depends on
libopenblas.dll and a couple of mingw runtime DLLs that are not on a
plain end-user's Windows install, so a wheel built here is not yet
redistributable as-is -- .github/workflows/wheels.yml's Windows job runs
`delvewheel repair` as a post-build step to bundle them into the wheel.
"""

import subprocess
import sys
from pathlib import Path

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
        command = (
            ["bash", str(build_script), str(out_path)]
            if sys.platform == "win32"
            else [str(build_script), str(out_path)]
        )
        subprocess.run(command, check=True, cwd=REPO_ROOT)

    def get_ext_filename(self, fullname):
        if fullname == "sedumipy._libsedumi_placeholder":
            return str(Path(*fullname.split(".")[:-1]) / LIB_NAME)
        return super().get_ext_filename(fullname)


setup(
    ext_modules=[Extension("sedumipy._libsedumi_placeholder", sources=[])],
    cmdclass={"build_ext": BuildLibsedumi},
)

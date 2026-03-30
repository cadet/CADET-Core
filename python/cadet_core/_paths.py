from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Optional


def _package_root() -> Path:
    return Path(__file__).resolve().parent


def get_cadet_library_dir() -> Path:
    system = platform.system().lower()
    if system == "windows":
        return _package_root() / "bin"
    if system == "darwin":
        return _package_root() / ".dylibs"
    return _package_root() / "lib"


def get_cadet_binary() -> Path:
    bin_dir = _package_root() / "bin"
    exe = "cadet-cli.exe" if os.name == "nt" else "cadet-cli"
    p = bin_dir / exe
    if p.exists():
        return p
    raise FileNotFoundError(
        f"CADET-Core executable not found at {p}. "
        "Your CMake install step must place it into cadet_core/bin/"
    )


def get_cadet_library_path(prefer: Optional[str] = None) -> Path:
    lib_dir = get_cadet_library_dir()
    if prefer:
        p = lib_dir / prefer
        if p.exists():
            return p
        raise FileNotFoundError(f"Requested library '{prefer}' not found in {lib_dir}")

    system = platform.system().lower()
    if system == "windows":
        candidates = ["cadet.dll"]
    elif system == "darwin":
        candidates = ["libcadet.0.dylib", "libcadet.dylib"]
    else:
        candidates = ["libcadet.so"]

    for name in candidates:
        p = lib_dir / name
        if p.exists():
            return p

    listing = sorted([x.name for x in lib_dir.iterdir() if x.is_file()]) if lib_dir.exists() else []
    raise FileNotFoundError(
        "CADET-Core shared library not found in wheel. "
        f"Searched: {candidates}. Directory: {lib_dir}. Files present: {listing}"
    )


def _env_with_library_path(env: dict[str, str]) -> dict[str, str]:
    pkg_root = _package_root()
    bin_dir = pkg_root / "bin"

    if os.name == "nt":
        paths = [str(bin_dir)]

        # delvewheel vendors here by default: <site-packages>/<distribution>.libs
        # In your case it will be "cadet_core.libs" (distribution name can vary),
        # so search for any sibling directory ending with ".libs".
        site_packages = pkg_root.parent
        libs_dirs = sorted([p for p in site_packages.glob("*.libs") if p.is_dir()])

        # Prefer the one matching the package name if present
        preferred = site_packages / "cadet_core.libs"
        if preferred.is_dir():
            paths.insert(0, str(preferred))
        else:
            # otherwise add all .libs dirs (usually just one)
            paths = [*paths, *[str(p) for p in libs_dirs]]

        env["PATH"] = os.pathsep.join(paths) + os.pathsep + env.get("PATH", "")
        return env

    lib_dir = str(get_cadet_library_dir())
    key = "DYLD_LIBRARY_PATH" if platform.system().lower() == "darwin" else "LD_LIBRARY_PATH"
    env[key] = lib_dir + os.pathsep + env.get(key, "")
    return env



def run_cadet() -> int:
    exe = get_cadet_binary()
    env = _env_with_library_path(dict(os.environ))
    proc = subprocess.run([str(exe), *sys.argv[1:]], env=env)
    return int(proc.returncode)

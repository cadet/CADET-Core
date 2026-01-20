from __future__ import annotations

import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Optional


def _package_root() -> Path:
    # cadet_core package directory
    return Path(__file__).resolve().parent


def get_cadet_library_dir() -> Path:
    """
    Directory inside the installed wheel where CADET-Core shared libraries
    are expected to live.
    """
    # Convention: bundle runtime libs here
    return _package_root() / ".libs"


def get_cadet_binary() -> Path:
    """
    Path to the CADET-Core CLI executable shipped inside the wheel.

    Expected layout:
      cadet_core/bin/cadet (or cadet.exe on Windows)
    """
    bin_dir = _package_root() / "bin"
    exe = "cadet.exe" if os.name == "nt" else "cadet"
    return bin_dir / exe


def get_cadet_library_path(prefer: Optional[str] = None) -> Path:
    """
    Path to the main CADET-Core shared library shipped inside the wheel.

    prefer:
      Optional filename override if your library name differs from the defaults.
      Example: "libcadet.so" or "cadet-core.dll".
    """
    lib_dir = get_cadet_library_dir()
    if prefer:
        p = lib_dir / prefer
        if p.exists():
            return p
        raise FileNotFoundError(f"Requested library '{prefer}' not found in {lib_dir}")

    system = platform.system().lower()

    # These names are guesses. You should set the correct one once you know it.
    candidates: list[str]
    if system == "windows":
        candidates = ["cadet.dll", "CADET-Core.dll", "cadet-core.dll"]
    elif system == "darwin":
        candidates = ["libcadet.dylib", "libCADET-Core.dylib", "libcadet-core.dylib"]
    else:
        candidates = ["libcadet.so", "libCADET-Core.so", "libcadet-core.so"]

    for name in candidates:
        p = lib_dir / name
        if p.exists():
            return p

    # If nothing matched, provide a helpful error with directory listing
    listing = []
    if lib_dir.exists():
        listing = sorted([x.name for x in lib_dir.iterdir() if x.is_file()])
    raise FileNotFoundError(
        "CADET-Core shared library not found in wheel. "
        f"Searched: {candidates}. "
        f"Directory: {lib_dir}. "
        f"Files present: {listing}"
    )


def _env_with_library_path(env: dict[str, str]) -> dict[str, str]:
    """
    Prepare environment variables so the embedded cadet executable can find
    its bundled shared libraries.
    """
    lib_dir = str(get_cadet_library_dir())

    if os.name == "nt":
        # On Windows, PATH is used for DLL lookup
        env["PATH"] = lib_dir + os.pathsep + env.get("PATH", "")
        return env

    if platform.system().lower() == "darwin":
        key = "DYLD_LIBRARY_PATH"
    else:
        key = "LD_LIBRARY_PATH"

    env[key] = lib_dir + os.pathsep + env.get(key, "")
    return env


def run_cadet() -> int:
    """
    Console entry point. Runs the embedded CADET-Core CLI binary.
    """
    exe = get_cadet_binary()
    if not exe.exists():
        raise FileNotFoundError(
            f"CADET-Core executable not found at {exe}. "
            "Your CMake install step must place it into cadet_core/bin/"
        )

    env = _env_with_library_path(dict(os.environ))
    # Pass through all CLI args
    proc = subprocess.run([str(exe), *sys.argv[1:]], env=env)
    return int(proc.returncode)

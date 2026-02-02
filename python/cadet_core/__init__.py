from __future__ import annotations

from importlib.metadata import version as _pkg_version

from ._paths import (
    get_cadet_binary,
    get_cadet_library_dir,
    get_cadet_library_path,
)

__all__ = [
    "__version__",
    "get_cadet_binary",
    "get_cadet_library_dir",
    "get_cadet_library_path",
]

__version__ = _pkg_version("cadet-core")

# cadet_core/_ci/inspect_wheel.py
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
import glob
import importlib.util


def norm(p: str) -> str:
    return p.replace("\\", "/").lower()


def expand_args(args: list[str]) -> list[str]:
    expanded: list[str] = []
    for a in args:
        if any(ch in a for ch in ["*", "?", "["]):
            matches = glob.glob(a)
            if matches:
                expanded.extend(matches)
            else:
                expanded.append(a)
        else:
            expanded.append(a)
    return expanded


def inspect_installed() -> int:
    spec = importlib.util.find_spec("cadet_core")
    if spec is None or not spec.submodule_search_locations:
        print("ERROR: cadet_core is not importable in this environment")
        return 2

    pkg_root = Path(list(spec.submodule_search_locations)[0])
    bin_dir = pkg_root / "bin"

    print("\n" + "=" * 80)
    print("installed package inspection")
    print("=" * 80)
    print(f"cadet_core package root: {pkg_root}")
    print(f"cadet_core/bin exists: {bin_dir.exists()}")

    if bin_dir.exists():
        entries = sorted(p.name for p in bin_dir.iterdir())
        print("\n== cadet_core/bin entries (first 200) ==")
        for n in entries[:200]:
            print(" ", n)
        if len(entries) > 200:
            print(f"  ... ({len(entries) - 200} more)")

    cli_candidates = [bin_dir / "cadet-cli", bin_dir / "cadet-cli.exe"]
    has_cli = any(p.exists() for p in cli_candidates)

    print("\n== sanity ==")
    if not bin_dir.exists():
        print("WARNING: cadet_core/bin directory missing")
    if not has_cli:
        print("WARNING: cadet-cli not found (cadet-cli or cadet-cli.exe)")

    print("\nInspection OK")
    return 0


def inspect_wheel_file(wheel: str) -> int:
    wheel_path = Path(wheel)
    print("\n" + "=" * 80)
    print(f"wheel: {wheel_path}")
    print("=" * 80)

    if not wheel_path.exists():
        print(f"ERROR: wheel does not exist: {wheel_path}")
        return 2

    try:
        with zipfile.ZipFile(wheel_path) as z:
            names = z.namelist()
    except zipfile.BadZipFile:
        print(f"ERROR: not a valid wheel/zip file: {wheel_path}")
        return 2

    names_norm = [norm(n) for n in names]
    name_by_norm = dict(zip(names_norm, names))

    bin_like = [name_by_norm[n] for n in names_norm if n.startswith("cadet_core/bin/")]
    print("\n== cadet_core/bin-like entries (first 200) ==")
    for n in bin_like[:200]:
        print(" ", n)
    if len(bin_like) > 200:
        print(f"  ... ({len(bin_like) - 200} more)")

    cadet_dll = [name_by_norm[n] for n in names_norm if n.endswith("cadet_core/bin/cadet.dll")]
    cadet_cli = [
        name_by_norm[n]
        for n in names_norm
        if n.endswith("cadet_core/bin/cadet-cli.exe") or n.endswith("cadet_core/bin/cadet-cli")
    ]


    print("\n== core artifacts ==")
    print("cadet.dll entries:", cadet_dll if cadet_dll else "(none)")
    print("cadet-cli entries:", cadet_cli if cadet_cli else "(none)")

    exts = (".dll", ".pyd", ".so", ".dylib")
    dep_hits_norm = [n for n in names_norm if n.endswith(exts)]
    dep_hits = [name_by_norm[n] for n in dep_hits_norm]
    dep_hits_sorted = sorted(dep_hits, key=lambda s: norm(s))

    print("\n== binary payload in wheel ==")
    counts = defaultdict(int)
    for n in dep_hits_norm:
        for e in exts:
            if n.endswith(e):
                counts[e] += 1
                break
    for e in exts:
        print(f"{e}: {counts[e]}")

    print("\nFirst 200 binaries (by path):")
    for n in dep_hits_sorted[:200]:
        print(" ", n)
    if len(dep_hits_sorted) > 200:
        print(f"  ... ({len(dep_hits_sorted) - 200} more)")

    print("\n== sanity ==")
    if not bin_like:
        print("WARNING: no cadet_core/bin entries found at all")
    if not cadet_cli:
        print("WARNING: cadet-cli (exe or no-ext) not found under cadet_core/bin/")
    if any(n.endswith(".dll") for n in names_norm) and not cadet_dll:
        print("WARNING: wheel contains DLLs but cadet_core/bin/cadet.dll not found")

    print("\nInspection OK")
    return 0


def main() -> int:
    # No args: inspect installed package (works in CIBW_TEST_COMMAND)
    if len(sys.argv) == 1:
        return inspect_installed()

    wheels = expand_args(sys.argv[1:])
    rc = 0
    for w in wheels:
        this_rc = inspect_wheel_file(w)
        if this_rc != 0:
            rc = this_rc
    return rc


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import base64
import csv
import hashlib
import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


def b64u(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).rstrip(b"=").decode("ascii")


def sha256_urlsafe(data: bytes) -> str:
    return "sha256=" + b64u(hashlib.sha256(data).digest())


def norm(p: str) -> str:
    return p.replace("\\", "/")


def verify_record(whl: Path) -> None:
    with zipfile.ZipFile(whl, "r") as z:
        names = [norm(n) for n in z.namelist()]
        record_paths = [n for n in names if n.endswith(".dist-info/RECORD")]
        if not record_paths:
            raise SystemExit("RECORD not found during verification")
        record_path = record_paths[0]

        rows = {}
        for r in csv.reader(z.read(record_path).decode("utf-8").splitlines()):
            if not r:
                continue
            path = norm(r[0])
            h = r[1] if len(r) > 1 else ""
            size = r[2] if len(r) > 2 else ""
            rows[path] = (h, size)

        missing = [n for n in names if n not in rows]
        extra = [p for p in rows.keys() if p not in names]
        if missing or extra:
            raise SystemExit(
                f"RECORD path mismatch: missing={missing[:10]} extra={extra[:10]}"
            )

        bad = []
        for n in names:
            h, size = rows[n]
            data = z.read(n)
            if n == record_path:
                if h or size:
                    bad.append(f"{n}: RECORD entry must be empty hash/size")
                continue
            want_h = sha256_urlsafe(data)
            want_s = str(len(data))
            if h != want_h or size != want_s:
                bad.append(f"{n}: got ({h},{size}) want ({want_h},{want_s})")
                if len(bad) >= 20:
                    break

        if bad:
            raise SystemExit("RECORD mismatches:\n" + "\n".join(bad))

    print("RECORD verification OK")


def main() -> None:
    dest_dir = Path(os.environ["CIBW_DEST_DIR"].strip().strip('"')).resolve()
    wheels = sorted(dest_dir.glob("*.whl"), key=lambda p: p.stat().st_mtime)
    if not wheels:
        raise SystemExit(f"No wheels found in {dest_dir}")
    whl = wheels[-1]

    # Ensure wheel tooling is present
    subprocess.check_call([os.sys.executable, "-m", "pip", "install", "-U", "wheel"])

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # 1) Unpack the RECORD-consistent wheel created by delvewheel
        subprocess.check_call([os.sys.executable, "-m", "wheel", "unpack", str(whl), "-d", str(td)])
        unpacked = next(td.iterdir())  # e.g. td/cadet_core-5.1.5...

        # 2) Locate cadet_core/bin and *.libs in the unpacked tree
        cadet_pkg = unpacked / "cadet_core"
        bin_dir = cadet_pkg / "bin"
        if not bin_dir.is_dir():
            raise SystemExit(f"Expected bin dir not found: {bin_dir}")

        libs_dirs = [p for p in unpacked.iterdir() if p.is_dir() and p.name.endswith(".libs")]
        if not libs_dirs:
            raise SystemExit(f"No *.libs directory found in unpacked wheel at {unpacked}")
        libs_dir = libs_dirs[0]

        # 3) Copy DLLs from *.libs into cadet_core/bin
        dlls = sorted(libs_dir.glob("*.dll"))
        if not dlls:
            print(f"Warning: no DLLs found in {libs_dir}")
        for dll in dlls:
            target = bin_dir / dll.name
            if not target.exists():
                shutil.copy2(dll, target)

        # Optional: remove the *.libs directory to avoid duplicate DLLs in the wheel
        # If you still need it for Python imports, comment this out.
        shutil.rmtree(libs_dir, ignore_errors=True)

        # 4) Pack back into dest_dir (wheel pack regenerates RECORD)
        subprocess.check_call([os.sys.executable, "-m", "wheel", "pack", str(unpacked), "-d", str(dest_dir)])

    # Replace old wheel with newest packed one
    new_whl = sorted(dest_dir.glob("*.whl"), key=lambda p: p.stat().st_mtime)[-1]
    if new_whl != whl:
        whl.unlink()
        whl = new_whl

    print("Repacked OK:", whl)

    # 5) Verify RECORD matches final wheel contents
    verify_record(whl)


if __name__ == "__main__":
    main()


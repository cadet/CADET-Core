import os
import subprocess
import sys
import tempfile
from pathlib import Path

import cadet_core


def run(cmd, env):
    print("running:", " ".join(map(str, cmd)))
    r = subprocess.run(list(map(str, cmd)), env=env, capture_output=True, text=True)
    print("rc:", r.returncode)
    if r.stdout:
        print("stdout:\n", r.stdout[:4000])
    if r.stderr:
        print("stderr:\n", r.stderr[:4000])
    r.check_returncode()
    return r


def main():
    exe = Path(cadet_core.get_cadet_binary())
    bindir = exe.parent

    env = os.environ.copy()
    env["PATH"] = str(bindir) + os.pathsep + env.get("PATH", "")

    print("python:", sys.executable)
    print("platform:", sys.platform)
    print("cadet-core version:", cadet_core.__version__)
    print("cadet binary:", exe)
    print("bin dir:", bindir)

    # 1) basic version check
    run([exe, "--version"], env)

    # 2) generate a minimal input file with createLWE
    if sys.platform.startswith("win"):
        create_lwe = bindir / "createLWE.exe"
    else:
        create_lwe = bindir / "createLWE"

    if not create_lwe.exists():
        raise SystemExit(f"createLWE not found at {create_lwe}")

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        inp = td / "createLWE_smoke.h5"

        # createLWE requires -o <file>
        run([create_lwe, "-o", inp], env)

        if not inp.exists():
            raise SystemExit("createLWE did not produce input file")

        print("createLWE output file size:", inp.stat().st_size)

        # 3) run cadet-cli on the generated file
        run([exe, str(inp)], env)

    print("Smoke test completed successfully.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
EOC convergence study for radial LRM and LRMP DG models.
Runs polyDeg 1-5, nCells 2^0 - 2^7 (1 to 128).
"""

import json
import tempfile
import os
import subprocess
import sys

import h5py
import numpy as np

CADET_CLI = "/Users/yuvj/CADET/build/bin/cadet-cli"
TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def json_to_hdf5(json_data, filename):
    """Write JSON dict to HDF5 with /input/ prefix, matching CADET's expected format."""

    def write_group(h5grp, d):
        for key, val in d.items():
            if isinstance(val, dict):
                subgrp = h5grp.require_group(key)
                write_group(subgrp, val)
            elif isinstance(val, bool):
                h5grp.create_dataset(key, data=np.bool_(val))
            elif isinstance(val, (list, np.ndarray)):
                arr = np.array(val)
                h5grp.create_dataset(key, data=arr)
            elif isinstance(val, int):
                h5grp.create_dataset(key, data=np.int64(val))
            elif isinstance(val, float):
                h5grp.create_dataset(key, data=np.float64(val))
            elif isinstance(val, str):
                h5grp.create_dataset(key, data=np.bytes_(val))

    with h5py.File(filename, "w") as f:
        inp = f.create_group("input")
        write_group(inp, json_data)


def run_simulation(json_path, polydeg, nelem, n_time_points=200):
    """Run a single CADET simulation, return outlet time series or None on failure."""
    with open(json_path) as f:
        cfg = json.load(f)

    # Set discretization
    unit_key = "unit_001"
    cfg["model"][unit_key]["discretization"]["POLYDEG"] = polydeg
    cfg["model"][unit_key]["discretization"]["NELEM"] = nelem

    # Set output times
    t_end = cfg["solver"]["sections"]["SECTION_TIMES"][-1]
    times = np.linspace(t_end / n_time_points, t_end, n_time_points).tolist()
    cfg["solver"]["USER_SOLUTION_TIMES"] = times

    # Add required solver fields that the C++ test framework provides by default
    cfg["solver"].setdefault("CONSISTENT_INIT_MODE", 1)
    cfg["solver"].setdefault("NTHREADS", 1)

    # Ensure outlet is written
    cfg.setdefault("return", {}).setdefault(unit_key, {})["WRITE_SOLUTION_OUTLET"] = True

    # Write to temp HDF5
    tmpfile = tempfile.mktemp(suffix=".h5")
    try:
        json_to_hdf5(cfg, tmpfile)

        result = subprocess.run(
            [CADET_CLI, tmpfile],
            capture_output=True, text=True, timeout=600
        )

        if result.returncode != 0:
            print(f"    FAILED (exit {result.returncode}): {result.stderr.strip()[:200]}", file=sys.stderr)
            return None

        # Read outlet solution (may be split by component)
        with h5py.File(tmpfile, "r") as f:
            sol_grp = f[f"/output/solution/{unit_key}"]
            if "SOLUTION_OUTLET" in sol_grp:
                outlet = sol_grp["SOLUTION_OUTLET"][:]
            else:
                # Components split into SOLUTION_OUTLET_COMP_000, etc.
                comp_keys = sorted([k for k in sol_grp if k.startswith("SOLUTION_OUTLET_COMP_")])
                outlet = np.column_stack([sol_grp[k][:] for k in comp_keys])
        return outlet.flatten()
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT (polydeg={polydeg}, nelem={nelem})", file=sys.stderr)
        return None
    finally:
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)


def run_eoc_study(model_name, json_path, poly_degs, n_cells_list):
    """Run full EOC study and print results."""
    print(f"\n{'=' * 80}")
    print(f"  EOC Convergence Study: {model_name}")
    print(f"{'=' * 80}")

    for p in poly_degs:
        print(f"\n--- polyDeg = {p} (expected order ~ {p + 1}) ---")

        header = f"{'nCells':>8s} | {'nDOFs':>8s} | {'L2 error':>14s} | {'EOC':>8s}"
        print(header)
        print("-" * len(header))

        prev_outlet = None
        errors = []

        for nc in n_cells_list:
            ndofs = nc * (p + 1)
            outlet = run_simulation(json_path, p, nc)

            if outlet is None:
                print(f"{nc:8d} | {ndofs:8d} | {'FAILED':>14s} | {'---':>8s}")
                prev_outlet = None
                errors.append(None)
                continue

            if prev_outlet is not None and len(outlet) == len(prev_outlet):
                diff = prev_outlet - outlet
                err = np.sqrt(np.mean(diff ** 2))
                errors.append(err)

                # Compute EOC
                if len(errors) >= 2 and errors[-2] is not None and err > 1e-15:
                    eoc = np.log2(errors[-2] / err)
                    print(f"{nc:8d} | {ndofs:8d} | {err:14.6e} | {eoc:8.3f}")
                else:
                    print(f"{nc:8d} | {ndofs:8d} | {err:14.6e} | {'---':>8s}")
            else:
                errors.append(None)
                print(f"{nc:8d} | {ndofs:8d} | {'(ref)':>14s} | {'---':>8s}")

            prev_outlet = outlet

        print()


if __name__ == "__main__":
    poly_degs = list(range(1, 6))           # 1 to 5
    n_cells_list = [2**k for k in range(8)] # 1, 2, 4, 8, 16, 32, 64, 128

    models = [
        ("Radial LRM DG", os.path.join(TEST_DIR, "model_radLRM_DG_gaussianPulse_1comp_eocbenchmark.json")),
        ("Radial LRMP DG", os.path.join(TEST_DIR, "model_radLRMP_DG_gaussianPulse_1comp_eocbenchmark.json")),
    ]

    for name, path in models:
        run_eoc_study(name, path, poly_degs, n_cells_list)

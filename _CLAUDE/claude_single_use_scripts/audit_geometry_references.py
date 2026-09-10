"""Audit the references of every physical setting of scripts/verify_geometries.py.

Resolves, through the very code the study uses, which reference each setting would
take and whether that file exists. Nothing is generated and nothing is simulated:
geometry_references.resolve is called with generate=False.

Run from the CADET-Verification root.
"""

import os
from pathlib import Path

from src import analytical
from src import bench_func
from src import geometry_references
from src.benchmark_models import setting_Col1D_pureTransport_1comp_benchmark1 as transport_setting
from src.benchmark_models import setting_Col1D_lin_1comp_benchmark1 as lin_setting
from src.benchmark_models.setting_Col1D_pureTransport_1comp_benchmark1 import (
    create_convergence_object
)
from functools import partial

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
os.environ['PATH'] = (
    str(Path(CADET_PATH) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

DATA = str(Path.cwd() / 'data')

DEFAULT_N_LEVELS = {0: 13, 1: 13, 2: 12, 3: 11, 4: 10}

# setting name, kind, geometry, extra transport kwargs, ax_methods, active in the script
SETTINGS = [
    ('radialAdvDPFR_1comp_benchmark1', 'transport', 'RADIAL_FLOW_CYLINDER_SHELL',
     dict(advection=True, dispersion=False), [0, 1, 2, 3, 4], False),
    ('radialDispDPFR_1comp_benchmark1', 'transport', 'RADIAL_FLOW_CYLINDER_SHELL',
     dict(advection=False, dispersion=True), [0, 1, 2, 3, 4], False),
    ('radialDPFR_1comp_benchmark1', 'transport', 'RADIAL_FLOW_CYLINDER_SHELL',
     dict(advection=True, dispersion=True), [0, 1, 2, 3, 4], False),
    ('radialLRMP_dynLin_1comp_benchmark1', 'lrmp', 'RADIAL_FLOW_CYLINDER_SHELL',
     {}, [0, 1, 2, 3, 4], False),
    ('frustumAdvDPFR_1comp_benchmark1', 'transport', 'AXIAL_FLOW_FRUSTUM',
     dict(advection=True, dispersion=False), [0, 1, 2, 3, 4], False),
    ('frustumDispDPFR_1comp_benchmark1', 'transport', 'AXIAL_FLOW_FRUSTUM',
     dict(advection=False, dispersion=True), [0, 1, 2, 3, 4], False),
    ('frustumDPFR_1comp_benchmark1', 'transport', 'AXIAL_FLOW_FRUSTUM',
     dict(advection=True, dispersion=True), [0, 1, 2, 3, 4], False),
    ('frustumLRMP_dynLin_1comp_benchmark1', 'lrmp', 'AXIAL_FLOW_FRUSTUM',
     {}, [0, 1, 2, 3, 4], False),
    ('smoothlyVaryingAdvDPFR_1comp_benchmark1', 'transport', 'SMOOTHLY_VARYING',
     dict(advection=True, dispersion=False), [1, 2, 3, 4], True),
    ('smoothlyVaryingDispDPFR_1comp_benchmark1', 'transport', 'SMOOTHLY_VARYING',
     dict(advection=False, dispersion=True), [1, 2, 3, 4], True),
    ('smoothlyVaryingDPFR_1comp_benchmark1', 'transport', 'SMOOTHLY_VARYING',
     dict(advection=True, dispersion=True), [1, 2, 3, 4], True),
]

missing = []

for name, kind, geometry, extra, ax_methods, active in SETTINGS:
    tag = "ACTIVE  " if active else "disabled"
    analytic = analytical.load_reference(
        name, os.path.join(DATA, geometry_references.ANALYTICAL_SUBDIR))

    print(f"[{tag}] {name}")

    if analytic is not None:
        print("           analytical reference present -> used for every method")
        print()
        continue

    print("           no analytical reference -> needs one CADET reference per method")

    n_levels = [DEFAULT_N_LEVELS[m] for m in ax_methods]
    ax_discs = [bench_func.disc_list(1, n) for n in n_levels]

    for method, discs in zip(ax_methods, ax_discs):
        n_elem = int((8 * discs[-1]) if kind == 'transport' else discs[-1])
        path = geometry_references.cadet_reference_path(name, method, n_elem, DATA)
        exists = os.path.exists(path)
        mark = "found  " if exists else "MISSING"
        print(f"           {mark}  {os.path.basename(path)}")
        if not exists:
            missing.append((name, kind, geometry, extra, method, n_elem, discs, active))
    print()

print("=" * 78)
print(f"{len(missing)} missing reference files")
for name, kind, geometry, extra, method, n_elem, discs, active in missing:
    print(f"   {'ACTIVE  ' if active else 'disabled'} {name}  method {method}  {n_elem} elements")

# what is present but not requested by the current sweep lengths
print()
print("Present CADET references that no setting above asks for:")
requested = set()
for name, kind, geometry, extra, ax_methods, active in SETTINGS:
    if analytical.load_reference(name, os.path.join(DATA, geometry_references.ANALYTICAL_SUBDIR)) is not None:
        continue
    n_levels = [DEFAULT_N_LEVELS[m] for m in ax_methods]
    for method, n in zip(ax_methods, n_levels):
        discs = bench_func.disc_list(1, n)
        n_elem = int((8 * discs[-1]) if kind == 'transport' else discs[-1])
        requested.add(os.path.basename(
            geometry_references.cadet_reference_path(name, method, n_elem, DATA)))

ref_dir = Path(DATA) / geometry_references.CADET_SUBDIR
for f in sorted(p.name for p in ref_dir.glob('*.h5')):
    if f not in requested:
        print("   " + f)

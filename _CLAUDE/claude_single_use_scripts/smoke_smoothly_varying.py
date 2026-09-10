"""Smoke test of the SMOOTHLY_VARYING DG geometry.

Test 1 (equivalence): a SMOOTHLY_VARYING column whose nodal areas are the exact
frustum areas pi*r(x)^2 must reproduce the AXIAL_FLOW_FRUSTUM result, since the
quadratic area is represented exactly by the nodal interpolant for polyDeg >= 2.

Test 2 (input validation): wrong length, non-positive entries and a jump at an
element interface must all be rejected.

Run from the CADET-Verification root.
"""

import copy
import sys
from pathlib import Path

import numpy as np
from addict import Dict
from cadet import Cadet

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))

import src.benchmark_models.setting_Col1D_pureTransport_1comp_benchmark1 as setting
from src.benchmark_models.setting_Col1D_lin_1comp_benchmark1 import (
    get_column_geometry_configuration
)

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
OUT = Path(r"C:\Users\jmbr\AppData\Local\Temp\claude\smoke_smoothly_varying")
OUT.mkdir(parents=True, exist_ok=True)

POLYDEG = 3
NELEM = 8


def run(config, name):
    model = Cadet(install_path=CADET_PATH)
    model.root.input = copy.deepcopy(config['input'])
    model.filename = str(OUT / (name + '.h5'))
    model.save()
    data = model.run_simulation()
    if data.return_code != 0:
        raise RuntimeError(f"{name} failed: {data.error_message}\n{data.log}")
    model.load()
    return model


def bulk(model):
    return np.asarray(model.root.output.solution.unit_001.solution_bulk)


# ---------------------------------------------------------------- test 1
frustum_geom = get_column_geometry_configuration('AXIAL_FLOW_FRUSTUM')
r0 = np.sqrt(frustum_geom['cross_section_area_large_end'] / np.pi)
rL = np.sqrt(frustum_geom['cross_section_area_small_end'] / np.pi)


def frustum_area(xi):
    xi = np.asarray(xi, dtype=float)
    return np.pi * (r0 + (rL - r0) * xi) ** 2


common = dict(spatial_method_bulk=POLYDEG, refinement=1, advection=True,
              dispersion=True, write_solution_bulk=True)

cfg_frustum = setting.get_model(column_geometry='AXIAL_FLOW_FRUSTUM', **common)
cfg_smooth = setting.get_model(column_geometry='SMOOTHLY_VARYING',
                               area_profile=frustum_area, **common)

areas = cfg_smooth['input']['model']['unit_001']['cross_sectional_area_at_nodes']
print(f"n areas = {len(areas)} (expected {NELEM * (POLYDEG + 1)}), "
      f"first = {areas[0]:.6e}, last = {areas[-1]:.6e}")

sol_frustum = bulk(run(cfg_frustum, 'frustum'))
sol_smooth = bulk(run(cfg_smooth, 'smooth_as_frustum'))

diff = np.max(np.abs(sol_frustum - sol_smooth))
scale = np.max(np.abs(sol_frustum))
print(f"TEST 1 frustum equivalence: max abs diff = {diff:.3e}, "
      f"relative = {diff / scale:.3e}")

# ---------------------------------------------------------------- test 2
def expect_failure(label, mutate):
    cfg = setting.get_model(column_geometry='SMOOTHLY_VARYING', **common)
    areas = list(cfg['input']['model']['unit_001']['cross_sectional_area_at_nodes'])
    cfg['input']['model']['unit_001']['cross_sectional_area_at_nodes'] = mutate(areas)
    try:
        run(cfg, 'invalid')
    except RuntimeError as exc:
        first = [line for line in str(exc).splitlines() if 'CROSS_SECTIONAL' in line
                 or 'SMOOTHLY' in line]
        print(f"TEST 2 {label}: rejected -> {first[0].strip() if first else 'see log'}")
        return
    print(f"TEST 2 {label}: NOT REJECTED (problem)")


expect_failure('wrong length', lambda a: a[:-1])
expect_failure('non-positive entry', lambda a: [-1.0] + a[1:])
expect_failure('jump at interface',
               lambda a: a[:POLYDEG + 1] + [a[POLYDEG + 1] * 1.5] + a[POLYDEG + 2:])

# ---------------------------------------------------------------- profile info
sine = setting.get_cross_sectional_area_profile()
xi = np.linspace(0.0, 1.0, 11)
print("sine profile A(xi) / A(0):", np.array2string(sine(xi) / sine(0.0),
                                                    precision=4))

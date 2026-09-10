"""Dry run of the reference resolution for the geometry studies.

Calls the real benchmark builders with generate_references=False and prints which
reference each physical setting resolves to and whether the sweep is kept intact.

Run from the CADET-Verification root.
"""

import os
from pathlib import Path

from src import bench_configs

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
os.environ['PATH'] = (
    str(Path(CADET_PATH) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

DATA = str(Path.cwd() / 'data')
USER_TIMES = [5.0]

TRANSPORT = [
    ('radialAdvDPFR_1comp_benchmark1', 'RADIAL_FLOW_CYLINDER_SHELL', True, False, [0, 1, 2, 3, 4]),
    ('radialDispDPFR_1comp_benchmark1', 'RADIAL_FLOW_CYLINDER_SHELL', False, True, [0, 1, 2, 3, 4]),
    ('radialDPFR_1comp_benchmark1', 'RADIAL_FLOW_CYLINDER_SHELL', True, True, [0, 1, 2, 3, 4]),
    ('frustumAdvDPFR_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM', True, False, [0, 1, 2, 3, 4]),
    ('frustumDispDPFR_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM', False, True, [0, 1, 2, 3, 4]),
    ('frustumDPFR_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM', True, True, [0, 1, 2, 3, 4]),
    ('smoothlyVaryingAdvDPFR_1comp_benchmark1', 'SMOOTHLY_VARYING', True, False, [1, 2, 3, 4]),
    ('smoothlyVaryingDispDPFR_1comp_benchmark1', 'SMOOTHLY_VARYING', False, True, [1, 2, 3, 4]),
    ('smoothlyVaryingDPFR_1comp_benchmark1', 'SMOOTHLY_VARYING', True, True, [1, 2, 3, 4]),
]

LRMP = [
    ('radialLRMP_dynLin_1comp_benchmark1', 'RADIAL_FLOW_CYLINDER_SHELL'),
    ('frustumLRMP_dynLin_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM'),
]

for name, geometry, advection, dispersion, ax_methods in TRANSPORT:
    cfg = bench_configs.paper_geometry_transport_benchmark(
        setting_name=name, small_test=False, ref_filepath=DATA,
        generate_references=False, cadet_path=CADET_PATH,
        user_solution_times_unit_state=USER_TIMES,
        **{'advection': advection, 'dispersion': dispersion,
           'ax_methods': ax_methods, 'column_geometry': geometry}
    )
    n_refs = len({id(r) for r in cfg['ref_files'][0]})
    levels = [len(d) for d in cfg['ax_discs'][0]]
    print(f"    distinct references: {n_refs}, sweep levels per method: {levels}")
    print()

for name, geometry in LRMP:
    cfg = bench_configs.paper_geometry_LRMPdynLin_benchmark(
        setting_name=name, small_test=False, ref_filepath=DATA,
        generate_references=False, cadet_path=CADET_PATH,
        user_solution_times_unit_state=USER_TIMES,
        **{'column_geometry': geometry}
    )
    n_refs = len({id(r) for r in cfg['ref_files'][0]})
    levels = [len(d) for d in cfg['ax_discs'][0]]
    print(f"    distinct references: {n_refs}, sweep levels per method: {levels}")
    print()

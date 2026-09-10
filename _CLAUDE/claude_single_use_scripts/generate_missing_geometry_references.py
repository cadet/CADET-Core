"""Generate the missing CADET references of the varying cross section settings.

Uses the study's own reference resolution (bench_configs.paper_geometry_transport_benchmark
-> geometry_references.resolve with generate=True), so each physical setting gets exactly
one reference, at the highest DG degree of the study and one refinement step beyond its
finest sweep level, under the name the study looks up. Existing files are never
overwritten: resolve only generates what is missing.

Run from the CADET-Verification root.
"""

import os
import time
from pathlib import Path

from src import bench_configs

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
os.environ['PATH'] = (
    str(Path(CADET_PATH) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

DATA = str(Path.cwd() / 'data')

# same arguments as scripts/verify_geometries.py
SMALL_TEST = False
USER_SOLUTION_TIMES_UNIT_STATE = [5.0]
AX_METHODS = [1, 2, 3, 4]

SETTINGS = [
    ('smoothlyVaryingAdvDPFR_1comp_benchmark1', dict(advection=True, dispersion=False)),
    ('smoothlyVaryingDispDPFR_1comp_benchmark1', dict(advection=False, dispersion=True)),
    ('smoothlyVaryingDPFR_1comp_benchmark1', dict(advection=True, dispersion=True)),
]

for name, transport in SETTINGS:
    print(f"=== {name} ===", flush=True)
    start = time.time()
    bench_configs.paper_geometry_transport_benchmark(
        setting_name=name,
        small_test=SMALL_TEST, ref_filepath=DATA,
        generate_references=True, cadet_path=CADET_PATH,
        user_solution_times_unit_state=USER_SOLUTION_TIMES_UNIT_STATE,
        **{**transport, 'ax_methods': AX_METHODS, 'column_geometry': 'SMOOTHLY_VARYING'}
    )
    print(f"--- {name} done in {time.time() - start:.1f} s", flush=True)

print("all references present", flush=True)

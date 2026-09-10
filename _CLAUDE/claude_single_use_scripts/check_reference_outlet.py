"""Check that the SMOOTHLY_VARYING reference test compares a meaningful outlet signal.

Prints the magnitude of the stored reference outlet and, as a sensitivity check, the
deviation caused by a 1 permille change of the prescribed cross section area.
"""

import copy
import json
import os
from pathlib import Path

import h5py
import numpy as np
from cadet import Cadet

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
DATA_DIR = Path(r"C:\Users\jmbr\software\CADET-Core\test\data")
NAME = "smoothlyVaryingCOL1D_transport_1comp_benchmark1"
SCRATCH = Path(r"C:\Users\jmbr\AppData\Local\Temp\claude\refcheck")
SCRATCH.mkdir(parents=True, exist_ok=True)

os.environ['PATH'] = (
    str(Path(CADET_PATH) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

with h5py.File(DATA_DIR / ("ref_" + NAME + "_DG_P3Z8.h5"), 'r') as f:
    ref_outlet = np.array(f['/output/solution/unit_001/SOLUTION_OUTLET']).ravel()
    times = np.array(f['/output/solution/SOLUTION_TIMES']).ravel()

print(f"reference outlet: {len(ref_outlet)} points, "
      f"max = {np.max(ref_outlet):.6e} at t = {times[np.argmax(ref_outlet)]:.2f} s, "
      f"mean |.| = {np.mean(np.abs(ref_outlet)):.6e}")
print(f"points above 1e-6: {int(np.sum(np.abs(ref_outlet) > 1e-6))}")

# sensitivity: perturb the prescribed areas by 0.1 percent and rerun
setup = json.load(open(DATA_DIR / ("config_" + NAME + ".json")))
setup['model']['unit_001']['CROSS_SECTIONAL_AREA_AT_NODES'] = [
    1.001 * a for a in setup['model']['unit_001']['CROSS_SECTIONAL_AREA_AT_NODES']
]

model = Cadet(install_path=CADET_PATH)
model.root.input = copy.deepcopy(setup)
model.filename = str(SCRATCH / "perturbed.h5")
model.save()
data = model.run_simulation()
if data.return_code != 0:
    raise RuntimeError(data.error_message)
model.load()

perturbed = np.asarray(model.root.output.solution.unit_001.solution_outlet).ravel()
scale = np.max(np.abs(ref_outlet))
print(f"0.1% area change moves the outlet by max {np.max(np.abs(perturbed - ref_outlet)):.6e} "
      f"({np.max(np.abs(perturbed - ref_outlet)) / scale:.3e} relative to the peak)")

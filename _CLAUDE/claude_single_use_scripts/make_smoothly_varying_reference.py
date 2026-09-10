"""Generate the CADET-Core reference test files for the SMOOTHLY_VARYING geometry.

Writes
    test/data/config_smoothlyVaryingCOL1D_transport_1comp_benchmark1.json
    test/data/ref_smoothlyVaryingCOL1D_transport_1comp_benchmark1_DG_P3Z8.h5

using exactly the combined convection-dispersion setting of the EOC study, i.e.
the sine cross section area profile of
setting_Col1D_pureTransport_1comp_benchmark1.get_cross_sectional_area_profile.

Run from the CADET-Verification root.
"""

import copy
import json
import os
from pathlib import Path

from cadet import Cadet

import src.benchmark_models.setting_Col1D_pureTransport_1comp_benchmark1 as setting

CADET_PATH = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"
DATA_DIR = Path(r"C:\Users\jmbr\software\CADET-Core\test\data")
NAME = "smoothlyVaryingCOL1D_transport_1comp_benchmark1"
REF_NAME = NAME + "_DG_P3Z8"

POLYDEG = 3
REFINEMENT = 1  # -> NELEM = 8

os.environ['PATH'] = (
    str(Path(CADET_PATH) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)


def to_cadet_json(node):
    """Convert the addict configuration to the key convention of test/data configs.

    Leaf fields are upper case, groups keep their lower case name.
    """
    out = {}
    for key, value in node.items():
        if isinstance(value, dict):
            out[key] = to_cadet_json(value)
        else:
            out[key.upper()] = value.tolist() if hasattr(value, 'tolist') else value
    return out


config = setting.get_model(
    spatial_method_bulk=POLYDEG,
    refinement=REFINEMENT,
    column_geometry='SMOOTHLY_VARYING',
    advection=True,
    dispersion=True,
    write_solution_bulk=0,
)

disc = config['input']['model']['unit_001']['discretization']
print(f"POLYDEG = {disc['POLYDEG']}, NELEM = {disc['NELEM']}")
areas = config['input']['model']['unit_001']['cross_sectional_area_at_nodes']
print(f"{len(areas)} nodal areas, min = {min(areas):.6e}, max = {max(areas):.6e}")

# reference simulation: the same setting, run and stored with its input group
model = Cadet(install_path=CADET_PATH)
model.root.input = copy.deepcopy(config['input'])
model.filename = str(DATA_DIR / ("ref_" + REF_NAME + ".h5"))
model.save()
data = model.run_simulation()
if data.return_code != 0:
    raise RuntimeError(f"reference simulation failed: {data.error_message}\n{data.log}")
print("wrote " + model.filename)

# model setup file of the test
setup = to_cadet_json(config['input'])
with open(DATA_DIR / ("config_" + NAME + ".json"), 'w') as f:
    json.dump(setup, f, indent=4)
print("wrote " + str(DATA_DIR / ("config_" + NAME + ".json")))

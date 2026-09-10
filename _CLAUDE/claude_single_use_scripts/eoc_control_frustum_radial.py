"""Control run for the SMOOTHLY_VARYING EOC study.

Runs the pure dispersion and the combined transport settings of the already
verified frustum and radial geometries with exactly the same reduced sweeps,
the same self-convergence references (the temporary data directory hides the
stored analytical references) and the same error measure as
eoc_smoothly_varying.py. Any order reduction that shows up here as well is a
property of the DG operator, not of the new geometry.

Run from the CADET-Verification root.
"""

import os
from pathlib import Path

from src import bench_configs
from src import bench_func

cadet_path = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"

os.environ['PATH'] = (
    str(Path(cadet_path) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

scratch = Path(r"C:\Users\jmbr\AppData\Local\Temp\claude\eoc_control")
output_path = scratch / "output"
reference_data_path = str(scratch / "data")
os.makedirs(output_path / "chromatography", exist_ok=True)
os.makedirs(reference_data_path, exist_ok=True)

user_solution_times_unit_state = [5.0]

ax_methods = [1, 2, 3, 4]
n_levels = [8, 7, 6, 5]

cadet_configs = []
cadet_config_names = []
include_sens = []
ref_files = []
unit_IDs = []
which = []
idas_abstol = []
methods = []
discs = []
par_methods = []
par_discs = []
disc_refinement_functions = []

settings = [
    # the setting name has to carry the transport model token (DPFR), see
    # convergence.recalculate_results
    ('frustumDispSelfrefDPFR_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM', False, True),
    ('radialDispSelfrefDPFR_1comp_benchmark1', 'RADIAL_FLOW_CYLINDER_SHELL', False, True),
    ('frustumBothSelfrefDPFR_1comp_benchmark1', 'AXIAL_FLOW_FRUSTUM', True, True),
]

for setting_name, geometry, advection, dispersion in settings:
    addition = bench_configs.paper_geometry_transport_benchmark(
        setting_name=setting_name,
        small_test=False, ref_filepath=reference_data_path,
        generate_references=True, cadet_path=cadet_path,
        user_solution_times_unit_state=user_solution_times_unit_state,
        **{'advection': advection, 'dispersion': dispersion,
           'ax_methods': ax_methods, 'n_levels': n_levels,
           'column_geometry': geometry}
    )

    bench_configs.add_benchmark(
        cadet_configs, include_sens, ref_files, unit_IDs, which,
        methods, discs, par_methods, par_discs, idas_abstol=idas_abstol,
        cadet_config_names=cadet_config_names, addition=addition,
        disc_refinement_functions=disc_refinement_functions,
    )

bench_func.run_convergence_analysis(
    output_path=output_path / "chromatography",
    cadet_path=cadet_path,
    cadet_configs=cadet_configs,
    cadet_config_names=cadet_config_names,
    include_sens=include_sens,
    ref_files=ref_files,
    unit_IDs=unit_IDs,
    which=which,
    ax_methods=methods,
    ax_discs=discs,
    par_methods=par_methods,
    par_discs=par_discs,
    idas_abstol=idas_abstol,
    n_jobs=-1,
    rerun_sims=1,
    disc_refinement_functions=disc_refinement_functions,
    time_point=0,
)

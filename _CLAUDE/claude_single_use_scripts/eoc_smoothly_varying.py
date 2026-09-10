"""Reduced EOC study of the SMOOTHLY_VARYING DG geometry.

Same three settings that scripts/verify_geometries.py adds for this geometry
(pure advection, pure dispersion, combined), but with shortened refinement
sweeps so that the study runs in minutes rather than the paper-scale hours.
References and simulation output go to a temporary directory, so nothing is
written into the CADET-Verification repository.

Run from the CADET-Verification root.
"""

import os
from pathlib import Path

import src.utility.convergence as convergence
from src import bench_configs
from src import bench_func

cadet_path = r"C:\Users\jmbr\software\CADET-Core\out\install\aRELEASE"

# The benchmark helpers construct Cadet objects without an explicit install path and
# rely on cadet-cli being found on PATH.
os.environ['PATH'] = (
    str(Path(cadet_path) / 'bin') + os.pathsep + os.environ.get('PATH', '')
)

scratch = Path(r"C:\Users\jmbr\AppData\Local\Temp\claude\eoc_smoothly_varying")
output_path = scratch / "output"
reference_data_path = str(scratch / "data")
os.makedirs(output_path / "chromatography", exist_ok=True)
os.makedirs(reference_data_path, exist_ok=True)

user_solution_times_unit_state = [5.0]

ax_methods = [1, 2, 3, 4]
# Element counts are 8 * 2^level, i.e. the finest (reference) levels are
# 1024, 512, 256 and 128 elements for the polynomial degrees 1 to 4.
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
    ('smoothlyVaryingAdvDPFR_1comp_benchmark1', True, False),
    ('smoothlyVaryingDispDPFR_1comp_benchmark1', False, True),
    ('smoothlyVaryingDPFR_1comp_benchmark1', True, True),
]

for setting_name, advection, dispersion in settings:
    addition = bench_configs.paper_geometry_transport_benchmark(
        setting_name=setting_name,
        small_test=False, ref_filepath=reference_data_path,
        generate_references=True, cadet_path=cadet_path,
        user_solution_times_unit_state=user_solution_times_unit_state,
        **{'advection': advection, 'dispersion': dispersion,
           'ax_methods': ax_methods, 'n_levels': n_levels,
           'column_geometry': 'SMOOTHLY_VARYING'}
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

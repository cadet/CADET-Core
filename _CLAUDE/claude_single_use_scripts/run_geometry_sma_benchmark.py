"""Run one refinement level of the SMA-LRMP column geometry performance benchmark.

Reproduces the IDAS step counts and outlet L_inf errors of the CADET-Verification
sweep (scripts/verify_geometries.py, --run-performance-sma-tests with
--sma-particle-resolutions HOMOGENEOUS_PARTICLE) for a single spatial resolution,
without going through cadetrdm. Each geometry is simulated with the CADET-Core
build given on the command line and compared against the stored DG P5 Z256
reference of CADET-Verification.

Usage:
  python run_geometry_sma_benchmark.py <cadet-cli path> <output dir> [n_cells] [geometries]
"""

import os
import sys
from pathlib import Path

VERIFICATION_ROOT = Path(r"C:\Users\jmbr\software\CADET-Verification")

sys.path.insert(0, str(VERIFICATION_ROOT))
os.chdir(VERIFICATION_ROOT)

import src.utility.convergence as convergence  # noqa: E402
from src.bench_configs import SMA_performance_benchmark  # noqa: E402
from src.bench_configs import geometry_reference_path  # noqa: E402
from src.bench_func import create_object_from_config  # noqa: E402
from src.bench_func import run_simulation_in_verification  # noqa: E402

GEOMETRIES = {
    'axial': 'AXIAL_FLOW_CYLINDER',
    'frustum': 'AXIAL_FLOW_FRUSTUM',
    'radial': 'RADIAL_FLOW_CYLINDER_SHELL',
    # The frustum geometry with both ends of the same size: the same column as the
    # axial cylinder, but run through the cross-section weighted (weight exponent 2,
    # zero slope) reconstruction instead of the equidistant axial one
    'frustumConstR': 'AXIAL_FLOW_FRUSTUM',
}


def make_radius_constant(config_data):
    """Set the small end of the frustum to the size of its large end."""

    unit = config_data['input']['model']['unit_000']
    keys = {key.lower(): key for key in unit}
    unit[keys['cross_section_area_small_end']] = \
        unit[keys['cross_section_area_large_end']]


def main():
    cadet_path = sys.argv[1]
    output_path = sys.argv[2]

    # CADET-Python autodetects an installation when a Cadet object is created, before
    # run_simulation_in_verification gets to set the one given here
    os.environ['PATH'] = str(Path(cadet_path).parent) + os.pathsep + os.environ['PATH']
    n_cells = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    geometries = sys.argv[4].split(',') if len(sys.argv) > 4 else list(GEOMETRIES)

    os.makedirs(output_path, exist_ok=True)
    reference_data_path = str(VERIFICATION_ROOT / 'data')

    for name in geometries:
        geometry = GEOMETRIES[name]

        # The benchmark configuration carries the model and the name the sweep
        # gives its simulations; the tolerance is the one of the sweep (1e-12 for
        # the homogeneous particle)
        benchmark = SMA_performance_benchmark(
            column_geometry=geometry, particle_type='HOMOGENEOUS_PARTICLE')

        setting_name = benchmark['cadet_config_names'][0]
        if name == 'frustumConstR':
            make_radius_constant(benchmark['cadet_config_jsons'][0])
            setting_name = setting_name.replace('frustum_', 'frustumConstR_')

        simulation = create_object_from_config(
            config_data=benchmark['cadet_config_jsons'][0],
            setting_name=setting_name,
            unit_id=benchmark['unit_IDs'][0],
            ax_method=0,
            ax_cells=n_cells,
            par_method=None,
            par_cells=None,
            output_path=output_path,
            idas_abstol=benchmark['idas_abstol'][0][0],
            include_sens=False,
        )

        print(f"{name}: running {Path(simulation.filename).name}", flush=True)
        run_simulation_in_verification([simulation], cadet_path)

        steps = convergence.get_idas_timesteps(simulation.filename)
        seconds = convergence.get_compute_time(simulation.filename)

        # The constant radius column is the axial cylinder, so that is the reference
        # it is measured against
        reference = geometry_reference_path(
            'sma_lrmp',
            'AXIAL_FLOW_CYLINDER' if name == 'frustumConstR' else geometry,
            reference_data_path)
        if os.path.exists(reference):
            error = convergence.calculate_all_max_errors(
                [simulation.filename], reference,
                unit=benchmark['unit_IDs'][0], which='outlet')[0]
            error = f"{error:.4e}"
        else:
            error = f"reference {Path(reference).name} missing"

        print(f"{name}: {steps:.0f} IDAS steps, {seconds:.1f} s, "
              f"outlet Linf error {error}", flush=True)


if __name__ == '__main__':
    main()

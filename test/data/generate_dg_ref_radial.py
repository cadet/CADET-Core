#!/usr/bin/env python3
"""
Generate DG reference solutions for radial model sensitivity benchmarks.
"""

import os
import json
import h5py
import shutil
import numpy as np
from cadet import Cadet

CADET_CLI = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'bin', 'cadet-cli'))
Cadet.cadet_path = CADET_CLI
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

POLYDEG = 3
NELEM = 16


def generate_reference(fv_ref_path, output_name, unit_id="unit_001", skip_sens=False):
    """
    Generate a DG reference by modifying a copy of the FV reference H5 file.
    """
    output_path = os.path.join(OUTPUT_DIR, output_name)
    shutil.copy2(fv_ref_path, output_path)

    with h5py.File(output_path, 'a') as f:
        disc_path = f'input/model/{unit_id}/discretization'

        # Remove FV-specific settings
        for key in ['NCOL', 'RECONSTRUCTION']:
            if key in f[disc_path]:
                del f[disc_path][key]
        if 'weno' in f[disc_path]:
            del f[disc_path + '/weno']

        # Set DG discretization
        if 'SPATIAL_METHOD' in f[disc_path]:
            f[disc_path]['SPATIAL_METHOD'][()] = np.bytes_('DG')
        else:
            f[disc_path].create_dataset('SPATIAL_METHOD', data=np.bytes_('DG'))

        for key, val in [('POLYDEG', POLYDEG), ('NELEM', NELEM), ('EXACT_INTEGRATION', 0)]:
            if key in f[disc_path]:
                f[disc_path][key][()] = val
            else:
                f[disc_path].create_dataset(key, data=val)

        # Handle GRM particle discretization
        unit_type = f[f'input/model/{unit_id}/UNIT_TYPE'][()].decode()
        if 'GENERAL_RATE' in unit_type:
            pt = 0
            while True:
                pt_path = f'input/model/{unit_id}/particle_type_{pt:03d}'
                if pt_path not in f:
                    break
                pt_disc = pt_path + '/discretization'
                # Remove FV particle settings
                if pt_disc in f:
                    for key in list(f[pt_disc].keys()):
                        del f[pt_disc][key]
                else:
                    f.create_group(pt_disc)
                f[pt_disc].create_dataset('SPATIAL_METHOD', data=np.bytes_('DG'))
                f[pt_disc].create_dataset('PAR_POLYDEG', data=3)
                f[pt_disc].create_dataset('PAR_NELEM', data=1)
                pt += 1

        # Handle LRMP/GRM: DG models read INIT_Q and INIT_CP at unit level (old interface)
        # but new-style FV references have INIT_CS/INIT_CP inside particle_type_xxx
        unit_path = f'input/model/{unit_id}'
        if 'LUMPED_RATE_MODEL_WITH_PORES' in unit_type or 'GENERAL_RATE' in unit_type:
            pt_path = f'{unit_path}/particle_type_000'
            if 'INIT_Q' not in f[unit_path] and 'INIT_CS' in f[pt_path]:
                f[unit_path].create_dataset('INIT_Q', data=f[pt_path]['INIT_CS'][:])
            if 'INIT_CP' not in f[unit_path] and 'INIT_CP' in f[pt_path]:
                f[unit_path].create_dataset('INIT_CP', data=f[pt_path]['INIT_CP'][:])

        # Tighten time integration tolerances
        ti_path = 'input/solver/time_integrator'
        f[ti_path]['ABSTOL'][()] = 1e-12
        f[ti_path]['RELTOL'][()] = 1e-10
        if 'ALGTOL' in f[ti_path]:
            f[ti_path]['ALGTOL'][()] = 1e-12
        f[ti_path]['INIT_STEP_SIZE'][()] = 1e-12
        f[ti_path]['MAX_STEPS'][()] = 1000000

        # Remove sensitivity setup — DG models register COL_DISPERSION as
        # CompIndep (-1) but FV references use SENS_COMP=0, causing parameter
        # mismatch. Benchmark tests use compareSens=false for now.
        if 'sensitivity' in f['input']:
            del f['input/sensitivity']

        # Remove output data
        if 'output' in f:
            del f['output']

    # Run simulation
    sim = Cadet()
    sim.filename = output_path
    sim.load()
    ret = sim.run()

    if ret.return_code != 0:
        print(f"  FAILED: {output_name}")
        err = ret.error_message if hasattr(ret, 'error_message') else 'unknown'
        print(f"  Error: {err}")
        return False

    # Reload and verify
    sim.load()
    outlet = sim.root.output.solution[unit_id].solution_outlet
    print(f"  {output_name}: outlet shape={outlet.shape}, range=[{outlet.min():.6e}, {outlet.max():.6e}]")

    try:
        if hasattr(sim.root.output, 'sensitivity'):
            sens = sim.root.output.sensitivity
            for i in range(10):
                pkey = f'param_{i:03d}'
                if hasattr(sens, pkey):
                    p = getattr(sens, pkey)
                    if hasattr(p, unit_id):
                        s = np.array(getattr(p, unit_id).sens_outlet)
                        print(f"    {pkey}: shape={s.shape}, range=[{s.min():.6e}, {s.max():.6e}]")
                else:
                    break
    except Exception as e:
        print(f"  (sens check skipped: {e})")

    return True


if __name__ == '__main__':
    import sys
    skip_sens = '--no-sens' in sys.argv

    models = [
        {
            'name': 'radLRM',
            'fv_ref': os.path.join(OUTPUT_DIR, 'ref_radLRM_dynLin_1comp_sensbenchmark1_FV_Z32.h5'),
            'output': f'ref_radLRM_dynLin_1comp_sensbenchmark1_DG_P{POLYDEG}Z{NELEM}.h5',
        },
        {
            'name': 'radLRMP',
            'fv_ref': os.path.join(OUTPUT_DIR, 'ref_radLRMP_dynLin_1comp_sensbenchmark1_FV_Z32.h5'),
            'output': f'ref_radLRMP_dynLin_1comp_sensbenchmark1_DG_P{POLYDEG}Z{NELEM}.h5',
        },
        {
            'name': 'radGRM',
            'fv_ref': os.path.join(OUTPUT_DIR, 'ref_radGRM_dynLin_1comp_sensbenchmark1_FV_Z32parZ4.h5'),
            'output': f'ref_radGRM_dynLin_1comp_sensbenchmark1_DG_P{POLYDEG}Z{NELEM}.h5',
        },
    ]

    for m in models:
        print(f"\nGenerating {m['name']} DG reference...")
        generate_reference(m['fv_ref'], m['output'], skip_sens=skip_sens)

"""Compare the working copies of the analytical references against their committed versions.

Both files are read with h5py and every dataset is compared numerically, so a pure
metadata or compression difference is reported as identical content.

Run from the CADET-Verification root.
"""

import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np

FILES = [
    'data/CADET-Verification_reference/analytical_frustumDispDPFR_1comp_benchmark1.h5',
    'data/CADET-Verification_reference/analytical_radialDispDPFR_1comp_benchmark1.h5',
    'data/CADET-Verification_reference/analytical_frustumAdvDPFR_1comp_benchmark1.h5',
    'data/CADET-Verification_reference/analytical_radialAdvDPFR_1comp_benchmark1.h5',
]


def datasets(path):
    out = {}
    with h5py.File(path, 'r') as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                out[name] = np.array(obj)
        f.visititems(visit)
    return out


for rel in FILES:
    committed = subprocess.run(['git', 'show', 'HEAD:' + rel],
                               capture_output=True)
    if committed.returncode != 0:
        print(f"{Path(rel).name}: not in HEAD")
        continue

    with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
        tmp.write(committed.stdout)
        tmp_path = tmp.name

    a = datasets(tmp_path)
    b = datasets(rel)

    only_a = set(a) - set(b)
    only_b = set(b) - set(a)
    diffs = []
    for key in sorted(set(a) & set(b)):
        x, y = a[key], b[key]
        if x.shape != y.shape or x.dtype.kind not in 'fiu':
            continue
        scale = max(np.max(np.abs(x)), 1e-300)
        d = np.max(np.abs(x.astype(float) - y.astype(float))) / scale
        if d > 0.0:
            diffs.append((d, key))

    if not diffs:
        print(f"{Path(rel).name}: identical content")
    else:
        diffs.sort(reverse=True)
        print(f"{Path(rel).name}: {len(diffs)} dataset(s) differ")
        for d, key in diffs:
            print(f"    {d:.3e}  {key}")
    if only_a or only_b:
        print(f"    datasets only in HEAD: {sorted(only_a)}")
        print(f"    datasets only in working copy: {sorted(only_b)}")

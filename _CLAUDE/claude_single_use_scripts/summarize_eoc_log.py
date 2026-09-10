"""Print a compact EOC summary of a run_convergence_analysis log.

Usage: python summarize_eoc_log.py <logfile> [<logfile> ...]
"""

import ast
import io
import re
import sys


def summarize(path):
    text = io.open(path, encoding='utf-8-sig', errors='replace').read()
    lines = text.splitlines()

    setting = None
    method = None
    for line in lines:
        stripped = line.strip()
        if stripped.endswith('_benchmark1') and not stripped.startswith('{'):
            setting = stripped
        m = re.match(r'^Method:\s+(\d+)$', stripped)
        if m:
            method = int(m.group(1))
        if stripped.startswith('{0: {'):
            table = ast.literal_eval(stripped)[0]
            n_elem = table['$N_e^z$']
            l2 = table['$L^2$ EOC'][1:]
            print(f"{setting:45s} P{method}  "
                  f"Ne {n_elem[0]}..{n_elem[-1]}  "
                  f"L2 EOC " + " ".join(f"{v:5.2f}" for v in l2))


for path in sys.argv[1:]:
    print(f"===== {path} =====")
    summarize(path)

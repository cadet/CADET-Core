#!/usr/bin/env python3
"""
Generate CADET-Core FV reference (WENO3, 10000 cells) for:
  CADET-Julia/test/test12_radial_flow/test12_EOC_vs_FV_porosity1.jl

Physical parameters match the Julia test:
  rin=0.1, rout=1.1, VELOCITY_COEFF=0.01, COL_DISPERSION=1e-4, porosity=1.0
  Gaussian inlet: c(t) = exp(-((t - 50) / 10)^2)

Accuracy design (1e-12 target):
  Spatial  : WENO3 at NCOL=10000 -> h=1e-4, error ~ h^3 = 1e-12
  Inlet    : Hermite cubic sections, dt=0.02 s in [0, 100 s]
             -> error <= dt^4 / 384 * 1.2e-3 = 5e-13  (< 1e-12)
             Coarse sections of 10 s in [100, 500 s]
             -> G(t) < e^-25 ~ 1e-11 there; contribution negligible
  Temporal : ABSTOL = 1e-12

CADET note on section naming:
  PiecewiseCubicPoly.cpp uses std::setw(3) (minimum 3 digits),
  so sec_000 .. sec_999 then sec_1000 .. sec_N.  Python f'{i:03d}'
  matches this automatically.

Output CSV (overwrites existing file; same format):
  columns: index, Time_s, C_inlet_mol_m3, C_quarter_mol_m3,
           C_mid_mol_m3, C_outlet_mol_m3
  (columns 3-5 are 0.0; Julia test only reads column 6 = outlet)
  5001 rows: t = 0.0 to 500.0 s at dt = 0.1 s

Usage:
    /Users/yuvj/CADET/build/venv/bin/python3 \\
        /Users/yuvj/CADET/test/data/generate_fv_ref_gaussianPulse_porosity1.py
"""

import csv
import tempfile
import numpy as np
from pathlib import Path
from addict import Dict
from cadet import Cadet

# =========================================================================
# PATHS
# =========================================================================

CADET_PATH = "/Users/yuvj/CADET/build/src/cadet-cli/cadet-cli"
OUT_CSV = Path(
    "/Users/yuvj/CADET-Julia/test/test12_radial_flow/"
    "fv_reference_porosity1_koren_10000cells_chromatogram.csv"
)

Cadet.cadet_path = CADET_PATH

# =========================================================================
# PHYSICAL PARAMETERS  (must match test12_EOC_vs_FV_porosity1.jl)
# =========================================================================

RIN      = 0.1       # inner radius [m]
ROUT     = 1.1       # outer radius [m]
V        = 0.01      # VELOCITY_COEFF  (u(r) = V/r)
D        = 1.0e-4    # COL_DISPERSION  [m²/s]
POROSITY = 1.0
T_CENTER = 50.0      # Gaussian centre [s]
T_WIDTH  = 10.0      # Gaussian sigma  [s]

# =========================================================================
# OUTPUT TIME GRID  (match existing CSV: 0–500 s, dt=0.1 s, 5001 points)
# =========================================================================

T_FINAL  = 500.0
DT_SAVE  = 0.1
N_SAVE   = int(round(T_FINAL / DT_SAVE)) + 1   # 5001
t_save   = list(np.linspace(0.0, T_FINAL, N_SAVE))

# =========================================================================
# INLET SECTION GRID  (hybrid: fine in pulse region, coarse elsewhere)
#
#   Hermite cubic error bound: dt^4 / 384 * max|G''''|
#   max|G''''| = 12 / T_WIDTH^4 = 12 / 10000 = 1.2e-3  (at t = T_CENTER)
#
#   dt = 0.02 s -> error <= (0.02)^4 / 384 * 1.2e-3 = 5e-13  (< 1e-12)
#   dt = 10  s -> acceptable in [100, 500] where G < exp(-25) ~ 1e-11
# =========================================================================

# Fine region: [0, T_PULSE_END] s at dt = DT_FINE
#   G = 1e-12 when |t - T_CENTER| = T_WIDTH * sqrt(12 * ln(10)) ~ 52.5 s
#   -> G < 1e-12 for t > 102.5 s; use 100 s (G(100) = e^-25 ~ 1.4e-11)
#   chosen so that (T_FINAL - T_PULSE_END) is exactly divisible by DT_COARSE
T_PULSE_END = 100.0
DT_FINE     = 0.02
N_FINE      = int(round(T_PULSE_END / DT_FINE))   # 5200 sections

# Coarse region: (T_PULSE_END, T_FINAL] s at dt = DT_COARSE
DT_COARSE = 10.0
N_COARSE  = int(round((T_FINAL - T_PULSE_END) / DT_COARSE))  # 40 sections

# Build section boundary list
sec_times_fine   = [i * DT_FINE for i in range(N_FINE + 1)]          # 5201 pts
coarse_start     = sec_times_fine[-1]                                  # = 104.0
sec_times_coarse = [coarse_start + j * DT_COARSE for j in range(1, N_COARSE + 1)]
sec_times        = sec_times_fine + sec_times_coarse  # 5241 boundaries

N_SEC = len(sec_times) - 1  # 5240 sections total
assert abs(sec_times[-1] - T_FINAL) < 1e-9, "section grid must reach T_FINAL"

print(f"Inlet sections: {N_FINE} fine (dt={DT_FINE}s) + {N_COARSE} coarse (dt={DT_COARSE}s)"
      f" = {N_SEC} total")
print(f"Fine inlet error bound: {(DT_FINE)**4 / 384 * 1.2e-3:.1e}  (< 1e-12)")

# =========================================================================
# GAUSSIAN HELPERS
# =========================================================================

def gaussian(t):
    return float(np.exp(-((t - T_CENTER) / T_WIDTH) ** 2))

def gaussian_deriv(t):
    """dG/dt"""
    return float(-2.0 * (t - T_CENTER) / T_WIDTH**2 * gaussian(t))

def hermite_coeffs(t_i, h):
    """CADET piecewise-cubic coefficients for the Taylor expansion
    c(tau) = a0 + a1*tau + a2*tau^2 + a3*tau^3  (tau = t - t_i)
    that Hermite-interpolates G on [t_i, t_i + h]."""
    y0, y1 = gaussian(t_i), gaussian(t_i + h)
    d0, d1 = gaussian_deriv(t_i), gaussian_deriv(t_i + h)
    a0 = y0
    a1 = d0
    a2 = (3.0 * (y1 - y0) / h - 2.0 * d0 - d1) / h
    a3 = (2.0 * (y0 - y1) / h + d0 + d1) / h**2
    return a0, a1, a2, a3

# =========================================================================
# BUILD CADET MODEL
# =========================================================================

sim = Dict()

# --- network ---
sim.input.model.NUNITS = 3
sim.input.model.connections.NSWITCHES = 1
sim.input.model.connections.switch_000.CONNECTIONS = [
    0.0, 1.0, -1.0, -1.0, 1.0,
    1.0, 2.0, -1.0, -1.0, 1.0,
]
sim.input.model.connections.switch_000.SECTION = 0
sim.input.model.solver.GS_TYPE      = 1
sim.input.model.solver.MAX_KRYLOV   = 0
sim.input.model.solver.MAX_RESTARTS = 10
sim.input.model.solver.SCHUR_SAFETY = 1e-8

# --- unit 000: inlet with Gaussian pulse (Hermite cubic) ---
sim.input.model.unit_000.UNIT_TYPE  = 'INLET'
sim.input.model.unit_000.INLET_TYPE = 'PIECEWISE_CUBIC_POLY'
sim.input.model.unit_000.NCOMP      = 1

print(f"Building {N_SEC} inlet sections ...", flush=True)
for i in range(N_SEC):
    t_i = sec_times[i]
    h_i = sec_times[i + 1] - t_i
    a0, a1, a2, a3 = hermite_coeffs(t_i, h_i)
    sec = sim.input.model.unit_000[f'sec_{i:03d}']
    sec.CONST_COEFF = [a0]
    sec.LIN_COEFF   = [a1]
    sec.QUAD_COEFF  = [a2]
    sec.CUBE_COEFF  = [a3]
print("  done.")

# --- unit 001: radial LRM without pores ---
u1 = sim.input.model.unit_001
u1.UNIT_TYPE          = 'RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES'
u1.NCOMP              = 1
u1.NPARTYPE           = 1
u1.TOTAL_POROSITY     = POROSITY
u1.VELOCITY_COEFF     = V
u1.COL_DISPERSION     = D
u1.COL_RADIUS_INNER   = RIN
u1.COL_RADIUS_OUTER   = ROUT
u1.INIT_C             = [0.0]
u1.particle_type_000.HAS_FILM_DIFFUSION    = False
u1.particle_type_000.HAS_PORE_DIFFUSION    = False
u1.particle_type_000.HAS_SURFACE_DIFFUSION = False
u1.particle_type_000.NBOUND                = [0]

disc = u1.discretization
disc.SPATIAL_METHOD       = 'FV'
disc.NCOL                 = 100000
disc.USE_ANALYTIC_JACOBIAN = 1
disc.RECONSTRUCTION       = 'WENO'
disc.weno.WENO_ORDER      = 3
disc.weno.WENO_EPS        = 1e-10
disc.weno.BOUNDARY_MODEL  = 1   # ZeroWeights: 2nd-order at boundary vs 1st for BOUNDARY_MODEL=0
disc.koren.KOREN_EPS      = 1e-10

# --- unit 002: outlet ---
sim.input.model.unit_002.UNIT_TYPE = 'OUTLET'
sim.input.model.unit_002.NCOMP     = 1

# --- solver ---
sim.input.solver.sections.NSEC               = N_SEC
sim.input.solver.sections.SECTION_CONTINUITY = [1] * (N_SEC - 1)
sim.input.solver.sections.SECTION_TIMES      = sec_times
sim.input.solver.CONSISTENT_INIT_MODE        = 1
sim.input.solver.NTHREADS                    = 1
sim.input.solver.USER_SOLUTION_TIMES         = t_save

ti = sim.input.solver.time_integrator
ti.ABSTOL         = 1e-12
ti.RELTOL         = 1e-10
ti.ALGTOL         = 1e-12
ti.INIT_STEP_SIZE = 1e-6
ti.MAX_STEPS      = 10_000_000

# --- return (outlet only; SOLUTION_BULK not needed) ---
sim.input['return'].SPLIT_COMPONENTS_DATA            = 0
sim.input['return'].SPLIT_PORTS_DATA                 = 0
sim.input['return'].unit_001.WRITE_SOLUTION_OUTLET   = True
sim.input['return'].unit_001.WRITE_SOLUTION_INLET    = False
sim.input['return'].unit_001.WRITE_SOLUTION_BULK     = False

# =========================================================================
# RUN CADET
# =========================================================================

cadet_sim = Cadet()
cadet_sim.root = sim

with tempfile.NamedTemporaryFile(suffix='.h5', delete=False) as tmp:
    tmp_path = tmp.name

cadet_sim.filename = tmp_path
print("Writing HDF5 input ...", flush=True)
cadet_sim.save()

print(f"Running CADET FV (WENO3, NCOL={disc.NCOL}, {N_SEC} sections) ...", flush=True)
ret = cadet_sim.run_simulation()
if ret.return_code != 0:
    Path(tmp_path).unlink(missing_ok=True)
    raise RuntimeError(f"Simulation failed:\n{ret.error_message}")

cadet_sim.load_from_file()

# =========================================================================
# EXTRACT & WRITE CSV
# =========================================================================

outlet = np.array(cadet_sim.root.output.solution.unit_001.solution_outlet).ravel()
times  = np.array(t_save)

assert len(outlet) == N_SAVE, f"Expected {N_SAVE} time points, got {len(outlet)}"

i_peak = int(np.argmax(outlet))
print(f"Peak outlet: {outlet[i_peak]:.8e} at t = {times[i_peak]:.2f} s")

with open(OUT_CSV, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['', 'Time_s', 'C_inlet_mol_m3',
                     'C_quarter_mol_m3', 'C_mid_mol_m3', 'C_outlet_mol_m3'])
    for i, (t, c) in enumerate(zip(times, outlet)):
        writer.writerow([
            i,
            f"{t:.16e}",
            f"{0.0:.16e}",
            f"{0.0:.16e}",
            f"{0.0:.16e}",
            f"{c:.16e}",
        ])

Path(tmp_path).unlink(missing_ok=True)
print(f"Saved: {OUT_CSV}")
# Import libraries
import numpy as np
import matplotlib.pyplot as plt

from cadet import Cadet
Cadet.cadet_path = '/path/to/cadet-core/install/prefix' # (e.g., '/usr/local')

# Create model object
model = Cadet()

# Number of unit operations
model.root.input.model.nunits = 3

# Inlet
model.root.input.model.unit_000.unit_type = 'INLET'
model.root.input.model.unit_000.ncomp = 1
model.root.input.model.unit_000.inlet_type = 'PIECEWISE_CUBIC_POLY'

# General Rate Model
model.root.input.model.unit_001.unit_type = 'GENERAL_RATE_MODEL'
model.root.input.model.unit_001.ncomp = 1

## Geometry
model.root.input.model.unit_001.col_length = 0.1                # m
model.root.input.model.unit_001.cross_section_area = 0.01       # m
model.root.input.model.unit_001.col_porosity = 0.37             # -
model.root.input.model.unit_001.par_porosity = 0.33             # -
model.root.input.model.unit_001.par_radius = 1e-6               # m
                
## Bulk transport
model.root.input.model.unit_001.col_dispersion = 1e-8           # m^2 / s (interstitial volume)

## Particle transport
model.root.input.model.unit_001.particle_type_000.film_diffusion = [1e-5]         # m / s
model.root.input.model.unit_001.particle_type_000.pore_diffusion = [1e-10,]        # m^2 / s (mobile phase)  
model.root.input.model.unit_001.particle_type_000.surface_diffusion = [0.0,]      # m^2 / s (solid phase)

## Adsorption
model.root.input.model.unit_001.particle_type_000.adsorption_model = 'MULTI_COMPONENT_LANGMUIR'
model.root.input.model.unit_001.particle_type_000.nbound = [1]
model.root.input.model.unit_001.particle_type_000.adsorption.is_kinetic = True    # Kinetic binding
model.root.input.model.unit_001.particle_type_000.adsorption.mcl_ka = [1.0,]      # m^3 / (mol * s)   (mobile phase)
model.root.input.model.unit_001.particle_type_000.adsorption.mcl_kd = [1.0,]      # 1 / s (desorption)
model.root.input.model.unit_001.particle_type_000.adsorption.mcl_qmax = [100.0,]  # mol / m^3   (solid phase)

## Initial conditions
model.root.input.model.unit_001.init_c = [0.0,]
model.root.input.model.unit_001.particle_type_000.init_cp = [0.0,]
model.root.input.model.unit_001.particle_type_000.init_cs = [0.0,]

## Discretization
### Grid cells
model.root.input.model.unit_001.discretization.spatial_method = "DG" # "FV" recommended for strong gradients, where DG might oscillate
model.root.input.model.unit_001.discretization.polydeg = 3 # recommended to be increased up to 5 if higher resolution is desired for smooth solutions
model.root.input.model.unit_001.discretization.nelem = 3 # recommended to be increased if higher resolution is desired
model.root.input.model.unit_001.particle_type_000.discretization.spatial_method = "DG"
model.root.input.model.unit_001.particle_type_000.discretization.par_polydeg = 3 # recommended to be increased if higher resolution is desired
model.root.input.model.unit_001.particle_type_000.discretization.par_nelem = 1 #  recommended to stay fixed at 1 except for targeted resolution of particle (ie user defined par_disc_type)
### Other options
model.root.input.model.unit_001.discretization.particle_type_000.par_disc_type = 'EQUIDISTANT'    
model.root.input.model.unit_001.discretization.use_analytic_jacobian = 1

## Outlet
model.root.input.model.unit_002.unit_type = 'OUTLET'
model.root.input.model.unit_002.ncomp = 1

# Sections 
model.root.input.solver.sections.nsec = 1
model.root.input.solver.sections.section_times = [0.0, 1200,]   # s
model.root.input.solver.sections.section_continuity = []

# Inlet sections
model.root.input.model.unit_000.sec_000.const_coeff = [1.0e-3,] # mol / m^3
model.root.input.model.unit_000.sec_000.lin_coeff = [0.0,]
model.root.input.model.unit_000.sec_000.quad_coeff = [0.0,]
model.root.input.model.unit_000.sec_000.cube_coeff = [0.0,]

# Switches
model.root.input.model.connections.nswitches = 1
model.root.input.model.connections.switch_000.section = 0
model.root.input.model.connections.switch_000.connections = [
    0, 1, -1, -1, 60/1e6,  # [unit_000, unit_001, all components, all components, Q/ m^3*s^-1 
    1, 2, -1, -1, 60/1e6]  # [unit_001, unit_002, all components, all components, Q/ m^3*s^-1 

# Solver settings
model.root.input.model.solver.gs_type = 1
model.root.input.model.solver.max_krylov = 0
model.root.input.model.solver.max_restarts = 10
model.root.input.model.solver.schur_safety = 1e-8

# Number of cores for parallel simulation
model.root.input.solver.nthreads = 1

# Tolerances for the time integrator
model.root.input.solver.time_integrator.abstol = 1e-6
model.root.input.solver.time_integrator.algtol = 1e-10
model.root.input.solver.time_integrator.reltol = 1e-6
model.root.input.solver.time_integrator.init_step_size = 1e-6
model.root.input.solver.time_integrator.max_steps = 1000000

# Return data
model.root.input['return'].split_components_data = 0
model.root.input['return'].split_ports_data = 0
model.root.input['return'].unit_000.write_solution_bulk = 1
model.root.input['return'].unit_000.write_solution_inlet = 1
model.root.input['return'].unit_000.write_solution_outlet = 1

# Copy settings to the other unit operations
model.root.input['return'].unit_001 = model.root.input['return'].unit_000
model.root.input['return'].unit_002 = model.root.input['return'].unit_000

# Solution times
model.root.input.solver.user_solution_times = np.linspace(0, 1200, 1001)

# Save and run simulation
model.filename = 'model.h5'
model.save()

data = model.run()

if data.return_code == 0:
    print("Simulation completed successfully")
    model.load()   
else:
    print(data)
    raise Exception("Simulation failed")

# Plot results
plt.figure()

time = model.root.output.solution.solution_times
c = model.root.output.solution.unit_001.solution_outlet
plt.plot(time/60, c)
plt.xlabel('$time~/~min$')
plt.ylabel('$Outlet~concentration~/~mol \cdot m^{-3} $')
plt.show()

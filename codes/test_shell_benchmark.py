# -*- coding: utf-8 -*-
from parameters import ParameterHandler
from gravity_field import GravityType
from time_stepping import IMEXType
from buoyant_fluid_solver import SolverType
#------------------------------------------------------------------------------#
params = ParameterHandler()
# geometry parameters
params.dim = 2
params.radii = (0.35, 1.0)
# runtime parameters
params.t_end = 150.
params.n_steps = 20000
params.output_frequency = 25
params.checkpoint_frequency = 100
params.global_avg_frequency = 10
params.cfl_frequency = 10
# solver parameters
params.use_assembler_method = True
params.solver_type = SolverType.implicit_solver
# time stepping parameters
params.adaptive_time_stepping = True
params.timestep = 5e-5
params.min_cfl = 0.4
params.max_cfl = 0.6
params.max_timestep = 1e-4
params.min_timestep = 1e-9
# dimensionless constants
params.rotation = True
params.buoyancy = True
params.gravity_type = GravityType.radial_linear
params.ekman = 1e-3
params.rayleigh = 1e5
params.prandtl = 1.0
#------------------------------------------------------------------------------#
# auxiliary ids
inner_circle_id = 1
outer_circle_id = 2
#------------------------------------------------------------------------------#
# mesh generation
from grid_generator import spherical_shell
boundary_ids = (inner_circle_id, outer_circle_id)
mesh, boundary_marker = spherical_shell(params.dim, params.radii, boundary_ids)
#------------------------------------------------------------------------------#
# mesh generation
import boundary_conditions
bcs = boundary_conditions.Christensen2001Case0(boundary_ids)
#------------------------------------------------------------------------------#
import initial_conditions
ics = initial_conditions.Christensen2001Case0(2, params.radii)
#------------------------------------------------------------------------------#
from buoyant_fluid_solver import BuoyantFluidSolver
solver = BuoyantFluidSolver(mesh, boundary_marker, bcs, params, ics)
solver.run()
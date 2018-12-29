# -*- coding: utf-8 -*-
from parameters import ParameterHandler, GravityType
#------------------------------------------------------------------------------#
params = ParameterHandler()
# runtime parameters
params.t_end = 150.
params.n_steps = 50000
params.output_frequency = 50
params.checkpoint_frequency = 1000
params.rms_frequency = 5
# solver parameters
params.use_assembler_method = True
# time stepping parameters
params.adaptive_time_stepping = True
params.timestep = 1e-6
params.min_cfl = 0.4
params.max_cfl = 0.6
params.max_timestep = 5e-1
params.min_timestep = 1e-12
# dimensionless constants
params.rotation = True
params.buoyancy = True
params.gravity_type = GravityType.radial
params.ekman = 1e-3
params.rayleigh = 1e5
params.prandtl = 1.0
#------------------------------------------------------------------------------#
# model parameters
inner_temperature = 0.5
outer_temperature = -0.5
# geometry parameters
aspect_ratio = 0.35
outer_radius = 1. 
inner_radius = aspect_ratio * outer_radius
n_initial_refinements = 1
# auxiliary ids
inner_circle_id = 1
outer_circle_id = 2
#------------------------------------------------------------------------------#
# mesh creation
#
# TODO: I think this and the boundary ids can go to a grid generator function
#
# dolfin mshr module
from dolfin import Point
from mshr import Circle, generate_mesh
center = Point(0., 0.)
domain = Circle(center, outer_radius) \
       - Circle(center, inner_radius)
mesh = generate_mesh(domain, 50)
for i in range(n_initial_refinements):
    from dolfin import refine
    mesh = refine(mesh)
#------------------------------------------------------------------------------#
# subdomains for boundaries
from dolfin import MeshFunctionSizet, SubDomain
facet_marker = MeshFunctionSizet(mesh, mesh.topology().dim() - 1)
facet_marker.set_all(0)
# size of smallest element
hmin = mesh.hmin()
# inner circle boundary
class InnerCircle(SubDomain):
    def inside(self, x, on_boundary):
        # tolerance: half length of smallest element
        tol = hmin / 2. 
        from math import sqrt
        result = abs(sqrt(x[0]**2 + x[1]**2) - inner_radius) < tol
        return result and on_boundary
# outer cirlce boundary
class OuterCircle(SubDomain):
    def inside(self, x, on_boundary):
        # tolerance: half length of smallest element
        tol = hmin/2.
        from math import sqrt
        result = abs(sqrt(x[0]**2 + x[1]**2) - outer_radius) < tol
        return result and on_boundary
# mark boundaries
gamma_inner = InnerCircle()
gamma_inner.mark(facet_marker, inner_circle_id)
gamma_outer = OuterCircle()
gamma_outer.mark(facet_marker, outer_circle_id)
#------------------------------------------------------------------------------#
# initial condition for temperature
#
# TODO: I think this could go to an initial condition class
#
if mesh.topology().dim() == 2:
    radius_string = "sqrt( pow(x[0],2) + pow(x[1],2))"
else:
    radius_string = "sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2))"
from dolfin import Expression
ics = dict()
ics["temperature"] = Expression(
            "(To - Ti) / (ro - ri) * (sqrt( pow(x[0],2) + pow(x[1],2)) -ri) + Ti ",
            ri = inner_radius, ro = outer_radius,
            Ti = inner_temperature, To = outer_temperature,
            degree = 2)
#------------------------------------------------------------------------------#
from buoyant_fluid_solver import VelocityBCType, TemperatureBCType
bcs = dict()
bcs["velocity"] = ((VelocityBCType.no_slip , inner_circle_id, None), 
                   (VelocityBCType.no_slip , outer_circle_id, None))
bcs["temperature"] = ((TemperatureBCType.constant, inner_circle_id, inner_temperature), 
                      (TemperatureBCType.constant, outer_circle_id, outer_temperature))
#----------------------------------------------------------------------------#
from buoyant_fluid_solver import BuoyantFluidSolver
solver = BuoyantFluidSolver(mesh, facet_marker, bcs, params)
solver.run()
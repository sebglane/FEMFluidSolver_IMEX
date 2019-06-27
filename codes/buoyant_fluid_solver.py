#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import dolfin as dlfn
from time_stepping import IMEXCoefficients

from enum import Enum
class SolverType(Enum):
    imex_solver = 0
    implicit_solver = 1

# define auxiliary operators
def a_op(phi, psi):
    from dolfin import inner, grad
    return inner(grad(phi), grad(psi))
def b_op(phi, psi):
    from dolfin import div
    return div(phi) * psi
# 2D cross product
def cross_2D(phi):
    from dolfin import as_vector
    return as_vector([phi[1], -phi[0]])

class BuoyantFluidSolver:
    def __init__(self, mesh, facet_ids, bcs, params, ics = dict()):
        # input check: mesh
        assert isinstance(mesh, dlfn.Mesh)
        self._mesh = mesh
        self._space_dim = self._mesh.topology().dim()
        assert self._space_dim in (2,3)
        # input check: facet_ids
        assert isinstance(facet_ids, dlfn.MeshFunctionSizet)
        self._facet_markers = facet_ids
        # input check: boundary conditions
        self._check_boundary_conditions(bcs)
        self._bcs = bcs
        # input check: parameters
        from parameters import ParameterHandler
        assert isinstance(params, ParameterHandler)
        self._parameters = params
        # input check: initial conditions
        assert isinstance(ics, dict)
        self._ics = ics
        # equation coefficients
        tmp = self._parameters.coefficients()
        self._coefficients = tuple([dlfn.Constant(t) for t in tmp])
        # imex coefficients
        self._imex = IMEXCoefficients(self._parameters.imex_type)        
        # gravity field
        import gravity_field as gf
        if params.gravity_type not in (gf.GravityType.radial, gf.GravityType.radial_linear):
            self._gravity = gf.get_gravity_field(self._space_dim, params.gravity_type)
        else:
            self._gravity = gf.get_gravity_field(self._space_dim, params.gravity_type, radii=params.radii)
        # initialize timestep as dolfin constant
        self._timestep = dlfn.Constant(1.0)
        self._timestep.assign(self._parameters.timestep)
        # runtime flags
        self._rebuild_matrices = True
        # helpful constants
        self._one = dlfn.Constant(1.0)
        self._omega = dlfn.Constant(1.0)
        if self._space_dim == 2:
            self._null_vector = dlfn.Constant((0.0, 0.0))
        elif self._space_dim == 3:
            self._null_vector = dlfn.Constant((0.0, 0.0, 0.0))
        else:
            raise ValueError()
        
        print "The system's version of FEniCS is ", dlfn.dolfin_version(), "."
        
    def _check_boundary_conditions(self, bcs):
        from boundary_conditions import VelocityBCType, TemperatureBCType
        # input check: boundary conditions
        assert isinstance(bcs, dict)
        assert bcs.has_key("temperature")
        assert bcs.has_key("velocity")
        # check if structure of dictionary is correct
        ids_bc = set()
        for key, bc in bcs.iteritems():
            if key is "velocity":
                bc_types = VelocityBCType
                none_type = VelocityBCType.no_slip
            elif key is "temperature":
                bc_types = TemperatureBCType
                none_type = None
            else:
                raise ValueError()
            const_type = bc_types.constant
            assert isinstance(bc, tuple)
            assert len(bc) > 0
            # check sub boundary conditions
            for i in range(len(bc)):
                assert isinstance(bc[i], tuple)
                assert len(bc[i]) == 3
                # check if type of bc is correct
                assert bc[i][0] in bc_types
                # check if type of facet_id is correct
                assert isinstance(bc[i][1], int) and bc[i][1] > 0
                ids_bc.add(bc[i][1])
                # check if value type of bc is correct
                if none_type is VelocityBCType.no_slip:
                    assert bc[i][2] is None
                elif bc[i][2] is const_type:
                    if key is "velocity":
                        assert isinstance(bc[i][2], (list, tuple))
                        assert len(bc[i][2]) == self._space_dim
                        assert all(isinstance(x, float) for x in bc[i][2])
                    elif key is "temperature":
                        assert isinstance(bc[i][2], float)
                else:
                    isinstance(bc[i][2], dlfn.Expression)
        # check if facet_ids of bcs occur in markers
        ids_bc = tuple(ids_bc)
        ids_bc_found = [False, ] * len(ids_bc)
        for facet in dlfn.facets(self._mesh):
            if facet.exterior():
                if self._facet_markers[facet] in ids_bc:
                    i = ids_bc.index(self._facet_markers[facet])
                    ids_bc_found[i] = True
                    if all(ids_bc_found):
                        break

    def run(self):
        self._print_info_message()        
        
        self._setup_function_spaces()
        
        self._setup_boundary_conditions()
        
        self._update_imex_coefficients()
        
        if self._parameters.solver_type is SolverType.imex_solver:
            self._setup_imex_problem()
        elif self._parameters.solver_type is SolverType.implicit_solver:
            self._setup_implicit_problem()
        
        self._setup_initial_conditions()
        
        velocity, pressure, temperature = self._sol0.split()
        self._output_results(velocity, temperature, pressure, 0, 0, True)
        self._output_global_averages(velocity, temperature, 0, 0, True)
        
        
        # initialize timestep
        timestep = self._parameters.timestep
        omega = 1.0
        # initialize timestepping
        from time_stepping import TimestepControl
        timestep_control = TimestepControl((self._parameters.min_cfl,
                                            self._parameters.max_cfl,),
                                           (self._parameters.min_timestep,
                                            self._parameters.max_timestep,))
        # update flag
        timestep_modified = False
        # dolfin timer 
        dlfn.tic()
        # time loop
        time = 0.0
        step = 0
        while time < self._parameters.t_end and step < self._parameters.n_steps + 1:
            print "Iteration: {0:08d}, time = {1:10.5f}, time step = {2:5.4e}"\
                  .format(step, time, timestep)
            if timestep_modified:
                self._timestep.assign(timestep)
                self._omega.assign(omega)
                self._update_imex_coefficients(step, omega)
                self._rebuild_matrices = True
            # compute solution
            self._solve(step)
            # extract solution
            velocity, pressure, temperature = self._sol.split()
            # compute cfl number
            cfl = self._compute_cfl_number(velocity, timestep)
            if step % self._parameters.cfl_frequency == 0:
                print "   current cfl number:\t{0:03.2e}".format(cfl)
            #===================================================================
                # modify time step
            if self._parameters.adaptive_time_stepping:
                timestep_modified, timestep, omega \
                    = timestep_control.adjust_time_step(cfl, timestep)
            elif step == 0:
                # update coefficients after initial step to switch to SBF2
                timestep_modified = True
            elif step == 1:
                timestep_modified = False
            # write output
            if step % self._parameters.output_frequency == 0:
                self._output_results(velocity, temperature, pressure, step, time)
            # compute rms values
            if step % self._parameters.global_avg_frequency == 0:
                self._output_global_averages(velocity, temperature, step, time)
            # write checkpoint
            if step % self._parameters.checkpoint_frequency == 0 and step != 0:
                self._write_checkpoint(step, time)
            # update solutions for next iteration
            self._sol00.assign(self._sol0)
            self._sol0.assign(self._sol)
            # update time
            time += timestep
            # increase counter
            step += 1
        print "elapsed simulation time: ", dlfn.toc(), " sec"
        
    def _solve(self, step):
        assert hasattr(self, "_solver")
        if self._parameters.use_assembler_method and isinstance(self._solver, dlfn.LUSolver):
            assert hasattr(self, "_assembler")
            assert hasattr(self, "_system_matrix")
            assert hasattr(self, "_system_rhs")
            if step == 0 or self._rebuild_matrices:
                # assemble system matrix / rhs
                if self._parameters.use_assembler_method:
                    print "   assembling system..."
                    self._assembler.assemble(self._system_matrix,
                                             self._system_rhs)
                    # re-compute LU decomposition
                    self._solver.parameters['reuse_factorization'] = False
                self._rebuild_matrices = False
            else:
                assert hasattr(self, "_assembler")
                assert hasattr(self, "_system_rhs")
                print "   assembling rhs..."
                self._assembler.assemble(self._system_rhs)
            self._solver.solve(self._sol.vector(), self._system_rhs)
            self._solver.parameters['reuse_factorization'] = True
        else:
            self._solver.solve()

    def _setup_function_spaces(self):
        assert hasattr(self, "_mesh")
        # element and function space definition
        cell = self._mesh.ufl_cell()
        # taylor-hood element
        self._elemV = dlfn.VectorElement("CG", cell, self._parameters.velocity_degree)
        self._elemP = dlfn.FiniteElement("CG", cell, self._parameters.pressure_degree)
        self._elemT = dlfn.FiniteElement("CG", cell, self._parameters.temperature_degree)
        self._mixedElem = dlfn.MixedElement([self._elemV,
                                             self._elemP,
                                             self._elemT])
        self._Wh = dlfn.FunctionSpace(self._mesh, self._mixedElem)
        # print info message
        ndofs_velocity = self._Wh.sub(0).dim()
        ndofs_pressure = self._Wh.sub(1).dim()
        ndofs_temperature = self._Wh.sub(2).dim()
        print "DOFs velocity : ", ndofs_velocity, "\n" \
                "DOFs pressure : ", ndofs_pressure, "\n" \
                "DOFs temperature : ", ndofs_temperature
        # functions
        self._sol = dlfn.Function(self._Wh)
        self._sol0 = dlfn.Function(self._Wh)
        self._sol00 = dlfn.Function(self._Wh)
        self._v0, _, self._T0 = dlfn.split(self._sol0)
        self._v00, _, self._T00 = dlfn.split(self._sol00)
    
    def _setup_boundary_conditions(self):
        assert hasattr(self, "_bcs")
        assert hasattr(self, "_Wh")
        assert hasattr(self, "_facet_markers")
        self._dirichlet_bcs = []
        # temperature part
        print "   setup temperature boundary conditions..."
        temperature_space = self._Wh.sub(2)
        temperature_bcs = self._bcs["temperature"]
        from boundary_conditions import TemperatureBCType
        for bc in temperature_bcs:
            if bc[0] is TemperatureBCType.constant:
                const_function = dlfn.Constant(bc[2])
                self._dirichlet_bcs.append(
                    dlfn.DirichletBC(temperature_space, const_function,
                                     self._facet_markers, bc[1]))
            elif bc[0] is TemperatureBCType.function:
                self._dirichlet_bcs.append(
                    dlfn.DirichletBC(temperature_space, bc[2],
                                     self._facet_markers, bc[1]))
            else: 
                raise NotImplementedError()
        # velocity part
        print "   setup velocity boundary conditions..."
        velocity_space = self._Wh.sub(0)
        velocity_bcs = self._bcs["velocity"]
        from boundary_conditions import VelocityBCType
        for bc in velocity_bcs:
            if bc[0] is VelocityBCType.no_slip:
                const_function = self._null_vector
                self._dirichlet_bcs.append(
                    dlfn.DirichletBC(velocity_space, const_function,
                                     self._facet_markers, bc[1]))
            elif bc[0] is VelocityBCType.constant:
                assert len(bc[2]) == self._space_dim
                const_function = dlfn.Constant(bc[2])
                self._dirichlet_bcs.append(
                    dlfn.DirichletBC(velocity_space, const_function,
                                     self._facet_markers, bc[1]))
            elif bc[0] is VelocityBCType.function:
                self._dirichlet_bcs.append(
                    dlfn.DirichletBC(velocity_space, bc[2],
                                     self._facet_markers, bc[1]))
            else: 
                raise NotImplementedError()                
                
    def _setup_initial_conditions(self):
        assert hasattr(self, "_ics")
        assert hasattr(self, "_Wh")
        print "   setup initial conditions..."
        # velocity part
        if self._ics.has_key("velocity"):
            initial_velocity_expr = self._ics["velocity"]
            assert isinstance(initial_velocity_expr, dlfn.Expression)
            assert initial_velocity_expr.ufl_shape[0] == self._space_dim
            initial_velocity = dlfn.interpolate(initial_velocity_expr,
                                                self._Wh.sub(0).collapse())
        else:
            initial_velocity = dlfn.interpolate(self._null_vector,
                                                self._Wh.sub(0).collapse())
        # pressure part
        if self._ics.has_key("pressure"):
            initial_pressure_expr = self._ics["pressure"]
            assert isinstance(initial_pressure_expr, dlfn.Expression)
            assert len(initial_pressure_expr.ufl_shape) == 0
            initial_pressure= dlfn.interpolate(initial_pressure_expr,
                                               self._Wh.sub(1).collapse())
        else:
            initial_pressure = dlfn.interpolate(dlfn.Expression("0.0", degree=1),
                                                self._Wh.sub(1).collapse())
        # temperature part
        if self._ics.has_key("temperature"):
            initial_temperature_expr = self._ics["temperature"]
            assert isinstance(initial_temperature_expr, dlfn.Expression)
            assert len(initial_temperature_expr.ufl_shape) == 0
            initial_temperature = dlfn.interpolate(initial_temperature_expr,
                                                   self._Wh.sub(2).collapse())
        else:
            initial_temperature = dlfn.interpolate(dlfn.Constant(0.),
                                                   self._Wh.sub(2).collapse())
        # assign initial condition
        dlfn.assign(self._sol0, [initial_velocity,
                                 initial_pressure,
                                 initial_temperature])
    
    def _setup_implicit_problem(self):
        assert hasattr(self, "_parameters")
        assert hasattr(self, "_mesh")
        assert hasattr(self, "_Wh")
        assert hasattr(self, "_coefficients")
        assert hasattr(self, "_one")
        assert hasattr(self, "_omega")
        assert hasattr(self, "_v0")
        assert hasattr(self, "_v00")
        assert hasattr(self, "_T0")
        assert hasattr(self, "_T00")
        #=======================================================================
        # retrieve imex coefficients
        a = self._imex_alpha
        #=======================================================================
        # trial and test function
        (del_v, del_p, del_T) = dlfn.TestFunctions(self._Wh)
        v, p, T = self._sol.split()
        # volume element
        dV = dlfn.Measure("dx", domain = self._mesh)
        # reference to time step
        timestep = self._timestep
        #=======================================================================
        from dolfin import dot, grad
        # 1) momentum equation
        momentum_eqn = dot((a[0] * v + a[1] * self._v0 + a[2] * self._v00)/ timestep, del_v) *  dV \
                     + dot(dot(grad(v), v), del_v) * dV \
                     + self._coefficients[1] * a_op(v, del_v) * dV \
                     - b_op(del_v, p) * dV \
                     - b_op(v, del_p) * dV
        # 2d) momentum equation: coriolis term
        if self._parameters.rotation is True:
            assert self._coefficients[0] != 0.0
            print "   adding rotation to the model..."
            # add Coriolis term
            if self._space_dim == 2:
                coriolis_term = cross_2D(v)
                momentum_eqn += self._coefficients[0] * (-v[1] * del_v[0] + v[0] * del_v[1]) * dV
            elif self._space_dim == 3:
                from dolfin import cross
                coriolis_term = cross(self._rotation_vector, v)
                momentum_eqn += self._coefficients[0] * dot(coriolis_term, del_v) * dV
        # 2e) rhs momentum equation: buoyancy term
        if self._parameters.buoyancy is True:
            assert self._coefficients[2] != 0.0
            print "   adding buoyancy to the model..."
            # buoyancy term
            momentum_eqn -= self._coefficients[2] * T * dot(self._gravity, del_v) * dV
        #=======================================================================        
        # 3) energy equation
        energy_eqn = (a[0] * T + a[1] * self._T0 + a[2] * self._T00)/ timestep * del_T * dV \
                   + dot(v, grad(T)) * del_T * dV \
                   + self._coefficients[3] * a_op(T, del_T) * dV
        #=======================================================================        
        # full problem
        self._eqn = momentum_eqn + energy_eqn
        if not hasattr(self, "_dirichlet_bcs"):
            self._setup_boundary_conditons()
        self._jacobian = dlfn.derivative(self._eqn, self._sol,
                                         du = dlfn.TrialFunction(self._Wh))
        problem = dlfn.NonlinearVariationalProblem(self._eqn, self._sol,
                                                   J = self._jacobian,
                                                   bcs=self._dirichlet_bcs)
        self._solver = dlfn.NonlinearVariationalSolver(problem)
    
    def _setup_imex_problem(self):
        assert hasattr(self, "_parameters")
        assert hasattr(self, "_mesh")
        assert hasattr(self, "_Wh")
        assert hasattr(self, "_coefficients")
        assert hasattr(self, "_one")
        assert hasattr(self, "_omega")
        assert hasattr(self, "_v0")
        assert hasattr(self, "_v00")
        assert hasattr(self, "_T0")
        assert hasattr(self, "_T00")
        #=======================================================================
        # retrieve imex coefficients
        a, b, c = self._imex_alpha, self._imex_beta, self._imex_gamma
        #=======================================================================
        # trial and test function
        (del_v, del_p, del_T) = dlfn.TestFunctions(self._Wh)
        (dv, dp, dT) = dlfn.TrialFunctions(self._Wh)
        # volume element
        dV = dlfn.Measure("dx", domain = self._mesh)
        # reference to time step
        timestep = self._timestep
        #=======================================================================
        from dolfin import dot, grad, inner
        # 1) lhs momentum equation
        lhs_momentum = a[0] / timestep * dot(dv, del_v) *  dV \
                     + c[0] * self._coefficients[1] * a_op(dv, del_v) * dV\
                     - b_op(del_v, dp) * dV\
                     - b_op(dv, del_p) * dV
        # 2a) rhs momentum equation: time derivative
        rhs_momentum = - dot(a[1] / timestep * self._v0 + a[2] / timestep * self._v00, del_v) * dV
        # 2b) rhs momentum equation: nonlinear term
        nonlinear_term_velocity =  b[0] * dot(grad(self._v0), self._v0) \
                                 + b[1] * dot(grad(self._v00), self._v00)
        rhs_momentum -= dot(nonlinear_term_velocity, del_v) * dV
        # 2c) rhs momentum equation: linear term
        rhs_momentum -= self._coefficients[1] * inner(c[1] * grad(self._v0) + c[2] * grad(self._v00), grad(del_v)) * dV
        # 2d) rhs momentum equation: coriolis term
        if self._parameters.rotation is True:
            assert self._coefficients[0] != 0.0
            print "   adding rotation to the model..."
            # defining extrapolated velocity
            extrapolated_velocity = (self._one + self._omega) * self._v0 \
                                    - self._omega * self._v00
            # set Coriolis term
            if self._space_dim == 2:
                coriolis_term = cross_2D(extrapolated_velocity)
            elif self._space_dim == 3:
                from dolfin import cross
                coriolis_term = cross(self._rotation_vector,
                                      extrapolated_velocity)
            rhs_momentum -= self._coefficients[0] * dot(coriolis_term, del_v) * dV
        # 2e) rhs momentum equation: buoyancy term
        if self._parameters.buoyancy is True:
            assert self._coefficients[2] != 0.0
            print "   adding buoyancy to the model..."
            # defining extrapolated temperature
            extrapolated_temperature = (self._one + self._omega) * self._T0 - self._omega * self._T00
            # buoyancy term
            rhs_momentum += self._coefficients[2] * extrapolated_temperature * dot(self._gravity, del_v) * dV
        #=======================================================================        
        # 3) lhs energy equation
        lhs_energy = a[0] / timestep * dot(dT, del_T) * dV \
                   + self._coefficients[3] * a_op(dT, del_T) * dV
        # 4a) rhs energy equation: time derivative
        rhs_energy = -dot(a[1] / timestep * self._T0 + a[2] / timestep * self._T00, del_T) * dV
        # 4b) rhs energy equation: nonlinear term
        nonlinear_term_temperature =  b[0] * dot(self._v0, grad(self._T0)) \
                                    + b[1] * dot(self._v00, grad(self._T00))
        rhs_energy -= nonlinear_term_temperature * del_T * dV
        # 4c) rhs energy equation: linear term
        rhs_energy -= self._coefficients[3] \
                        * dot(c[1] * grad(self._T0) + c[2] * grad(self._T00), grad(del_T)) * dV
        #=======================================================================        
        # full problem
        self._lhs = lhs_momentum + lhs_energy
        self._rhs = rhs_momentum + rhs_energy
        if not hasattr(self, "_bcs"):
            self._setup_bcs()
        if self._parameters.use_assembler_method:
            # system assembler
            self._assembler = dlfn.SystemAssembler(self._lhs, self._rhs,
                                                   bcs=self._dirichlet_bcs)
            self._system_matrix = dlfn.Matrix()
            self._system_rhs = dlfn.Vector()
            self._solver = dlfn.LUSolver(self._system_matrix)
        else:
            # linear problem
            problem = dlfn.LinearVariationalProblem(self._lhs, self._rhs, self._sol,
                                                    bcs=self._dirichlet_bcs)
            self._solver = dlfn.LinearVariationalSolver(problem)
            
    def _normalize_pressure(self, pressure):
        assert hasattr(self, "_mesh")
        assert hasattr(self, "_elemP")
        dV = dlfn.Measure("dx", domain = self._mesh)
        V = dlfn.assemble( self._one * dV)
        mean_pressure = dlfn.assemble(pressure * dV) / V
        return dlfn.project(pressure - mean_pressure,
                            dlfn.FunctionSpace(self._mesh, self._elemP))

    def _update_imex_coefficients(self, step=0, omega=1.0):
        assert hasattr(self, "_imex")
        if not hasattr(self, "_imex_alpha"):
            self._imex_alpha = [dlfn.Constant(0.), dlfn.Constant(0.), dlfn.Constant(0.)]
        if not hasattr(self, "_imex_beta"):
            self._imex_beta = [dlfn.Constant(0.), dlfn.Constant(0.)]
        if not hasattr(self, "_imex_gamma"):
            self._imex_gamma = [dlfn.Constant(0.), dlfn.Constant(0.), dlfn.Constant(0.)]
        if step != 0:
            alpha = self._imex.alpha(omega)
            beta = self._imex.beta(omega)
            gamma = self._imex.gamma(omega)
            for i in xrange(3):
                self._imex_alpha[i].assign(alpha[i])
                self._imex_gamma[i].assign(gamma[i])
            for i in xrange(2):
                self._imex_beta[i].assign(beta[i])
        else:
            # Crank Nicholson scheme
            self._imex_alpha[0].assign(1.0)
            self._imex_alpha[1].assign(-1.0)
            self._imex_beta[1].assign(1.0)
            self._imex_gamma[0].assign(0.5)
            self._imex_gamma[1].assign(0.5)
    
    def _compute_cfl_number(self, velocity, time_step):
        assert hasattr(self, "_mesh")
        assert(isinstance(velocity, dlfn.Function))
        assert(isinstance(time_step, float) and time_step > 0.0)
        Vh = dlfn.FunctionSpace(self._mesh, "DG", 0)
        h_cell = dlfn.CellDiameter(self._mesh)
        from dolfin import sqrt, inner
        cfl_cell_expression = sqrt(inner(velocity, velocity)) * time_step / (self._parameters.velocity_degree * h_cell)
        local_cfl = dlfn.project(cfl_cell_expression, Vh)
        return local_cfl.vector().max()    
    
    def _compute_global_averages(self, velocity, temperature):
        assert hasattr(self, "_mesh")
        dV = dlfn.Measure("dx", domain = self._mesh)
        V = dlfn.assemble( self._one * dV)

        from dolfin import dot
        velocity_sqrd = dlfn.assemble(dot(velocity, velocity)* dV)
        temperature_sqrd = dlfn.assemble(temperature * temperature * dV)

        from math import sqrt
        kinetic_energy = 0.5 * velocity_sqrd / V
        velocity_rms = sqrt(velocity_sqrd / V)
        temperature_rms = sqrt(temperature_sqrd / V)
        
        return velocity_rms, kinetic_energy, temperature_rms

            
    def _output_global_averages(self, velocity, temperature, step, time, initial=False):
        if not hasattr(self, "_global_avg_log"):
            self._global_avg_log = []
        
        velocity_rms, kinetic_energy, temperature_rms \
        = self._compute_global_averages(velocity, temperature)
        
        if not initial:
            self._global_avg_log.append((step, time, velocity_rms, kinetic_energy, temperature_rms))
        
        print "   velocity rms:\t\t{0:.6f}".format(velocity_rms)
        print "   kinetic energy:\t{0:.6f}".format(kinetic_energy)
        print "   temperature rms:\t{0:.6f}".format(temperature_rms)
        
    
    def _output_results(self, velocity, temperature, pressure, step, time, initial = False):
        if not hasattr(self, "_pvd_velocity"):
            self._pvd_velocity = dlfn.File("pvd/solution-velocity.pvd")
        if not hasattr(self, "_pvd_pressure"):
            self._pvd_pressure = dlfn.File("pvd/solution-pressure.pvd")
        if not hasattr(self, "_pvd_temperature"):
            self._pvd_temperature = dlfn.File("pvd/solution-temperature.pvd")
        if not hasattr(self, "_pvd_logging"):
            self._pvd_logging = []
        if not initial:
            # normalize pressure to mean value pressure
            normalized_pressure = self._normalize_pressure(pressure)
            # write output to files
            print "   writing output..."
            self._pvd_velocity << (velocity, time)
            self._pvd_pressure << (normalized_pressure, time)
            self._pvd_temperature << (temperature, time)
            # write log
            self._pvd_logging.append((step, time))
        else:
            pvd_velocity = dlfn.File("pvd/initial-velocity.pvd")
            pvd_pressure = dlfn.File("pvd/initial-pressure.pvd")
            pvd_temperature = dlfn.File("pvd/initial-temperature.pvd")
            # write output to files
            print "   writing output..."
            pvd_velocity << velocity
            pvd_pressure << pressure
            pvd_temperature << temperature
            
    def _write_checkpoint(self, step, time):
        assert hasattr(self, "_mesh")
        assert hasattr(self, "_sol")
        
        if not hasattr(self, "_checkpoint_log"):
            self._checkpoint_log = []
            
        print "   writing checkpoint..."
        # write current solution for restart
        fname = "checkpoints/" + "solution-{:06d}.h5".format(step)
        output_file = dlfn.HDF5File(self._mesh.mpi_comm(), fname, "w")
        output_file.write(self._sol, "solution")
        output_file.close()
        # write old solution for restart
        fname = "checkpoints/" + "old_solution-{:06d}.h5".format(step)
        output_file = dlfn.HDF5File(self._mesh.mpi_comm(), fname, "w")
        output_file.write(self._sol0, "old_solution")
        output_file.close()
        # save checkpoint time
        self._checkpoint_log.append((step, time))

    
    def _print_info_message(self):
        print """------------------------------------------------------------------------
Boussinesq solver by S. Glane and B. Turan
        
This program solves the Navier-Stokes system with thermal convection.
A direct solver and the stable Taylor-Hood (P2-P1) element are used.

The governing equations are

\t- Incompressibility constraint:
\t\t div(v) = 0,
\t- Navier-Stokes equation:
\t\tdv/dt + v . grad(v) + C1 Omega .times. v 
\t\t\t\t= - grad(p) + C2 div(grad(v)) - C3 T g,

\t- Heat conduction equation:
\t\tdT/dt + v . grad(T) = C4 div(grad(T)).

The coefficients C1 to C4 depend on the normalization.

+-------------------+---------+---------------+---------+-------------------+
|   Normalization   |    C1   |      C2       |   C3    |        C4         |
+-------------------+---------+---------------+---------+-------------------+
| Non-rotating case |   0     | sqrt(Pr / Ra) | 1       | 1 / sqrt(Ra * Pr) |
| Rotating case     | 2 / Ek  | 1             | Ra / Pr | 1 / Pr            |
+-------------------+---------+---------------+---------+-------------------+"""
        ekman = self._parameters.ekman
        rayleigh = self._parameters.rayleigh
        prandtl = self._parameters.prandtl
        if self._parameters.rotation:
            print """You have chosen the rotating case with:
        Ek = {0:3.2e},
        Ra = {1:3.2e},
        Pr = {2:3.2e}.\n""".format(ekman,rayleigh,prandtl)
        else:
            print """You have chosen the non-rotating case with:
        Ra = {0:3.2e},
        Pr = {1:3.2e}.\n""".format(rayleigh, prandtl)
            
        tmp = self._parameters.coefficients()
        print """The related coefficients C1 to C4 are given by: \t
        C1 = {0:3.2e},
        C2 = {1:3.2e},
        C3 = {2:3.2e},
        C4 = {3:3.2e}.\n""".format(*tmp)
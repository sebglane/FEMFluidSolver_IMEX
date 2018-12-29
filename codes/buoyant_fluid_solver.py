#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import dolfin as dlfn
from time_stepping import IMEXCoefficients

from enum import Enum

class VelocityBCType(Enum):
    no_slip = 0
    slip = 1
    constant = 2
    function = 3

class TemperatureBCType(Enum):
    constant = 0
    function = 1

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

    def run(self):
        self._print_info_message()        
        
        self._setup_function_spaces()
        
        self._setup_boundary_conditions()
        
        self._update_imex_coefficients()
        
        self._setup_problem()
        
        self._setup_initial_condition()
        
        # initialize timestep
        timestep = self._parameters.timestep
        old_timestep = self._parameters.timestep
        # initialize timestepping
        from time_stepping import TimestepControl
        timestep_control = TimestepControl((self._parameters.min_cfl,
                                            self._parameters.max_cfl,),
                                           (self._parameters.min_timestep,
                                            self._parameters.max_timestep,))
        # update flag
        update_imex = False
        # dolfin timer 
        dlfn.tic()
        # time loop
        time = 0.0
        step = 0
        while time < self._paramters.t_end and step < self._parameters.n_steps:
            print "Iteration: {0:08d}, time = {1:10.5f}, time step = {2:5.4e}"\
                  .format(step, time, timestep)
            if update_imex:
                self._update_imex_coefficients(step, timestep / old_timestep)
            # compute solution
            self._solve()
            # extract solution
            velocity, pressure, temperature = self._sol.split()
            # compute cfl number
            cfl = self._compute_cfl_number(velocity, timestep)
            print "   current cfl number: {0:03.2e}".format(cfl)
            # reject step if clf number > 1.0
            if cfl > 1.0:
                print "   time step rejected".format(cfl)
                # decrease time step
                timestep /= 10.0
                # recompute time step
                continue
            #===========================================================================
            # modify time step
            if step != 0 and self._parameters.adaptive_time_stepping:
                old_timestep = timestep
                update_imex, timestep \
                    = timestep_control.adjust_time_step(cfl, timestep)
                self._timestep.assign(timestep)
            elif step == 0:
                # update coefficients after initial step to switch to SBF2
                old_timestep = timestep
                update_imex = True
            # write output
            if step % self._parameters.output_frequency == 0:
                self._output_results(velocity, temperature, step, time)
            # compute rms values
            if step % self._parameters.rms_frequency == 0:
                self._output_rms_values(velocity, temperature, step, time)
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
        if self._parameters.use_assembler_method:
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
            initial_pressure = dlfn.interpolate(0.,
                                                self._Wh.sub(1).collapse())
        # temperature part
        if self._ics.has_key("pressure"):
            initial_temperature_expr = self._ics["temperature"]
            assert isinstance(initial_temperature_expr, dlfn.Expression)
            assert len(initial_temperature_expr.ufl_shape) == 0
            initial_temperature = dlfn.interpolate(initial_temperature_expr,
                                                   self._Wh.sub(2).collapse())
        else:
            initial_temperature = dlfn.interpolate(0.,
                                                   self._Wh.sub(2).collapse())
        # assign initial condition
        dlfn.assign(self._sol0, [initial_velocity,
                                 initial_pressure,
                                 initial_temperature])


    def _update_imex_coefficients(self, step=0, omega=1.0):
        assert hasattr(self, "_imex")
        if not hasattr(self, "_imex_alpha"):
            self._imex_alpha = [dlfn.Constant(0.), ] * 3
        if not hasattr(self, "_imex_beta"):
            self._imex_beta = [dlfn.Constant(0.), ] * 2
        if not hasattr(self, "_imex_gamma"):
            self._imex_gamma = [dlfn.Constant(0.), ] * 3
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

    def _setup_problem(self):
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
        # define auxiliary operators
        from dolfin import inner, div, grad, dot
        def a_op(phi, psi):
            return inner(grad(phi), grad(psi)) * dV
        def b_op(phi, psi):
            return div(phi) * psi * dV
        # 2D cross product
        def cross_2D(phi):
            return dlfn.as_vector([phi[1], -phi[0]])
        #=======================================================================
        # trial and test function
        (del_v, del_p, del_T) = dlfn.TestFunctions(self._Wh)
        (dv, dp, dT) = dlfn.TrialFunctions(self._Wh)
        # volume element
        dV = dlfn.Measure("dx", domain = self._mesh)
        # old solutions
        v0 = self._v0; v00 = self._v00
        T0 = self._T0; T00 = self._T00
        # reference to time step
        timestep = self._timestep
        #=======================================================================
        # 1) lhs momentum equation
        lhs_momentum = a[0] * dot(dv, del_v) *  dV \
                     + timestep * c[0] * self._coefficients[1] * a_op(dv, del_v) \
                     - timestep * b_op(del_v, dp) \
                     - timestep * b_op(dv, del_p)
        # 2a) rhs momentum equation: time derivative
        time_derivative_term_velocity =  a[1] * v0 \
                                       + a[2] * v00
        rhs_momentum = - dot(time_derivative_term_velocity, del_v) * dV
        # 2b) rhs momentum equation: nonlinear term
        nonlinear_term_velocity =  b[0] * dot(grad(v0), v0) \
                                 + b[1] * dot(grad(v00), v00)
        rhs_momentum -= timestep * dot(nonlinear_term_velocity, del_v) * dV
        # 2c) rhs momentum equation: linear term
        linear_term_velocity =  c[1] * grad(v0) \
                              + c[2] * grad(v00)
        rhs_momentum -= timestep * self._coefficients[1] \
                            * inner(linear_term_velocity, grad(del_v)) * dV
        # 2d) rhs momentum equation: coriolis term
        if self._parameters.rotation is True:
            assert self._coefficients[0] != 0.0
            print "   adding rotation to the model...\n"
            # defining extrapolated velocity
            extrapolated_velocity = (self._one + self._omega) * v0 \
                                    - self._omega * v00
            # set Coriolis term
            if self._space_dim == 2:
                coriolis_term = cross_2D(extrapolated_velocity)
            elif self._space_dim == 3:
                from dolfin import cross
                coriolis_term = cross(self._rotation_vector,
                                      extrapolated_velocity)
            rhs_momentum -= timestep * self._coefficients[0] \
                                * dot(coriolis_term, del_v) * dV
        # 2e) rhs momentum equation: buoyancy term
        if self._parameters.buoyancy is True:
            assert self._coefficients[2] != 0.0
            print "   adding buoyancy to the model...\n"
            # defining extrapolated temperature
            extrapolated_temperature = (self._one + self._omega) * T0 \
                                       - self._omega * T00
            # gravity vector
            gravity = self._get_gravity_vector()
            #
            rhs_momentum += timestep * self._coefficients[2] \
                            * extrapolated_temperature * dot(gravity, del_v) * dV
        #=======================================================================        
        # 3) lhs energy equation
        lhs_energy = a[0] * dot(dT, del_T) * dV \
                   + timestep * self._coefficients[3] * a_op(dT, del_T)
        # 4a) rhs energy equation: time derivative
        time_derivative_term_temperature =  a[1] * T0 \
                                          + a[2] * T00
        rhs_energy = -dot(time_derivative_term_temperature, del_T) * dV
        # 4b) rhs energy equation: nonlinear term
        nonlinear_term_temperature =  b[0] * dot(v0, grad(T0)) \
                                    + b[1] * dot(v00, grad(T00))
        rhs_energy -= timestep * nonlinear_term_temperature * del_T * dV
        # 4c) rhs energy equation: linear term
        linear_term_temperature =  c[1] * grad(T0) \
                                 + c[2] * grad(T00)
        rhs_energy -= timestep * self._coefficients[3] \
                        * dot(linear_term_temperature, grad(del_T)) * dV
        #=======================================================================        
        # full problem
        lhs = lhs_momentum + lhs_energy
        rhs = rhs_momentum + rhs_energy
        if not hasattr(self, "_bcs"):
            self._setup_bcs()
        if self._parameters.use_assembler_method:
            # system assembler
            self._assembler = dlfn.SystemAssembler(lhs, rhs,
                                                   bcs=self._dirichlet_bcs)
            self._system_matrix = dlfn.Matrix()
            self._system_rhs = dlfn.Vector()
            self._solver = dlfn.LUSolver(self._system_matrix)
        else:
            # linear problem
            problem = dlfn.LinearVariationalProblem(lhs, rhs, self._sol,
                                                    bcs=self._dirichlet_bcs)
            self._solver = dlfn.LinearVariationalSolver(problem)
    
    def _compute_cfl_number(self, velocity, time_step):
        assert hasattr(self, "_mesh")
        assert(isinstance(velocity, dlfn.Function))
        assert(isinstance(time_step, float) and time_step > 0.0)
        Vh = dlfn.FunctionSpace(self._mesh, "DG", 0)
        h_cell = dlfn.CellDiameter(self._mesh)
        from dolfin import sqrt, inner
        cfl_cell_expression = sqrt(inner(velocity, velocity)) * time_step / h_cell
        local_cfl = dlfn.project(cfl_cell_expression, Vh)
        return local_cfl.vector().max()
    
    def _get_gravity_vector(self):
        # constant gravity in radial direction
        assert hasattr(self, "_space_dim")
        from parameters import GravityType
        if self._parameters.gravity_type == GravityType.radial:
            if self._space_dim == 2:
                gravity_vector = dlfn.Expression(
                                ('x[0] / sqrt( pow(x[0],2) + pow(x[1],2) )', 
                                 'x[1] / sqrt( pow(x[0],2) + pow(x[1],2) )'), degree=2)
            elif self._space_dim == 3:
                gravity_vector = dlfn.Expression(
                                ('x[0] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )', 
                                 'x[1] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )',
                                 'x[2] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )'),
                                degree=2)
            return gravity_vector
        else:
            raise NotImplementedError

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
        
        +-------------------+-------+---------------+---------+-------------------+
        |   Normalization   |    C1 |      C2       |   C3    |        C4         |
        +-------------------+-------+---------------+---------+-------------------+
        | Non-rotating case | 0     | sqrt(Pr / Ra) | 1       | 1 / sqrt(Ra * Pr) |
        | Rotating case     | Ek    | 1             | Ra / Pr | 1 / Pr            |
        +-------------------+-------+---------------+---------+-------------------+"""
        ekman = self._parameters.ekman
        rayleigh = self._parameters.rayleigh
        prandtl = self._parameters.prandtl
        if self._parameters.rotation:
            print """You have chosen the rotating case with: \t
            Ek = {0:3.2e},\t 
            Ra = {1:3.2e},\t
            Pr = {2:3.2e} .""".format(ekman,rayleigh,prandtl)
        else:
            print """You have chosen the non-rotating case with: \t
            Ra = {0:3.2e},\t Pr = {1:3.2e} .""".format(rayleigh, prandtl)
            
        tmp = self._parameters.coefficients()
        print """
        The related coefficients C1 to C4 are given by: \t
        C1 = {0:3.2e},\t
        C2 = {1:3.2e},\t 
        C3 = {2:3.2e},\t 
        C4 = {3:3.2e} .\n""".format(*tmp)
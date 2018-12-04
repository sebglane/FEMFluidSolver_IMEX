# -*- coding: utf-8 -*-
from enum import Enum

class GravityType(Enum):
    unit_vector_x = 0
    unit_vector_y = 1
    unit_vector_z = 2
    radial = -1
    
class ParameterHandler(object):
    def __init__(self):
        # physical parameters
        self.__rotation = False
        self.__buoyancy = True
        self.__gravity_type = None
        self.__ekman = 1.0
        self.__froude = 1.0
        self.__rayleigh = 1.0
        self.__reynolds = 1.0
        self.__prandtl = 1.0
        
        # discretization parameters
        self.__velocity_degree = 2
        self.__pressure_degree = 1
        self.__temperature_degree = 1
        
        # solver options
        self.__use_assembler_method = True
        
        # time stepping parameters
        from time_stepping import IMEXType
        self.__imex_type = IMEXType.CNAB
        self.__n_steps = 1
        self.__timestep = 1.0
        self.__t_end = 1.0
        
        # time stepping control parameters
        self.__adaptive_timestepping = True
        self.__max_timestep = None
        self.__min_timestep = None
        self.__max_cfl = None
        self.__min_cfl = None
        
        # monitoring frequencies
        self.__output_frequency = 100
        self.__rms_frequency = 10
        self.__checkpoint_frequency = 100
    
    @property
    def rotation(self):
        return self.__rotation
        
    @rotation.setter
    def rotation(self, x):
        assert isinstance(x, bool)
        self.__rotation = x

    @property
    def gravity_type(self):
        return self.__gravity_type
    
    @gravity_type.setter
    def gravity_type(self, x):
        assert isinstance(x, GravityType)
        self.__gravity_type = x
    
    @property
    def buoyancy(self):
        return self.__buoyancy
        
    @buoyancy.setter
    def buoyancy(self, x):
        assert isinstance(x, bool)
        self.__buoyancy = x

    @property
    def ekman(self):
        return self.__ekman

    @ekman.setter
    def ekman(self, x):
        assert self.__rotation is True
        assert isinstance(x, float) and x > 0.0
        self.__ekman = x
    
    @property
    def froude(self):
        return self.__froude
    
    @froude.setter
    def froude(self, x):
        assert isinstance(x, float) and x > 0.0
        self.__ekman = x

    @property
    def rayleigh(self):
        return self.__rayleigh
    
    @rayleigh.setter
    def rayleigh(self, x):
        assert self.__buoyancy is True
        assert isinstance(x, float) and x > 0.0
        self.__rayleigh = x

    @property
    def reynolds(self):
        return self.__reynolds
    
    @reynolds.setter
    def reynolds(self, x):
        assert isinstance(x, float) and x > 0.0
        self.__reynolds = x
    
    @property
    def prandtl(self):
        return self.__prandtl
    
    @prandtl.setter
    def prandtl(self, x):
        assert isinstance(x, float) and x > 0.0
        self.__prandtl = x
    
    @property
    def n_steps(self):
        return self.__n_steps
        
    @n_steps.setter
    def n_steps(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__n_steps = x

    @property
    def timestep(self):
        return self.__timestep
        
    @timestep.setter
    def timestep(self, x):
        assert isinstance(x, float)  and x >= 0.0
        self.__timestep = x

    @property
    def t_end(self):
        return self.__t_end
        
    @t_end.setter
    def t_end(self, x):
        assert isinstance(x, float)  and x >= 0.0
        self.__t_end = x
    
    @property
    def adaptive_timestepping(self):
        return self.__adaptive_timestepping
        
    @adaptive_timestepping.setter
    def adaptive_timestepping(self, x):
        assert isinstance(x, bool)
        self.__adaptive_timestepping = x
        
    @property
    def max_timestep(self):
        return self.__max_timestep
        
    @max_timestep.setter
    def max_timestep(self, x):
        assert self.__adaptive_timestepping is True
        assert isinstance(x, float) and x > 0.0
        if self.__min_timestep is not None:
            assert self.__min_timestep < x
        self.__max_timestep = x
    
    @property
    def min_timestep(self):
        return self.__min_timestep
        
    @min_timestep.setter
    def min_timestep(self, x):
        assert self.__adaptive_timestepping is True
        assert isinstance(x, float) and x > 0.0
        if self.__max_timestep is not None:
            assert x < self.__max_timestep
        self.__min_timestep = x
    
    @property
    def max_cfl(self):
        return self.__max_cfl
        
    @max_cfl.setter
    def max_cfl(self, x):
        assert self.__adaptive_timestepping is True
        assert isinstance(x, float) and x > 0.0
        if self.__min_cfl is not None:
            assert self.__min_cfl < x
        self.__max_cfl = x
    
    @property
    def min_cfl(self):
        return self.__min_cfl
        
    @min_cfl.setter
    def min_cfl(self, x):
        assert self.__adaptive_timestepping is True
        assert isinstance(x, float) and x > 0.0
        if self.__max_cfl is not None:
            assert x < self.__max_cfl
        self.__min_cfl = x

    @property
    def output_frequency(self):
        return self.__output_frequency
        
    @output_frequency.setter
    def output_frequency(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__output_frequency = x

    @property
    def rms_frequency(self):
        return self.__rms_frequency
        
    @rms_frequency.setter
    def rms_frequency(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__rms_frequency = x

    @property
    def checkpoint_frequency(self):
        return self.__checkpoint_frequency
        
    @checkpoint_frequency.setter
    def checkpoint_frequency(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__checkpoint_frequency = x

    @property
    def use_assembler_method(self):
        return self.__use_assembler_method
        
    @use_assembler_method.setter
    def use_assembler_method(self, x):
        assert isinstance(x, bool)
        self.__use_assembler_method = x

    @property
    def velocity_degree(self):
        return self.__velocity_degree
        
    @velocity_degree.setter
    def velocity_degree(self, x):
        assert isinstance(x, int)  and x >= 1 and \
            self.__velocity_degree >= self.__pressure_degree + 1
        self.__velocity_degree = x
    
    @property
    def pressure_degree(self):
        return self.__pressure_degree
        
    @pressure_degree.setter
    def pressure_degree(self, x):
        assert isinstance(x, int)  and x >= 1 and \
            self.__velocity_degree >= self.__pressure_degree + 1
        self.__pressure_degree = x

    @property
    def temperature_degree(self):
        return self.__temperature_degree
        
    @temperature_degree.setter
    def temperature_degree(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__temperature_degree = x
        
    def coefficients(self):
        # equation coefficients
        if self.__rotation is True:
            tmp = (2.0 / self.__ekman,
                   1.0,
                   self.__rayleigh / self.__prandtl,
                   1.0 / self.__prandtl)
        else:
            from math import sqrt
            tmp = (0.0,
                   sqrt(self.__prandtl/ self.__rayleigh),
                   1.0,
                   1.0 / sqrt(self.__rayleigh * self.__prandtl))
        return tmp
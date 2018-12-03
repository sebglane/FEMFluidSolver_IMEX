# -*- coding: utf-8 -*-
class ParameterHandler(object):
    def __init__(self):
        # physical parameters
        self.__rotation = False
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
        self.__n_steps = 1
        self.__t_end = 1.0
        self.__adaptive_timestepping = True
        
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
    def ekman(self):
        return self.__rotation

    @ekman.setter
    def ekman(self, x):
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
    def adaptive_timestepping(self):
        return self.__adaptive_timestepping
        
    @adaptive_timestepping.setter
    def adaptive_timestepping(self, x):
        assert isinstance(x, bool)
        self.__adaptive_timestepping = x

    @property
    def n_steps(self):
        return self.__n_steps
        
    @n_steps.setter
    def n_steps(self, x):
        assert isinstance(x, int)  and x >= 1
        self.__n_steps = x

    @property
    def t_end(self):
        return self.__t_end
        
    @t_end.setter
    def t_end(self, x):
        assert isinstance(x, float)  and x >= 0.0
        self.__t_end = x

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
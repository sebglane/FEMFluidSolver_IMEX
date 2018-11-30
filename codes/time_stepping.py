#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from enum import Enum

class IMEXType(Enum):
    CNAB = 0
    MCNAB = 1
    CNLF = 2
    SBDF = 3

class IMEXCoefficients():
    def __init__(self, type_string = "CNAB"):
        imex_types = ("CNAB", "MCNAB", "CNLF", "SBDF")
        assert type_string in imex_types
        self._type = type_string

        self._omega = -1.
        
        self._alpha = [0., 0., 0.]
        self._beta = [0., 0.]
        self._gamma = [0., 0., 0.]
        
        self._update_alpha = True
        self._update_beta = True
        self._update_gamma = True
    
    def _compute_alpha(self):
        if not self._update_alpha:
            return
        if self._type is "SBDF":
            self._alpha[0] = (1. + 2. * self._omega) / (1. + self._omega)
            self._alpha[1] = -(1. + self._omega)
            self._alpha[2] = (self._omega * self._omega) / (1. + self._omega)
        elif self._type is "CNAB" or self._type is "MCNAB":
            self._alpha[0] = 1.0            
            self._alpha[1] = -1.0
            self._alpha[2] = 0.0
        elif self._type is "CNLF":
            self._alpha[0] = 1. / (1. + self._omega);
            self._alpha[1] = self._omega - 1.;
            self._alpha[2] = -(self._omega * self._omega) / (1. + self._omega);
        self._update_alpha = False

    def _compute_beta(self):
        if not self._update_beta:
            return
        if self._type is "SBDF":
            self._beta[0] = (1. + self._omega)
            self._beta[1] = -self._omega
        elif self._type is "CNAB" or self._type is "MCNAB":
            self._beta[0] = (1. + 0.5 * self._omega)
            self._beta[1] = -0.5 * self._omega
        elif self._type is "CNLF":
            self._beta[0] = 1.0
            self._beta[1] = 0.0
        self._update_gamma = False

    def _compute_gamma(self):
        if not self._update_gamma:
            return
        if self._type is "SBDF":
            self._gamma[0] = 1.0
            self._gamma[1] = 0.0
            self._gamma[2] = 0.0
        elif self._type is "CNAB":
            self._gamma[0] = 0.5
            self._gamma[1] = 0.5
            self._gamma[2] = 0.0
        elif self._type is "MCNAB":
            self._gamma_[0] = (8. * self._omega + 1.)/ (16. * self._omega);
            self._gamma_[1] = (7. * self._omega - 1.)/ (16. * self._omega);
            self._gamma_[2] = self._omega / (16. * self._omega);
        elif self._type is "CNLF":
            self._gamma[0] = 0.5 / self._omega
            self._gamma[1] = 0.5 * (1. - 1./self._omega)
            self._gamma[2] = 0.5
        self._update_gamma = False

    def alpha(self, timestep_ratio):
        assert isinstance(timestep_ratio, float)
        assert timestep_ratio > 0.0
        
        if timestep_ratio != self._omega:
            self._omega = timestep_ratio

            self._update_alpha = True
            self._update_beta = True
            self._update_gamma = True

        self._compute_alpha()

        return self._alpha

    def beta(self, timestep_ratio):
        assert isinstance(timestep_ratio, float)
        assert timestep_ratio > 0.0
        
        if timestep_ratio != self._omega:
            self._omega = timestep_ratio

            self._update_alpha = True
            self._update_beta = True
            self._update_gamma = True

        self._compute_beta()

        return self._beta
    

    def gamma(self, timestep_ratio):
        assert isinstance(timestep_ratio, float)
        assert timestep_ratio > 0.0
        
        if timestep_ratio != self._omega:
            self._omega = timestep_ratio

            self._update_alpha = True
            self._update_beta = True
            self._update_gamma = True

        self._compute_gamma()

        return self._gamma

class TimestepControl:
    def __init__(self, cfl_interval, timestep_interval,
                 save_history = False, print_on_update = True):
        # input check for cfl interval
        assert isinstance(cfl_interval, tuple) \
            and all(isinstance(x, float) for x in cfl_interval)
        cfl_min, cfl_max = cfl_interval
        assert cfl_min < cfl_max and cfl_min > 0. and cfl_max > 0.
        self._cfl_min = cfl_min
        self._cfl_max = cfl_max
        # input check for timestep interval
        assert isinstance(timestep_interval, tuple) \
            and all(isinstance(x, float) for x in timestep_interval)
        dt_min, dt_max = timestep_interval
        assert dt_min < dt_max and dt_min > 0. and dt_max > 0.
        self._dt_min = dt_min
        self._dt_max = dt_max
        # input check for history flag
        assert isinstance(save_history, bool)
        self._save_history = save_history
        # input check for printing flag
        assert isinstance(print_on_update, bool)
        self._print = print_on_update
        
        self._old_timestep = 0.
        self._timestep_history = []
        self._cfl_history = []
        
    def adjust_time_step(self, cfl, timestep):
        # input check
        assert isinstance(cfl, float) and cfl > 0.0
        assert isinstance(timestep, tuple) and timestep > 0.0
        # set boolean flags
        timestep_modified = False
        # compute intended timestep
        dt_cfl = 0.5 * (self._cfl_min + self._cfl_max) / cfl * self._old_timestep
        # initialize timestep variable
        new_timestep = -1.0
        # check if time step should be modified
        # case 1: cfl criterion
        if dt_cfl < self._dt_max and\
            (cfl > self._cfl_max or cfl < self._cfl_min):
            # case 1a: lower limit is violated
            if dt_cfl < self._dt_min:
                raise ValueError("time step is too small, aborting this run")
            # case 1b: timestep is modified
            elif dt_cfl != self._old_timestep:
                new_timestep = dt_cfl
                timestep_modified = True
            # case 1c: timestep is not modified
            elif dt_cfl == self._old_timestep:
                new_timestep = self._old_timestep
            else:
                raise RuntimeError()
        # case 2: cfl criterion gives timestep larger than dt_max
        elif dt_cfl > self._dt_max and self._old_timestep != self._dt_max:
            new_timestep = self._dt_max
            timestep_modified = True
        # case 3: cfl criterion gives timestep larger than dt_max
        elif (dt_cfl > self._dt_max and self._old_timestep == self._dt_max):
            new_timestep = self._old_timestep
        else:
            raise RuntimeError()
        assert new_timestep > 0.
        # update old timestep
        if timestep_modified:
            self._old_timestep = new_timestep
        # print update information
        if timestep_modified and self._print:
            print "   time step modified from " + \
                    "{0:03.2e} to {1:03.2e}".format(self._old_timestep, new_timestep)
        if self._save_history:
            self._timestep_history.append(new_timestep)
            self._cfl_history.append(cfl)
        return timestep_modified, new_timestep
    
    def get_timestep_history(self):
        assert self._save_history
        assert len(self._timestep_history) == len(self._cfl_history)
        
        import numpy as np
        timestep_history = np.array(self._timestep_history, dtype=np.float)
        
        return timestep_history

    def get_cfl_history(self):
        assert self._save_history
        assert len(self._timestep_history) == len(self._cfl_history)
        
        import numpy as np
        cfl_history = np.array(self._cfl_history, dtype=np.float)
        
        return cfl_history
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from enum import Enum

class IMEXType(Enum):
    CNAB = 0
    MCNAB = 1
    CNLF = 2
    SBDF = 3

class IMEXCoefficients():
    def __init__(self, imex_type= IMEXType.CNAB):
        assert isinstance(imex_type, IMEXType)
        self._type = imex_type

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
        if self._type is IMEXType.SBDF:
            self._alpha[0] = (1. + 2. * self._omega) / (1. + self._omega)
            self._alpha[1] = -(1. + self._omega)
            self._alpha[2] = (self._omega * self._omega) / (1. + self._omega)
        elif self._type is IMEXType.CNAB or self._type is IMEXType.MCNAB:
            self._alpha[0] = 1.0            
            self._alpha[1] = -1.0
            self._alpha[2] = 0.0
        elif self._type is IMEXType.CNLF:
            self._alpha[0] = 1. / (1. + self._omega);
            self._alpha[1] = self._omega - 1.;
            self._alpha[2] = -(self._omega * self._omega) / (1. + self._omega);
        self._update_alpha = False

    def _compute_beta(self):
        if not self._update_beta:
            return
        if self._type is IMEXType.SBDF:
            self._beta[0] = (1. + self._omega)
            self._beta[1] = -self._omega
        elif self._type is IMEXType.CNAB or self._type is IMEXType.MCNAB:
            self._beta[0] = (1. + 0.5 * self._omega)
            self._beta[1] = -0.5 * self._omega
        elif self._type is IMEXType.CNLF:
            self._beta[0] = 1.0
            self._beta[1] = 0.0
        self._update_beta = False

    def _compute_gamma(self):
        if not self._update_gamma:
            return
        if self._type is IMEXType.SBDF:
            self._gamma[0] = 1.0
            self._gamma[1] = 0.0
            self._gamma[2] = 0.0
        elif self._type is IMEXType.CNAB:
            self._gamma[0] = 0.5
            self._gamma[1] = 0.5
            self._gamma[2] = 0.0
        elif self._type is IMEXType.MCNAB:
            self._gamma_[0] = (8. * self._omega + 1.)/ (16. * self._omega);
            self._gamma_[1] = (7. * self._omega - 1.)/ (16. * self._omega);
            self._gamma_[2] = self._omega / (16. * self._omega);
        elif self._type is IMEXType.CNLF:
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
            and all(isinstance(x, float) for x in cfl_interval) \
            and len(cfl_interval) == 2
        cfl_min, cfl_max = cfl_interval
        assert cfl_min < cfl_max and cfl_min > 0. and cfl_max > 0.
        self._cfl_min = cfl_min
        self._cfl_max = cfl_max
        # input check for timestep interval
        assert isinstance(timestep_interval, tuple) \
            and all(isinstance(x, float) for x in timestep_interval) \
            and len(timestep_interval) == 2
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
        assert isinstance(timestep, float) and timestep > 0.0
        # compute intended timestep
        dt_cfl = 0.5 * (self._cfl_min + self._cfl_max) / cfl * timestep
        # check if time step should be modified
        # case 1: cfl criterion
        if dt_cfl < self._dt_max and\
            (cfl > self._cfl_max or cfl < self._cfl_min):
            # case 1a: lower limit is violated
            if dt_cfl < self._dt_min:
                raise ValueError("time step is too small, aborting this run\n"
                    + "dt_cfl = {0:3.2e}, dt_min = {1:3.2e}\n".format(dt_cfl, self._dt_min)
                    + "cfl = {0:3.2e}, old_timestep = {1:3.2e}".format(cfl, self._old_timestep))
            # case 1b: timestep is modified
            elif dt_cfl != self._old_timestep:
                self._finalize(True, cfl, dt_cfl)
                return True, dt_cfl, dt_cfl / timestep
            # case 1c: timestep is not modified
            elif dt_cfl == self._old_timestep:
                self._finalize(False, cfl)
                return False, self._old_timestep, None
            else:
                raise RuntimeError()
        # case 2: cfl criterion gives timestep larger than dt_max
        elif dt_cfl > self._dt_max and self._old_timestep != self._dt_max:
            self._finalize(True, cfl, self._dt_max)
            return True, self._dt_max, self._dt_max / timestep
        # case 3: cfl criterion gives timestep larger than dt_max
        elif (dt_cfl > self._dt_max and self._old_timestep == self._dt_max):
            self._finalize(False, cfl)
            return False, self._old_timestep, None
        elif dt_cfl < self._dt_max and cfl < self._cfl_max and cfl > self._cfl_min:
            self._finalize(False, cfl)
            return False, self._old_timestep, None
        else:
            raise RuntimeError()
            
    def _finalize(self, timestep_modified, cfl, new_timestep = None):
        # print update information
        if timestep_modified and self._print:
            print "   time step modified from " + \
                    "{0:03.2e} to {1:03.2e}".format(self._old_timestep, new_timestep)
        if self._save_history:
            self._timestep_history.append(new_timestep)
            self._cfl_history.append(cfl)
        if timestep_modified:
            self._old_timestep = new_timestep
    
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
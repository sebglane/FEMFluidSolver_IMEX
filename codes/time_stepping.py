#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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
        elif self._type is "CNAB" or self._type is "MCNAB":
            self._gamma[0] = 0.5
            self._gamma[1] = 0.5
            self._gamma[2] = 0.0
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

        
        
            
        
    
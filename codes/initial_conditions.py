#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:38:51 2019

@author: sg
"""

from math import pi

def Christensen2001Case0(dim, radii):
    assert isinstance(dim, int)
    assert dim == 2 or dim == 3
    assert isinstance(radii, (list, tuple)) and len(radii) == 2
    ri, ro = radii
    assert isinstance(ri, float) and ri > 0.
    assert isinstance(ro, float) and ro > 0.
    assert ri < ro
        
    # define coordinates
    if dim == 2:
        radius_string = "sqrt( pow(x[0],2) + pow(x[1],2))"
    elif dim == 3:
        radius_string = "sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2))"

    xi_string = "(2.0 * {0} - ro - ri)/(ro - ri)".format(radius_string)
    azimuthal_fun_string = "((pow(x[0],4) - 6.0 * pow(x[0],2) *  pow(x[1],2) + pow(x[1],4) ) \
            / ( pow(x[0],2) + pow(x[1],2) ))"

    # define functions
    if dim == 2:
        adiabat_string = "(log(ro) - log({0})) / (log(ro) - log(ri))".format(radius_string)
        perturbation_string = "21. / sqrt(17920. * pi) * \
            (1. - 3. * pow({0}, 2) + 3. * pow({0}, 4) - pow({0}, 6)) * {1}".format(xi_string, azimuthal_fun_string)
    elif dim == 3:
        adiabat_string = "(ro - {0}) / (ro - ri) * ri / {0}".format(radius_string)
        polar_fun_string = "( pow(pow(x[0],2) + pow(x[1],2), 2) \
            / pow(pow(x[0],2) + pow(x[1],2) + + pow(x[2],2), 2) )"
        perturbation_string = "21. / sqrt(17920. * pi) * \
            (1. - 3. * pow({0}, 2) + 3. * pow({0}, 4) - pow({0}, 6)) *\
            {1} * {2}".format(xi_string, polar_fun_string, azimuthal_fun_string)
    
    import string    
    expression_string = string.join((adiabat_string, perturbation_string), sep=" + ")

    
    from dolfin import Expression
    ics = dict()
    ics["temperature"] = Expression(expression_string, ri = ri, ro = ro,
                                    pi = pi, degree = 2)
    return ics
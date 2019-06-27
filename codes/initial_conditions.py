#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:38:51 2019

@author: sg
"""

from math import pi
from dolfin import Expression


cpp_code = """
#include <cmath>
#include <cassert>

class TemperatureProfile : public Expression
{
public:
  TemperatureProfile() : Expression(), ri(0.1), ro(1.0) { }

  void eval(Array<double>& values, const Array<double>& x) const
  {
    // spherical coordinates
    const double radius = sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2));
    assert (radius > 0.);

    const double r_cylinder = sqrt(x[0]*x[0] + x[1]*x[1]);
    assert(r_cylinder >= 0.);
    const double theta = atan2(r_cylinder, x[2]);
    assert(theta >= 0. && theta <= 3.14159265358979323846);
    
    const double phi = atan2(x[1], x[0]);
    assert(phi >= -3.14159265358979323846 && phi <= 3.14159265358979323846);
        
    // transformed radial coordinate
    const double xi = (2.0 * radius - ro - ri)/(ro - ri);

    // adiabat    
    const double adiabat_profile = (ro - radius) / (ro - ri) * ri / radius;
    
    // perturbation
    const double perturbation = 21. / sqrt(17920. * 3.14159265358979323846)
                    * (1. - 3. * pow(xi, 2) + 3. * pow(xi, 4) - pow(xi, 6))
                    * pow(sin(theta), 4) * cos(4.*phi);
    
    values[0] = adiabat_profile + perturbation;
  }
public:
    double ri;
    double ro;
};"""
        

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

    # define functions
    ics = dict()
    if dim == 2:
        adiabat_string = "(log(ro) - log({0})) / (log(ro) - log(ri))".format(radius_string)
        azimuthal_fun_string = "( pow(x[0],4) - 6.0 * pow(x[0],2) *  pow(x[1],2) + pow(x[1],4) ) \
            / pow(pow(x[0],2) + pow(x[1],2), 2)"
        perturbation_string = "21. / sqrt(17920. * pi) * \
            (1. - 3. * pow({0}, 2) + 3. * pow({0}, 4) - pow({0}, 6)) * {1}".format(xi_string, azimuthal_fun_string)

        import string    
        expression_string = string.join((adiabat_string, perturbation_string), sep=" + ")
        
        ics["temperature"] = Expression(expression_string, ri = ri, ro = ro, pi = pi, degree = 2)
        
        return ics
    
    elif dim == 3:
        ics["temperature"] = Expression(cpp_code, degree=2)
        ics["temperature"].ri = ri
        ics["temperature"].ro = ro
        
        return ics    
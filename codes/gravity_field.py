#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:56:45 2019

@author: sg
"""

from enum import Enum

class GravityType(Enum):
    unit_vector_x = 0
    unit_vector_y = 1
    unit_vector_z = 2
    radial = 3
    radial_linear = 4

def get_gravity_field(dim, gravity_type, radii = None):
    assert isinstance(dim, int)
    assert dim == 2 or dim == 3
    assert isinstance(gravity_type, GravityType)
    
    # constant gravity in radial direction
    from dolfin import Expression
    if gravity_type == GravityType.radial:
        if dim == 2:
            gravity_vector = Expression(
                            ('-x[0] / sqrt( pow(x[0],2) + pow(x[1],2) )', 
                             '-x[1] / sqrt( pow(x[0],2) + pow(x[1],2) )'), degree=2)
        elif dim == 3:
            gravity_vector = Expression(
                            ('-x[0] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )', 
                             '-x[1] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )',
                             '-x[2] / sqrt( pow(x[0],2) + pow(x[1],2) + pow(x[2],2) )'),
                            degree=2)
        return gravity_vector
    elif gravity_type == GravityType.radial_linear:
        assert radii is not None
        assert isinstance(radii, (list, tuple)) and len(radii) == 2
        ri, ro = radii
        assert isinstance(ri, float) and ri > 0.
        assert isinstance(ro, float) and ro > 0.
        assert ri < ro
        
        if dim == 2:
            gravity_vector = Expression(('-x[0] / ro', '-x[1] / ro'), degree=2, ro = ro)
        elif dim == 3:
            gravity_vector = Expression(('-x[0] / ro', '-x[1] / ro', '-x[2] / ro'),
                                        degree=2, ro = ro)
        return gravity_vector
    else:
        raise NotImplementedError


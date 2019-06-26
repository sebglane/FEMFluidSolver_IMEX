#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:34:51 2019

@author: sg
"""

from enum import Enum

class VelocityBCType(Enum):
    no_slip = 0
    slip = 1
    constant = 2
    function = 3

class TemperatureBCType(Enum):
    constant = 0
    function = 1
    
def Christensen2001Case0(boundary_ids):
    assert isinstance(boundary_ids, (list, tuple)) and len(boundary_ids) == 2
    inner_boundary_id, outer_boundary_id = boundary_ids
    assert isinstance(inner_boundary_id, int) and inner_boundary_id >= 0
    assert isinstance(outer_boundary_id, int) and outer_boundary_id >= 0
    assert inner_boundary_id != outer_boundary_id
    
    bcs = dict()
    bcs["velocity"] = (
            (VelocityBCType.no_slip , inner_boundary_id, None), 
            (VelocityBCType.no_slip , outer_boundary_id, None))
    bcs["temperature"] = (
            (TemperatureBCType.constant, inner_boundary_id, 1.0), 
            (TemperatureBCType.constant, outer_boundary_id, 0.0))
    
    return bcs
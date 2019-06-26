#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:58:09 2019

@author: sg
"""
from enum import Enum
class GeometryType(Enum):
    spherical_annulus = 0
    rectangle = 1
    other = 2

from dolfin import SubDomain
class CircleBoundary(SubDomain):
    def __init__(self, **kwargs):
        self._hmin = kwargs["hmin"]
        self._radius = kwargs["radius"]
    def inside(self, x, on_boundary):
        # tolerance: half length of smallest element
        tol = self._hmin / 2. 
        from math import sqrt
        result = abs(sqrt(x[0]**2 + x[1]**2) - self._radius) < tol
        return result and on_boundary

def spherical_shell(dim, radii, boundary_ids, n_refinements = 0):
    assert isinstance(dim, int)
    assert dim == 2 or dim == 3
    assert isinstance(radii, (list, tuple)) and len(radii) == 2
    ri, ro = radii
    assert isinstance(ri, float) and ri > 0.
    assert isinstance(ro, float) and ro > 0.
    assert ri < ro
    assert isinstance(boundary_ids, (list, tuple)) and len(boundary_ids) == 2
    inner_boundary_id, outer_boundary_id = boundary_ids
    assert isinstance(inner_boundary_id, int) and inner_boundary_id >= 0
    assert isinstance(outer_boundary_id, int) and outer_boundary_id >= 0
    assert inner_boundary_id != outer_boundary_id
    assert isinstance(n_refinements, int) and n_refinements >= 0
    
    # mesh generation
    from dolfin import Point
    center = Point(0., 0.)
    
    from mshr import Circle, generate_mesh
    domain = Circle(center, ro) \
           - Circle(center, ri)
    mesh = generate_mesh(domain, 120)
    
    # mesh refinement
    from dolfin import refine
    for i in range(n_refinements):
        mesh = refine(mesh)
        
    # subdomains for boundaries
    from dolfin import MeshFunctionSizet
    facet_marker = MeshFunctionSizet(mesh, mesh.topology().dim() - 1)
    facet_marker.set_all(0)
    
    
    # size of smallest element
    hmin = mesh.hmin()
    # inner circle boundary
    class InnerCircle(SubDomain):
        def inside(self, x, on_boundary):
            # tolerance: half length of smallest element
            tol = hmin / 2. 
            from math import sqrt
            result = abs(sqrt(x[0]**2 + x[1]**2) - ri) < tol
            return result and on_boundary
    # outer cirlce boundary
    class OuterCircle(SubDomain):
        def inside(self, x, on_boundary):
            # tolerance: half length of smallest element
            tol = hmin/2.
            from math import sqrt
            result = abs(sqrt(x[0]**2 + x[1]**2) - ro) < tol
            return result and on_boundary
    # mark boundaries
    gamma_inner = InnerCircle()
    gamma_inner.mark(facet_marker, inner_boundary_id)
    gamma_outer = OuterCircle()
    gamma_outer.mark(facet_marker, outer_boundary_id)
    
    
#    # mark boundaries
#    gamma_inner = CircleBoundary(hmin=mesh.hmin(), radius=ri)
#    gamma_inner.mark(facet_marker, inner_boundary_id)
#    gamma_outer = CircleBoundary(hmin=mesh.hmin(), radius=ro)
#    gamma_outer.mark(facet_marker, outer_boundary_id)
    
    return mesh, facet_marker
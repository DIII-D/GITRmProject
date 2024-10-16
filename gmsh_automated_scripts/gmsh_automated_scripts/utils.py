#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:06:17 2024

@author: cappellil
"""
import gmsh
#%%
""" miscellaneous functions """

def rectangle_def(x, y, z, width, height):
    # Create points for the rectangle corners
    p1 = gmsh.model.occ.addPoint(x, y, z)          # Bottom-left corner
    p2 = gmsh.model.occ.addPoint(x + width, y, z)  # Bottom-right corner
    p3 = gmsh.model.occ.addPoint(x + width, y + height, z)  # Top-right corner
    p4 = gmsh.model.occ.addPoint(x, y + height, z)  # Top-left corner
    
    # Create lines for the rectangle edges
    l1 = gmsh.model.occ.addLine(p1, p2)  # Bottom edge
    l2 = gmsh.model.occ.addLine(p2, p3)  # Right edge
    l3 = gmsh.model.occ.addLine(p3, p4)  # Top edge
    l4 = gmsh.model.occ.addLine(p4, p1)  # Left edge# Define a rectangle in the XY plane
    
    loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    surface = gmsh.model.occ.addPlaneSurface([loop])
    
    return p1, p2, p3, p4, l1, l2, l3, l4, loop, surface

#%%
def create_loops(input_dict, z_dimes, volumes_surfaces, dot_loops, ax, ay, az, theta_dimes):

    for elem_def in input_dict.values():
    
        if elem_def["shape"] == "circle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            r = elem_def['radius']
            
            dot_shape = gmsh.model.occ.addCircle(x, y, z, r)
            dot_loop = gmsh.model.occ.addCurveLoop([dot_shape])
            dot_surface = gmsh.model.occ.addPlaneSurface([dot_loop])
            gmsh.model.occ.rotate([(2 , dot_surface)], x, y, z, ax, ay, az, theta_dimes)
                    
        elif elem_def["shape"] == "rectangle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            width = elem_def['width']
            height = elem_def['height']
            
            dot_loop, dot_surface = rectangle_def(x, y, z, width, height)[-2:]
            gmsh.model.occ.rotate([(2 , dot_surface)], x + width / 2, y + height / 2, z, ax, ay, az, theta_dimes)
        
        # store dots loops IDs
        dot_loops.append(dot_loop)
        
        # store surfaces enclosing volume in a variable
        volumes_surfaces.append(gmsh.model.occ.addPlaneSurface([dot_loop]))
        

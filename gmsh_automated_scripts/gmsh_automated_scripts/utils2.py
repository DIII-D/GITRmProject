#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:07:27 2024

@author: cappellil
"""
import gmsh
from .utils import rectangle_def , create_loops

def make_dimes_geom(input_dict, l_radial=12, l_toroidal=12, l_vertical=5, 
                                 x_center_dimes = 0, y_center_dimes =0 , z_center_dimes=0, r_dimes = 1):
    #%% 
    
    """ Plasma volume geometry """
    # Cartesian coordinates of bottom plasma volume surface lower left corner
    
    x_center = 0
    y_center = 0
    z_center = 0
    
    x_plasma_volume_ll = -l_radial / 2  + x_center# toroidal
    y_plasma_volume_ll = -l_toroidal / 2+ y_center # radial
    z_plasma_volume_ll = z_center# vertical
    
    # Width and height of plasma volume base surface 
    width = l_radial
    height = l_toroidal
    dz_plasma_volume = l_vertical  # This variable is defined but not used in this surface creation
    
    # Create a curve loop for plasma volume base
    plasma_base_rectangle_loop, p1, p2, p3, p4, l1, l2, l3, l4 = \
        rectangle_def(x_plasma_volume_ll, y_plasma_volume_ll, z_plasma_volume_ll, width, height) 
    
    # Create a curve loop for plasma_volume top
    
    plasma_top_rectangle_loop, p5, p6, p7, p8, l5, l6, l7, l8 = rectangle_def(x_plasma_volume_ll, y_plasma_volume_ll, \
                                                              z_plasma_volume_ll + dz_plasma_volume, width, height) 
        
    # Create vertical lines connecting bottom and top
    l9 = gmsh.model.occ.addLine(p1, p5)
    l10 = gmsh.model.occ.addLine(p2, p6)
    l11 = gmsh.model.occ.addLine(p3, p7)
    l12 = gmsh.model.occ.addLine(p4, p8)
    
    # list of surfaces IDs enclosing the plasma volume
    volumes_surfaces = []
    
    # Create surfaces for the sides (lateral surfaces)
    plasma_side1 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l1, l10, -l5, -l9])]) # minus when instead of going for extremity A to B you go from B to A close a loop
    plasma_side2 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l2, l11, -l6, -l10])])
    plasma_side3 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l3, l12, -l7, -l11])])
    plasma_side4 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l4, l9, -l8, -l12])])
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_side1)
    volumes_surfaces.append(plasma_side2)
    volumes_surfaces.append(plasma_side3)
    volumes_surfaces.append(plasma_side4)
    
    # Create the top surface
    plasma_top = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])])
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_top)
    
    # Synchronize the GMSH model
    gmsh.model.occ.synchronize()
    
    #%% 
    """ DiMES geometry """
    
    DiMES_circle = gmsh.model.occ.addCircle(x_center_dimes , y_center_dimes , 0 , r_dimes)
    
    DiMES_circle_loop = gmsh.model.occ.addCurveLoop([DiMES_circle])
    
    plasma_base = gmsh.model.occ.addPlaneSurface([plasma_base_rectangle_loop, DiMES_circle_loop]) 
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_base)
    
    if z_center_dimes > 0:                  
    
        DiMES_edge_surface = gmsh.model.occ.extrude([(1,DiMES_circle)], 0, 0, z_center_dimes)
            
        # store surfaces enclosing volume in a variable
        for tup in DiMES_edge_surface:
            if tup[0] == 2:
                DiMES_edge_surface_id = tup[1]
                break
        
        volumes_surfaces.append(DiMES_edge_surface_id)
        
        #identify top_circle loop
        DiMES_top_circle_loop = gmsh.model.occ.addCurveLoop([DiMES_edge_surface[0][1]])
        
    else:
        DiMES_top_circle_loop = DiMES_circle_loop
        
    #%%
    
    """ geometry Dots (coatings)"""
        
    dot_loops = []
    
    create_loops(input_dict, z_center_dimes, volumes_surfaces, dot_loops)

    # Generate DiMES top surface
    gmsh.model.occ.synchronize()
    DiMES_top_surface = gmsh.model.occ.addPlaneSurface([DiMES_top_circle_loop] + dot_loops)
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(DiMES_top_surface)
    
    #%% generate the volume
    plasma_volume = gmsh.model.occ.addVolume([gmsh.model.occ.addSurfaceLoop(volumes_surfaces)])
    

#%% 
""" function """

def make_dimes_mesh(filename = "test.msh" , GUI_geo=False, GUI_msh=True):
    #%% Generate the mesh and visualize the result
    if GUI_geo:
    # Launch the GUI to see the results:
    # Optionally, run the GUI to visualize
        gmsh.fltk.run()
    # Final synchronization of the CAD model
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)
    
    # # Prevent very small elements in small dots
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
    # Set maximum mesh characteristic length for the whole model
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2) 
    gmsh.model.mesh.generate(3)
    
    if GUI_msh:
    # Launch the GUI to see the results:
    # Optionally, run the GUI to visualize
        gmsh.fltk.run()
    
    gmsh.write(filename)
    

#%%
""" function """

def generate_dimes_mesh(input_dict, **kw): 
    try:
        gmsh.finalize()
    except: 
        pass
    finally:
        gmsh.initialize()
        
    make_dimes_geom(input_dict, **kw)
    make_dimes_mesh(**kw)
    
    try:
        gmsh.finalize()
    except:
        pass
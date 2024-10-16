#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:08:13 2024

@author: cappellil
"""
import gmsh

def mesh_generation(filename = "test.msh" , GUI=True):
    #%% Generate the mesh and visualize the result

    # Final synchronization of the CAD model
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)
    
    # # Prevent very small elements in small dots
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
    # Set maximum mesh characteristic length for the whole model
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2) 
    gmsh.model.mesh.generate(3)
    
    if GUI:
    # Launch the GUI to see the results:
    # Optionally, run the GUI to visualize
        gmsh.fltk.run()
        
    gmsh.write(filename)
    
    # Finalize GMSH
    gmsh.finalize()
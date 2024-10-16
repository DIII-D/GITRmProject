#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:54:25 2024

@author: cappellil
"""
#%%

# set working directory
import os
os.chdir('/home/cappellil/GITRmProject/gmsh_automated_scripts')

from gmsh_automated_scripts.utils2 import generate_dimes_mesh 
# Initialize GMSH


#%% 

""" INPUT dots in structure """

N_dots = 3
Names = ["Dot_1","Dot_2","Dot_3"]

input_dict = {
    Names[0]: {
        "shape": "circle",
        "x": 0,
        "y": 0.75,
        "radius": 0.05  # For cylinder
    },
    Names[1]: {
        "shape": "rectangle",
        "x": -0.5,
        "y": -0.25,
        "width": 1,
        "height": 0.5
    },
    Names[2]: {
        "shape": "circle",
        "x": 0,
        "y": -0.75,
        "radius": 0.05
    }
    # Add more dots as needed
}

generate_dimes_mesh(input_dict)
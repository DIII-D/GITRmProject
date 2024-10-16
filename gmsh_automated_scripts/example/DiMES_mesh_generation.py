#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:54:25 2024

@author: cappellil
"""
#%%

# set working directory equal to gmsh_automated_scripts directory
import sys
import os

# Add the path to the script directory
script_path = os.path.abspath(os.environ['PWD'])
sys.path.append("../")

from gmsh_automated_scripts.utils2 import generate_dimes_mesh 
# Initialize GMSH


#%% 

""" INPUT dots in structure """

Names = ["Dot_1","Dot_2","Dot_3"]

input_dict = {
    Names[0]: {
        "shape": "circle",
        "x": 0,
        "y": 0.75,
        "radius": 0.1  # For cylinder
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
        "radius": 0.1
    }
    # Add more dots as needed
}

generate_dimes_mesh(input_dict, theta_dimes=25, msh_dim=2)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:54:25 2024

@author: cappellil
"""
# %%
# it is recommended to use this script in a conda/mamba environment, .e.g. with:
 # mamba create -n GITRMProject python=3.11 matplotlib numpy spyder=5.5`
 # mamba activate GITRMProject`
 # cd GITRMProject/gmsh_automated_scripts
 # poetry install
# In case of liblglui missin pleasae install sudo apt-get -y install libglu1

from gmsh_automated_scripts import generate_dimes_mesh


""" INPUT dots in structure """

input_dict = {
    # "Sample_holder": {
    #     "shape": "circle",
    #     "x": 0.,
    #     "y": 0.,
    #     "radius": 1.55, 
    #     "theta_dot": 0.
    # },
    "Dot_1": {
        "shape": "circle",
        "x": -1,
        "y": 0,
        "radius": 0.3,  # For cylinder
        "theta_dot": 0.
    },
    "Rectangle_2": {
        "shape": "rectangle",
        "x":  0,
        "y":  -0.6,
        "width": 1.2,
        "height": 1.2,
        "theta_dot": 0.
    }
    # "Dot_3": {
    #     "shape": "circle",
    #     "x": 0,
    #     "y": -1.5,
    #     "radius": 0.5,
    #     "theta_dot": 10
    # }
    # Add more dots as needed
}

path = "/home/cappellil/GITRmProject/mesh_repo/"

generate_dimes_mesh(input_dict, msh_dim=2, ay=-1, r_dimes = 2.5, theta_dimes = 5,
                    filename = path + "DiMES_3D.msh", save_msh = True)

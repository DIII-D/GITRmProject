#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:54:25 2024

@author: cappellil
"""
# %%

from gmsh_automated_scripts import generate_dimes_mesh


""" INPUT dots in structure """

input_dict = {
    "Dot_1": {
        "shape": "circle",
        "x": 0,
        "y": 0.75,
        "radius": 0.1  # For cylinder
    },
    "Rectangle_2": {
        "shape": "rectangle",
        "x": -0.5,
        "y": -0.25,
        "width": 1,
        "height": 0.5
    },
    "Dot_3": {
        "shape": "circle",
        "x": 0,
        "y": -0.75,
        "radius": 0.1
    }
    # Add more dots as needed
}

generate_dimes_mesh(input_dict, theta_dimes=25, msh_dim=2)

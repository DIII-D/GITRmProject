#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from gmsh_automated_scripts import generate_dimes_mesh

scale = 0.01 # from meters to cm

input_dict = {
    "samples_holder": {
        "shape": "circle",
        "x": 0.,
        "y": 0.,
        "radius": 1.5,
        "deposits": {
                "std_button": {
                    "shape": "circle",
                    "x": -1.,
                    "y": 0.,
                    "radius": 0.3,
                },
                "WEST_sample": {
                    "shape": "rectangle",
                    "x": 0.,
                    "y": -0.6,
                    "width": 1.2,
                    "height": 1.2,
                    "theta_dot": 0.
                },
            }
        },
    # W_coating sits between r=2 and r=2.5, which is outside samples_holder (r=1.5),
    # so it is a direct deposit on DiMES (r_dimes=2.5), not on samples_holder.
    "W_coating": {
        "shape":     "annulus",
        "x":         0.,
        "y":         0.,
        "r_inner":   2.,
        "r_outer":   2.45,   # must be < r_dimes (2.5); equal causes degenerate surface
        "phi_start": 90.,
        "phi_end":   270.,
    },
}

generate_dimes_mesh(input_dict, msh_dim=2, ay=-1, r_dimes=2.5, theta_dimes=5,
                    save_msh=False, GUI_msh=True, 
                    save_output = True,
                    save_ply = True, ply_normals=True, ply_ascii=True, 
                    save_npz = False, scale = scale)

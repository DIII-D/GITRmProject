#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: annulus shapes on DiMES.

Demonstrates all three annulus variants:
  - full annulus  (angle = 360, default)
  - partial sector (angle < 360)
  - annulus with a circular deposit inside

Run:
    python DiMES_mesh_annulus.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gmsh_automated_scripts import generate_dimes_mesh

input_dict = {

    # Full annulus: a ring centred on the DiMES axis
    "ring": {
        "shape":   "annulus",
        "x":       0.,
        "y":       0.,
        "r_inner": 0.5,
        "r_outer": 1.0,
    },

    # 90-degree sector specified via angle
    "sector_90": {
        "shape":     "annulus",
        "x":         0.,
        "y":         0.,
        "r_inner":   1.2,
        "r_outer":   1.8,
        "phi_start": 0.,
        "angle":     90.,
    },

    # sector specified via phi_start / phi_end (equivalent to angle=180)
    "sector_180": {
        "shape":     "annulus",
        "x":         0.,
        "y":         0.,
        "r_inner":   1.2,
        "r_outer":   1.8,
        "phi_start": 180.,
        "phi_end":   360.,
    },

    # Annulus with a circular deposit punched through it
    "ring_with_hole": {
        "shape":   "annulus",
        "x":       0.,
        "y":       0.,
        "r_inner": 2.0,
        "r_outer": 2.3,
        "deposits": {
            "notch": {
                "shape":  "circle",
                "x":      2.15,   # sits inside the ring band
                "y":      0.,
                "radius": 0.07,
            }
        },
    },
}

generate_dimes_mesh(
    input_dict,
    msh_dim=2,
    r_dimes=2.5,
    theta_dimes=0,
    save_msh=False,
    save_ply=False,
    save_npz = False,
    GUI_msh=True
)

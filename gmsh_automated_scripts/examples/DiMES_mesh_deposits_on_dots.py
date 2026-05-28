#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example: deposits on dots (substrate/deposit hierarchy)

Demonstrates nesting material deposits inside a larger coating on DiMES.
The parent coating acts as a substrate and its surface is perforated where
each child deposit sits, exactly as the DiMES surface is perforated for
top-level coatings.

Geometry (all coordinates in absolute global frame, all coplanar at z=z_top_dimes):

  DiMES head  r=2.5, center (0, 0)
  └─ Substrate      r=1.2, center (0, 0)        -- large coating on DiMES
       ├─ Inner_deposit   r=0.40, center (0, +0.3)  -- sits inside Substrate
       │    └─ Core_deposit  r=0.15, center (0, +0.3) -- sits inside Inner_deposit
       └─ Side_deposit    r=0.30, center (0, -0.6)  -- sits inside Substrate

Validity checks (each child must fit entirely within its parent):
  Substrate inside DiMES      : max dist 0+1.2 = 1.2 < 2.5  ok
  Inner_deposit inside Substrate: max dist 0.3+0.4 = 0.7 < 1.2  ok
  Core_deposit inside Inner_deposit: 0.15 < 0.4  ok
  Side_deposit inside Substrate : max dist 0.6+0.3 = 0.9 < 1.2  ok
  Inner_deposit vs Side_deposit : centre dist 0.9 > 0.4+0.3 = 0.7  ok (no overlap)
"""

from gmsh_automated_scripts import generate_dimes_mesh

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
        "r_outer":   2.4,   # must be < r_dimes (2.5); equal causes degenerate surface
        "phi_start": 90.,
        "phi_end":   270.,
    },
}

generate_dimes_mesh(input_dict, msh_dim=2, ay=-1, r_dimes=2.5, theta_dimes=5,
                    save_msh=False, GUI_msh=True, 
                    save_ply = True, ply_normals=True, ply_ascii=True, 
                    save_npz = False)

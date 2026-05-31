#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 30 2026

@author: cappellil
"""
# %%

from gmsh_automated_scripts.data_structures import *
from gmsh_automated_scripts import build_dimes_domain
import numpy as np
import os

scale = 0.01

WEST_sample = Rectangle(1.2, 1.2, [0, -0.6, 0.]).scale_position(scale=scale).scale_size(scale=scale)
W_std_button = Disk(0.3, [-1, 0., 0.]).scale_position(scale=scale).scale_size(scale=scale)
W_coating = Annulus(2.1, 2.48, (90, 270), [0., 0., 0.]).scale_position(scale=scale).scale_size(scale=scale)

shapes = [WEST_sample, W_std_button, W_coating]
labels = ["WEST_sample", "W_std_button", "W_coating"]

base_path = os.path.expanduser("~")

objects = build_dimes_domain(
    shapes=shapes,
    labels=labels,
    rot_axis=[0., -1., 0.],
    angle=np.deg2rad(5),
    GUI_msh=True,
    GUI_geo=False,
    save_output=True,
    save_ply=True,
    base_path = base_path,
    scale = scale
)

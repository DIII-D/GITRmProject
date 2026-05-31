#!/usr/bin/env python3
"""
Three AngledSample objects on a flat DiMES head (dimes tilt angle = 0).
Samples are placed symmetrically along the x-axis, all within the 2.5 cm DiMES disk.
Run from the repo root:  python -m examples.angled_samples_flat_dimes
"""

import math
import numpy as np
from gmsh_automated_scripts.data_structures import AngledSample, MeshConfig
from gmsh_automated_scripts import build_dimes_domain

scale = 0.01  # cm -> m

# Three angled samples: centre, right, left
sample_r      = 0.3   # cm
sample_height = 0.05  # cm
sample_angle  = math.radians(10)
sample_z_cut  = 0.015 # cm  (must be < 2 * r * tan(angle) ≈ 0.106 cm)

x_tl, y_tl = -math.cos(math.radians(30)), math.sin(math.radians(30))
x_tr, y_tr = math.cos(math.radians(30)), math.sin(math.radians(30))

x_bl, y_bl = -math.cos(math.radians(30)), -math.sin(math.radians(30))
x_br, y_br = math.cos(math.radians(30)), -math.sin(math.radians(30))

shapes = [
    AngledSample(center=[ 0., 0., 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[ 0., 1., 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[ 0., -1., 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[x_tl, y_tl, 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[x_tr, y_tr, 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[x_br, y_br, 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),
    AngledSample(center=[x_bl, y_bl, 0.], r=sample_r, height=sample_height, angle=sample_angle, z_cut=sample_z_cut).scale_position(scale).scale_size(scale),

]
labels = ["sample_centre", "sample_top", "sample_bottom", "sample_top_left", "sample_top_right",
           "sample_bottom_left", "sample_bottom_right"]

mesh = MeshConfig().scale_parameters(scale)

objects = build_dimes_domain(
    shapes=shapes,
    labels=labels,
    angle=0.,           # flat DiMES head
    mesh=mesh,
    GUI_geo=True,
    GUI_msh=True,
    scale=scale,
)

for obj in objects:
    print(f"{obj.label:20s}  surface_tag={obj.surface_tag}")

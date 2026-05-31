#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:06:17 2024

@author: cappellil
"""

from .data_structures import *
from .helpers import *
from .export import *
import gmsh
import math
import numpy as np


def add_sample(sample):
    """Returns (hole_loop, sample_surface). Caller punches hole_loop into the parent
       and rotates the surfaces; boundary curves follow their surface automatically."""
    if isinstance(sample, Disk):
        x, y, z = sample.center
        loop = gmsh.model.occ.addCurveLoop([gmsh.model.occ.addCircle(x, y, z, sample.r)])
        return loop, gmsh.model.occ.addPlaneSurface([loop])

    elif isinstance(sample, Rectangle):
        x, y, z = sample.ll
        *_, loop = rectangle_loop(x, y, z, sample.width, sample.height)
        return loop, gmsh.model.occ.addPlaneSurface([loop])

    elif isinstance(sample, Annulus):
        x, y, z = sample.center
        phi_start, phi_end = sample.angular_sector
        dimes_outer_loop, inner_loops = annulus_loop(
            x, y, z, sample.r_inner, sample.r_outer, phi_start, phi_end - phi_start
        )
        return dimes_outer_loop, gmsh.model.occ.addPlaneSurface([dimes_outer_loop] + inner_loops)

    else:
        raise TypeError(f"Sample type {type(sample)} is not allowed.")


def build_dimes_domain(
    shapes: list,
    labels: list[str],
    rot_axis: list[float] = None,
    angle: float = 0.,
    dimes_shape: Disk = None,
    L_tile: float = 6.,
    mesh: MeshConfig = None,
    GUI_geo: bool = False,
    GUI_msh: bool = True,
    filename: str = "dimes.msh",
    save_msh: bool = False,
    save_output: bool = False,
    save_npz: bool = None,
    save_ply: bool = False,
    ply_normals: bool = True,
    ply_ascii: bool = True,
    base_path: str = None,
    scale: float = 1.
) -> list[Object2D]:
    """Build, mesh, and return all domain components with populated surface tags."""
    if rot_axis is None:
        rot_axis = [0., -1., 0.]
    if dimes_shape is None:
        dimes_shape = Disk(2.5, [0., 0., 0.]).scale_position(scale=scale).scale_size(scale=scale)

    gmsh.initialize()

    ax, ay, az = rot_axis

    # Tile
    tile_obj = Object2D(Rectangle(L_tile, L_tile, [-L_tile / 2, -L_tile / 2, 0.]).scale_position(scale=scale).scale_size(scale=scale),
    label = "tile")
    *_, tile_loop = rectangle_loop(*tile_obj.shape.ll, tile_obj.shape.width, tile_obj.shape.height)

    # DiMES top-edge ellipse (projects the tilted disk footprint onto z=0)
    curve = truncated_cylinder_Ellipse(*dimes_shape.center, dimes_shape.r, ax, ay, angle)
    dimes_outer_loop = gmsh.model.occ.addCurveLoop([curve])

    tile_obj.surface_tag = gmsh.model.occ.addPlaneSurface([tile_loop, dimes_outer_loop])

    # Sample holes
    inner_loops, samples = [], []
    for shape, label in zip(shapes, labels):
        loop, surface = add_sample(shape)
        inner_loops.append(loop)
        obj = Object2D(shape=shape, label=label)
        obj.surface_tag = surface
        samples.append(obj)

    # DiMES top surface with sample holes punched in
    dimes_top_obj = Object2D(shape=dimes_shape, label="dimes_top")
    dimes_top_obj.surface_tag = gmsh.model.occ.addPlaneSurface([dimes_outer_loop] + inner_loops)

    # Rotate all DiMES-plane surfaces together
    if angle != 0.:
        cx, cy, cz = dimes_shape.center
        all_surfaces = [(2, dimes_top_obj.surface_tag)] + [(2, s.surface_tag) for s in samples]
        gmsh.model.occ.rotate(all_surfaces, cx, cy, cz, ax, ay, az, angle)
        gmsh.model.occ.synchronize()
        _, _, z_min, _, _, _ = gmsh.model.getBoundingBox(2, dimes_top_obj.surface_tag)
        gmsh.model.occ.translate(all_surfaces, 0, 0, cz - z_min)
        gmsh.model.occ.synchronize()

    objects = [tile_obj, dimes_top_obj] + samples

    # DiMES side surface — only exists when the head is tilted
    if angle != 0.:
        dimes_base_loop = gmsh.model.occ.addCurveLoop(
            [gmsh.model.occ.addCircle(*dimes_shape.center, dimes_shape.r)]
        )
        result = gmsh.model.occ.addThruSections(
            [dimes_base_loop, dimes_outer_loop], makeSolid=False
        )
        gmsh.model.occ.synchronize()
        dimes_side_obj = Object2D(shape=dimes_shape, label="dimes_side")
        dimes_side_obj.surface_tag = result[0][1]
        objects.append(dimes_side_obj)

    component_surfaces = {
        obj.label: [obj.surface_tag] if isinstance(obj.surface_tag, int) else obj.surface_tag
        for obj in objects
    }

    make_dimes_mesh(
        mesh=mesh,
        component_surfaces=component_surfaces,
        GUI_geo=GUI_geo,
        GUI_msh=GUI_msh,
        filename=filename,
        save_msh=save_msh,
        save_output=save_output,
        save_npz=save_npz,
        save_ply=save_ply,
        ply_normals=ply_normals,
        ply_ascii=ply_ascii,
        base_path=base_path,
        scale=scale
    )

    return objects
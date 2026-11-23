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


def make_flush_sample(sample):
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


def connect_surfaces(objects: list[Object2D]):
    """Fragment the objects' surfaces together so touching surfaces share curves and
       points (conformal mesh). Fragment may renumber or split surfaces, so each
       object's surface_tag is remapped in place from the fragment output map."""
    if len(objects) < 2:
        return
    tags_per_obj = [[obj.surface_tag] if isinstance(obj.surface_tag, int) else list(obj.surface_tag)
                    for obj in objects]
    dimtags = [(2, t) for tags in tags_per_obj for t in tags]

    _, out_map = gmsh.model.occ.fragment(dimtags[:1], dimtags[1:])
    gmsh.model.occ.synchronize()

    # out_map[i] lists the surfaces that replaced dimtags[i] (same order)
    k = 0
    for obj, tags in zip(objects, tags_per_obj):
        new_tags = [t for i in range(k, k + len(tags)) for _, t in out_map[i]]
        k += len(tags)
        obj.surface_tag = new_tags[0] if isinstance(obj.surface_tag, int) and len(new_tags) == 1 else new_tags


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
    save_all: bool = False,
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

    # DiMES top-edge ellipse: once rotated, it lies on the vertical cylinder of radius r
    curve = truncated_cylinder_Ellipse(*dimes_shape.center, dimes_shape.r, ax, ay, angle)
    dimes_outer_loop = gmsh.model.occ.addCurveLoop([curve])

    # DiMES footprint on the tile is the cylinder base circle (same as the top when flat).
    # The same curve is reused by the side wall, so tile and side share it.
    if angle != 0.:
        base_curve = truncated_cylinder_base(*dimes_shape.center, dimes_shape.r, ax, ay, angle)
        dimes_base_loop = gmsh.model.occ.addCurveLoop([base_curve])
    else:
        dimes_base_loop = dimes_outer_loop

    tile_obj.surface_tag = gmsh.model.occ.addPlaneSurface([tile_loop, dimes_base_loop])

    # Sample holes
    inner_loops, samples = [], []
    for shape, label in zip(shapes, labels):

        obj = Object2D(shape=shape, label=label)
        
        if isinstance(shape, AngledSample):
            loop, surface_tags = make_angled_sample(obj)
            obj.surface_tag = surface_tags
        else:
            loop, surface = make_flush_sample(shape)
            obj.surface_tag = surface

        inner_loops.append(loop)
        samples.append(obj)

    # DiMES top surface with sample holes punched in
    dimes_top_obj = Object2D(shape=dimes_shape, label="dimes_top")
    dimes_top_obj.surface_tag = gmsh.model.occ.addPlaneSurface([dimes_outer_loop] + inner_loops)

    # Rotate all DiMES-plane surfaces together
    if angle != 0.:
        cx, cy, cz = dimes_shape.center
        # lowest point of the rotated ellipse drops by (r/cos)*sin = r*tan; lift it back to cz.
        # Exact value instead of getBoundingBox, which OCC computes loosely.
        dz = dimes_shape.r * abs(math.tan(angle))
        all_surfaces = [(2, t) for o in [dimes_top_obj] + samples
                        for t in ([o.surface_tag] if isinstance(o.surface_tag, int) else o.surface_tag)]
        gmsh.model.occ.rotate(all_surfaces, cx, cy, cz, ax, ay, az, angle)
        gmsh.model.occ.translate(all_surfaces, 0, 0, dz)
        gmsh.model.occ.synchronize()

    # Connect DiMES top and samples (e.g. angled samples whose hole loop is a separate curve).
    # When tilted, the tile and side wall are left out: they share curves by construction, and
    # fragment would split the rim curves at the tangent touching point and break the ruled side.
    connect_surfaces(([tile_obj] if angle == 0. else []) + [dimes_top_obj] + samples)

    objects = [tile_obj, dimes_top_obj] + samples

    # DiMES side surface — only exists when the head is tilted
    if angle != 0.:
        top_curve = find_rim_curve(dimes_top_obj.surface_tag, cx, cy, dimes_shape.r)
        dimes_side_obj = Object2D(shape=dimes_shape, label="dimes_side")
        dimes_side_obj.surface_tag = truncated_cylinder_side(base_curve, top_curve)
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
        save_all=save_all,
        save_npz=save_npz,
        save_ply=save_ply,
        ply_normals=ply_normals,
        ply_ascii=ply_ascii,
        base_path=base_path,
        scale=scale
    )

    return objects
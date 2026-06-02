#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mesh export: gmsh meshes once and writes a single .msh; PLY/NPZ are derived
from that file per physical group. The .msh already stores one conformal global
node block, so per-group extraction shares seam nodes automatically -- no
includeBoundary duplication and no coordinate welding to tune."""

import os
import tempfile
import numpy as np
import meshio
import gmsh
from .data_structures import MeshConfig


def compute_vertex_normals(nodes, triangles):
    """Per-vertex normals as the area-weighted average of adjacent face normals.

    Parameters
    ----------
    nodes : (N, 3) float array
    triangles : (M, 3) int array

    Returns
    -------
    normals : (N, 3) float64 array, unit-length per vertex.
    """
    v0 = nodes[triangles[:, 0]]
    v1 = nodes[triangles[:, 1]]
    v2 = nodes[triangles[:, 2]]
    face_normals = np.cross(v1 - v0, v2 - v0)

    vertex_normals = np.zeros_like(nodes, dtype=np.float64)
    for i in range(3):
        np.add.at(vertex_normals, triangles[:, i], face_normals)

    norms = np.linalg.norm(vertex_normals, axis=1, keepdims=True)

    degenerate = (norms < np.finfo(np.float64).eps).ravel()
    if np.any(degenerate):
        import warnings
        warnings.warn(
            f"{degenerate.sum()} vertex/vertices with near-zero normals detected "
            f"(indices: {np.where(degenerate)[0].tolist()}). "
            "These are likely boundary or unreferenced vertices. "
            "Falling back to (0, 0, 1).",
            UserWarning,
            stacklevel=2,
        )
        vertex_normals[degenerate] = [0.0, 0.0, 1.0]
        norms[degenerate] = 1.0

    return vertex_normals / norms


def write_ply(path, nodes, triangles, normals=None, ascii_format=False):
    """Write a PLY file, binary little-endian by default."""
    verts = nodes.astype(np.float32)
    faces = triangles.astype(np.uint32)
    has_normals = normals is not None
    norms = normals.astype(np.float32) if has_normals else None

    normal_props = (
        "property float nx\n"
        "property float ny\n"
        "property float nz\n"
    ) if has_normals else ""

    fmt_line = "format ascii 1.0" if ascii_format else "format binary_little_endian 1.0"

    header = (
        "ply\n"
        f"{fmt_line}\n"
        f"element vertex {len(verts)}\n"
        "property float x\n"
        "property float y\n"
        "property float z\n"
        f"{normal_props}"
        f"element face {len(faces)}\n"
        "property list uchar uint vertex_indices\n"
        "end_header\n"
    )

    if ascii_format:
        with open(path, "w") as f:
            f.write(header)
            for i, v in enumerate(verts):
                line = f"{v[0]} {v[1]} {v[2]}"
                if has_normals:
                    n = norms[i]
                    line += f" {n[0]} {n[1]} {n[2]}"
                f.write(line + "\n")
            for tri in faces:
                f.write(f"3 {tri[0]} {tri[1]} {tri[2]}\n")
    else:
        vert_data = np.hstack([verts, norms]) if has_normals else verts
        n_faces = len(faces)
        face_buf = np.empty(n_faces * 13, dtype=np.uint8)
        face_view = face_buf.reshape(n_faces, 13)
        face_view[:, 0] = 3
        face_view[:, 1:] = faces.view(np.uint8).reshape(n_faces, 12)
        with open(path, "wb") as f:
            f.write(header.encode("ascii"))
            f.write(vert_data.tobytes())
            f.write(face_buf.tobytes())


def _read_surface_groups(msh_path):
    """Read a .msh, return (global_points (N,3), {physical_name: triangles}).

    Triangles are (M, 3) indices into the shared global point array, so nodes
    shared between groups (the conformal seam) carry the same index.
    """
    m = meshio.read(msh_path)
    points = m.points
    tag_to_name = {tag: name for name, (tag, dim) in m.field_data.items() if dim == 2}

    groups = {}
    for cb, phys in zip(m.cells, m.cell_data.get("gmsh:physical", [])):
        if cb.type != "triangle":
            continue
        for tag in np.unique(phys):
            name = tag_to_name.get(int(tag))
            if name is None:
                continue
            groups.setdefault(name, []).append(cb.data[phys == tag])

    return points, {name: np.vstack(blocks) for name, blocks in groups.items()}


def _compact(points, tris_global):
    """Reduce a global-indexed triangle set to a compact local (nodes, triangles)
    pair containing only the referenced (deduplicated) vertices."""
    used = np.unique(tris_global)
    remap = {int(t): i for i, t in enumerate(used)}
    nodes = points[used]
    triangles = np.vectorize(remap.__getitem__)(tris_global).astype(np.int64)
    return nodes, triangles


def save_component_meshes(msh_path, base_path, save_npz=True, save_ply=False,
                          ply_normals=True, ply_ascii=False,
                          weld_dimes=True, dimes_groups=("dimes_top", "dimes_side")):
    """Derive one .npz (and optionally .ply) per component from a .msh.

    Each physical group is read against the shared global node block and then
    compacted, so a component's surfaces (e.g. an angled sample's top + side,
    which live under one physical group) come out as a single watertight mesh
    with the seam already shared -- no welding step needed.

    Parameters
    ----------
    msh_path : str
        Path to the .msh written by gmsh.
    weld_dimes : bool
        If True, the `dimes_groups` (top + side, two separate physical groups)
        are fused into one 'dimes' output with the rim deduplicated. If False
        they are written as separate components.
    dimes_groups : tuple[str, str]
        Names of the DiMES top and side physical groups to fuse.
    """
    points, groups = _read_surface_groups(msh_path)

    # decide the output grouping
    outputs = {}
    if weld_dimes and all(g in groups for g in dimes_groups):
        outputs["dimes"] = np.vstack([groups.pop(g) for g in dimes_groups])
    outputs.update(groups)                      # remaining groups as-is

    for name, tris_global in outputs.items():
        if len(tris_global) == 0:
            continue
        nodes_out, tris_out = _compact(points, tris_global)

        saved = []
        if save_npz:
            npz_path = f"{base_path}/{name}.npz"
            np.savez(npz_path, nodes=nodes_out, triangles=tris_out)
            saved.append(npz_path)
        if save_ply:
            ply_path = f"{base_path}/{name}.ply"
            normals_out = compute_vertex_normals(nodes_out, tris_out) if ply_normals else None
            write_ply(ply_path, nodes_out, tris_out, normals=normals_out, ascii_format=ply_ascii)
            saved.append(ply_path)
        if saved:
            print(f"  saved {name:30s}  {len(nodes_out):6d} nodes  "
                  f"{len(tris_out):6d} triangles  ->  {',  '.join(saved)}")


def make_dimes_mesh(mesh: MeshConfig = None, filename="test.msh", save_msh=False,
                    GUI_geo=False, GUI_msh=True, component_surfaces=None,
                    save_output=False, save_npz=None, save_ply=False,
                    ply_normals=True, ply_ascii=True, scale=1., base_path: str = None,
                    weld_dimes=True):
    """Finalize geometry, mesh it per `mesh`, label components, show GUIs, export.

    The mesh is written to a single .msh and PLY/NPZ are derived from it per
    physical group (see save_component_meshes). `weld_dimes` controls whether
    the DiMES top and side groups are fused into one component.
    """
    if mesh is None:
        mesh = MeshConfig().scale_parameters(scale=scale)

    gmsh.model.occ.synchronize()

    if GUI_geo:
        gmsh.fltk.run()

    # Label each component as a named Physical Surface (2-D meshes only)
    if mesh.dim == 2 and component_surfaces is not None:
        for comp_name, surf_tags in component_surfaces.items():
            pg = gmsh.model.addPhysicalGroup(2, surf_tags)
            gmsh.model.setPhysicalName(2, pg, comp_name)

    mesh.generate()

    if GUI_msh:
        gmsh.fltk.run()

    # save_npz=None means "follow save_msh"; explicit True/False overrides
    _save_npz = save_msh if save_npz is None else save_npz
    want_export = (mesh.dim == 2 and (_save_npz or save_ply)
                   and save_output and component_surfaces is not None)

    # Write the .msh once: keep it at `filename` if requested, else to a temp
    # file used only to derive the per-component outputs and then removed.
    msh_path = None
    if save_msh:
        gmsh.write(filename)
        msh_path = filename
    elif want_export:
        fd, msh_path = tempfile.mkstemp(suffix=".msh")
        os.close(fd)
        gmsh.write(msh_path)

    gmsh.finalize()

    if want_export:
        _base_path = base_path if base_path is not None else os.path.splitext(filename)[0]
        save_component_meshes(msh_path, _base_path,
                              save_npz=_save_npz, save_ply=save_ply,
                              ply_normals=ply_normals, ply_ascii=ply_ascii,
                              weld_dimes=weld_dimes)
        if not save_msh:
            os.remove(msh_path)
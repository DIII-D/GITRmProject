#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gmsh
import numpy as np
import os
from .data_structures import MeshConfig


def extract_surface_mesh(surf_tag):
    """Return (nodes_xyz, triangles) for one OCC surface tag.

    Nodes is (N, 3) float64.  Triangles is (M, 3) int64 of local 0-based indices.
    Returns empty arrays if the surface carries no 2-D elements.
    """
    node_tags, coords, _ = gmsh.model.mesh.getNodes(
        dim=2, tag=surf_tag, includeBoundary=True)

    if len(node_tags) == 0:
        return np.empty((0, 3)), np.empty((0, 3), dtype=np.int64)

    nodes_xyz = coords.reshape(-1, 3)
    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}

    elem_types, _, node_conns = gmsh.model.mesh.getElements(dim=2, tag=surf_tag)

    tri_blocks = []
    for etype, conn in zip(elem_types, node_conns):
        if etype == 2:  # 3-node triangle
            tri_blocks.append(
                np.array([tag_to_idx[int(t)] for t in conn],
                         dtype=np.int64).reshape(-1, 3))

    triangles = np.vstack(tri_blocks) if tri_blocks else np.empty((0, 3), dtype=np.int64)
    return nodes_xyz, triangles


def compute_vertex_normals(nodes, triangles):
    """Compute per-vertex normals as the area-weighted average of adjacent face normals.

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

    # Identify degenerate vertices and warn
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
    """Write a PLY file, binary little-endian by default.

    Parameters
    ----------
    path : str
        Output file path (should end with .ply).
    nodes : (N, 3) float array
        Vertex coordinates.
    triangles : (M, 3) int array
        Triangle vertex indices (0-based).
    normals : (N, 3) float array or None
        Per-vertex normals.  If None, normals are omitted.
    ascii_format : bool
        Write ASCII format instead of binary little-endian.
    """
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

        # face buffer: 1 uint8 (count=3) + 3 uint32s = 13 bytes per face
        n_faces = len(faces)
        face_buf = np.empty(n_faces * 13, dtype=np.uint8)
        face_view = face_buf.reshape(n_faces, 13)
        face_view[:, 0] = 3
        face_view[:, 1:] = faces.view(np.uint8).reshape(n_faces, 12)

        with open(path, "wb") as f:
            f.write(header.encode("ascii"))
            f.write(vert_data.tobytes())
            f.write(face_buf.tobytes())


def save_component_meshes(component_surfaces, base_path, save_npz=True, save_ply=False,
                          ply_normals=True, ply_ascii=False):
    """Save one .npz (and optionally one .ply) per component.

    Parameters
    ----------
    component_surfaces : dict
        {name: [surface_tag, ...]} as returned by make_dimes_geom.
    base_path : str
        Path prefix without extension, e.g. "/path/to/DiMES_3D".
        Output files are written as  <base_path>/<name>.npz / .ply.
    save_npz : bool
        Write a numpy .npz file for each component (default True).
    save_ply : bool
        Write a PLY file for each component (default False).
    ply_normals : bool
        Include per-vertex normals in PLY output (default True).
    ply_ascii : bool
        Write PLY in ASCII format instead of binary little-endian (default False).
    """
    for name, surf_tags in component_surfaces.items():
        all_nodes = []
        all_tris = []
        offset = 0

        for tag in surf_tags:
            try:
                nodes, tris = extract_surface_mesh(tag)
            except Exception:
                continue
            if len(nodes) == 0:
                continue
            all_nodes.append(nodes)
            all_tris.append(tris + offset)
            offset += len(nodes)

        if not all_nodes:
            continue

        nodes_out = np.vstack(all_nodes)
        tris_out = (np.vstack(all_tris) if all_tris
                    else np.empty((0, 3), dtype=np.int64))

        saved = []

        if save_npz:
            npz_path = f"{base_path}/{name}.npz"
            np.savez(npz_path, nodes=nodes_out, triangles=tris_out)
            saved.append(npz_path)

        if save_ply and len(tris_out) > 0:
            ply_path = f"{base_path}/{name}.ply"
            normals_out = compute_vertex_normals(nodes_out, tris_out) if ply_normals else None
            write_ply(ply_path, nodes_out, tris_out, normals=normals_out, ascii_format=ply_ascii)
            saved.append(ply_path)

        if saved:
            print(f"  saved {name:30s}  {len(nodes_out):6d} nodes  "
                  f"{len(tris_out):6d} triangles  →  {',  '.join(saved)}")

def make_dimes_mesh(mesh: MeshConfig = None, filename="test.msh", save_msh=False,
                    GUI_geo=False, GUI_msh=True, component_surfaces=None,
                    save_output=False, save_npz=None, save_ply=False,
                    ply_normals=True, ply_ascii=True, scale=1., base_path:str=None):
    """Finalize geometry, mesh it per `mesh`, label components, show GUIs, export.

    Mesh options now live in `mesh` (a MeshConfig). `mesh.dim` replaces the old
    `msh_dim` argument. Pass mesh=None to use defaults that match the original.
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

    # Apply mesh options and generate
    mesh.generate()

    if GUI_msh:
        gmsh.fltk.run()

    if save_msh:
        gmsh.write(filename)

    # save_npz=None means "follow save_msh"; explicit True/False overrides
    _save_npz = save_msh if save_npz is None else save_npz

    # Save per-component files for 2-D meshes (independent of save_msh)
    if mesh.dim == 2 and (_save_npz or save_ply) and save_output and component_surfaces is not None:
        _base_path = base_path if base_path is not None else os.path.splitext(filename)[0]
        save_component_meshes(component_surfaces, _base_path, save_npz=_save_npz,
                              save_ply=save_ply, ply_normals=ply_normals, ply_ascii=ply_ascii)

    # Finalize GMSH
    gmsh.finalize()


# def save_component_meshes(component_surfaces, base_path, save_npz=True, save_ply=False,
#                           ply_normals=True, ply_ascii=False):
#     """Save one .npz (and optionally one .ply) per component.

#     Parameters
#     ----------
#     component_surfaces : dict
#         {name: [surface_tag, ...]} as returned by make_dimes_geom.
#     base_path : str
#         Path prefix without extension, e.g. "/path/to/DiMES_3D".
#         Output files are written as  <base_path>_<name>.npz / .ply.
#     save_npz : bool
#         Write a numpy .npz file for each component (default True).
#     save_ply : bool
#         Write a PLY file for each component (default False).
#     ply_normals : bool
#         Include per-vertex normals in PLY output (default True).
#     ply_ascii : bool
#         Write PLY in ASCII format instead of binary little-endian (default False).
#     """
#     for name, surf_tags in component_surfaces.items():
#         all_nodes = []
#         all_tris = []
#         offset = 0

#         for tag in surf_tags:
#             try:
#                 nodes, tris = extract_surface_mesh(tag)
#             except Exception:
#                 continue
#             if len(nodes) == 0:
#                 continue
#             all_nodes.append(nodes)
#             all_tris.append(tris + offset)
#             offset += len(nodes)

#         if not all_nodes:
#             continue

#         nodes_out = np.vstack(all_nodes)
#         tris_out = (np.vstack(all_tris) if all_tris
#                     else np.empty((0, 3), dtype=np.int64))

#         saved = []

#         if save_npz:
#             npz_path = f"{base_path}_{name}.npz"
#             np.savez(npz_path, nodes=nodes_out, triangles=tris_out)
#             saved.append(npz_path)

#         if save_ply and len(tris_out) > 0:
#             ply_path = f"{base_path}_{name}.ply"
#             normals_out = compute_vertex_normals(nodes_out, tris_out) if ply_normals else None
#             write_ply(ply_path, nodes_out, tris_out, normals=normals_out, ascii_format=ply_ascii)
#             saved.append(ply_path)

#         if saved:
#             print(f"  saved {name:30s}  {len(nodes_out):6d} nodes  "
#                   f"{len(tris_out):6d} triangles  →  {',  '.join(saved)}")

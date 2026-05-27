#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gmsh
import numpy as np


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


def write_ply(path, nodes, triangles):
    """Write a binary little-endian PLY file (no extra dependencies required).

    Parameters
    ----------
    path : str
        Output file path (should end with .ply).
    nodes : (N, 3) float array
        Vertex coordinates.
    triangles : (M, 3) int array
        Triangle vertex indices (0-based).
    """
    verts = nodes.astype(np.float32)
    faces = triangles.astype(np.int32)

    header = (
        "ply\n"
        "format binary_little_endian 1.0\n"
        f"element vertex {len(verts)}\n"
        "property float x\n"
        "property float y\n"
        "property float z\n"
        f"element face {len(faces)}\n"
        "property list uchar int vertex_indices\n"
        "end_header\n"
    ).encode("ascii")

    # face data: each row is  [3, i0, i1, i2]  (count byte + 3 int32s)
    counts = np.full((len(faces), 1), 3, dtype=np.uint8)
    face_data = np.hstack([counts.view(np.uint8),
                           faces.view(np.uint8).reshape(len(faces), 12)])

    with open(path, "wb") as f:
        f.write(header)
        f.write(verts.tobytes())
        f.write(face_data.tobytes())


def save_component_meshes(component_surfaces, base_path, save_ply=False):
    """Save one .npz (and optionally one .ply) per component.

    Parameters
    ----------
    component_surfaces : dict
        {name: [surface_tag, ...]} as returned by make_dimes_geom.
    base_path : str
        Path prefix without extension, e.g. "/path/to/DiMES_3D".
        Output files are written as  <base_path>_<name>.npz / .ply.
    save_ply : bool
        Also write a binary PLY file for each component (default False).
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

        npz_path = f"{base_path}_{name}.npz"
        np.savez(npz_path, nodes=nodes_out, triangles=tris_out)

        extra = ""
        if save_ply and len(tris_out) > 0:
            ply_path = f"{base_path}_{name}.ply"
            write_ply(ply_path, nodes_out, tris_out)
            extra = f"  +  {ply_path}"

        print(f"  saved {name:30s}  {len(nodes_out):6d} nodes  "
              f"{len(tris_out):6d} triangles  →  {npz_path}{extra}")

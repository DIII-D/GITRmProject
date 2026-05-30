import os
from dataclasses import dataclass, field
from typing import Optional
import gmsh


@dataclass
class MeshConfig:
    """All mesh-generation settings in one validated place.

    Apply with .apply() (before generate) or .generate() to apply + mesh.
    The defaults reproduce the original function's *effective* behaviour, so an
    existing pipeline meshes identically when called with MeshConfig().
    """
    dim: int = 2                          # dimension passed to gmsh.model.mesh.generate()

    # --- element sizing -----------------------------------------------------
    size_min: float = 0.1                # Mesh.MeshSizeMin (prevents slivers in small dots)
    size_max: float = 0.2                 # Mesh.MeshSizeMax

    # --- curvature-based refinement -----------------------------------------
    elements_per_2pi: float = 10.0        # Mesh.MeshSizeFromCurvature (== MinimumElementsPerTwoPi); 0 = off

    # --- algorithm / quality (optional; None = leave gmsh default) ----------
    algorithm_2d: Optional[int] = None    # Mesh.Algorithm    (default 6 = Frontal-Delaunay)
    algorithm_3d: Optional[int] = None    # Mesh.Algorithm3D  (default 1 = Delaunay)
    element_order: int = 1                # Mesh.ElementOrder (2 = quadratic)
    optimize: bool = True                 # Mesh.Optimize

    # --- escape hatch for any other Mesh.* option ---------------------------
    extra_options: dict = field(default_factory=dict)   # e.g. {"Mesh.Smoothing": 5}

    def __post_init__(self):
        if self.dim not in (1, 2, 3):
            raise ValueError(f"dim must be 1, 2 or 3; got {self.dim}")
        if self.size_min <= 0 or self.size_max <= 0:
            raise ValueError("size_min and size_max must be positive")
        if self.size_min > self.size_max:
            raise ValueError(f"size_min ({self.size_min}) > size_max ({self.size_max})")
        if self.elements_per_2pi < 0:
            raise ValueError("elements_per_2pi must be >= 0 (0 disables curvature sizing)")
        if self.element_order not in (1, 2):
            raise ValueError("element_order must be 1 (linear) or 2 (quadratic)")

    def apply(self) -> None:
        """Push the settings onto the current gmsh model (call before generate())."""
        opt = gmsh.option.setNumber
        opt("Mesh.MeshSizeMin", self.size_min)
        opt("Mesh.MeshSizeMax", self.size_max)
        opt("Mesh.MeshSizeFromCurvature", self.elements_per_2pi)
        opt("Mesh.ElementOrder", self.element_order)
        opt("Mesh.Optimize", int(self.optimize))
        if self.algorithm_2d is not None:
            opt("Mesh.Algorithm", self.algorithm_2d)
        if self.algorithm_3d is not None:
            opt("Mesh.Algorithm3D", self.algorithm_3d)
        for name, value in self.extra_options.items():
            opt(name, value)

    def generate(self) -> None:
        """Apply settings and mesh in one call."""
        self.apply()
        gmsh.model.mesh.generate(self.dim)


def make_dimes_mesh(mesh: MeshConfig = None, filename="test.msh", save_msh=False,
                    GUI_geo=False, GUI_msh=True, component_surfaces=None,
                    save_output=False, save_npz=None, save_ply=False,
                    ply_normals=True, ply_ascii=True):
    """Finalize geometry, mesh it per `mesh`, label components, show GUIs, export.

    Mesh options now live in `mesh` (a MeshConfig). `mesh.dim` replaces the old
    `msh_dim` argument. Pass mesh=None to use defaults that match the original.
    """
    if mesh is None:
        mesh = MeshConfig()

    # Remove duplicates (coherence) and synchronize the CAD model
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    # Label each component as a named Physical Surface (2-D meshes only)
    if mesh.dim == 2 and component_surfaces is not None:
        for comp_name, surf_tags in component_surfaces.items():
            pg = gmsh.model.addPhysicalGroup(2, surf_tags)
            gmsh.model.setPhysicalName(2, pg, comp_name)

    if GUI_geo:
        gmsh.fltk.run()

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
        from .export import save_component_meshes
        base_path = os.path.splitext(filename)[0]
        save_component_meshes(component_surfaces, base_path, save_npz=_save_npz,
                              save_ply=save_ply, ply_normals=ply_normals, ply_ascii=ply_ascii)

    # Finalize GMSH
    gmsh.finalize()


# ----------------------------------------------------------------------------
# Usage
# ----------------------------------------------------------------------------
# Default (== original effective settings):
#     make_dimes_mesh(component_surfaces=comps)
#
# Custom mesh:
#     cfg = MeshConfig(dim=2, size_min=0.02, size_max=0.15,
#                      elements_per_2pi=30, element_order=2)
#     make_dimes_mesh(cfg, save_msh=True, filename="dimes.msh",
#                     component_surfaces=comps)
#
# Reuse one config across runs, tweak per call:
#     base = MeshConfig(size_max=0.2)
#     make_dimes_mesh(replace(base, size_max=0.1), ...)   # from dataclasses import replace
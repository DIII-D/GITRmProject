from dataclasses import dataclass, field
from typing import Union, Optional
import gmsh
import math

@dataclass
class Cube:
    L_R: float      # radial length
    L_phi: float    # toroidal length
    L_Z: float      # z-axis length
    center: list[float]  # [R, phi, Z]

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def scale_size(self, scale: float):
        self.L_R *= scale
        self.L_phi *= scale
        self.L_Z *= scale
        return self
    
    def scale_position(self, scale: float):
        self.center = [c * scale for c in self.center]
        return self

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if (self.L_R <= 0.) or (self.L_phi <= 0.) or (self.L_Z <= 0.):
            raise ValueError("One of the cube dimensions is less than or equal to zero")

@dataclass
class Disk:
    r : float # radius length
    center: list[float] # [R, phi, Z]

    def scale_size(self, scale: float) -> None:
        self.r *= scale
        return self
    
    def scale_position(self, scale: float):
        self.center = [c * scale for c in self.center]
        return self

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if self.r <= 0.:
            raise ValueError("Disk radius must be greater than 0")
        
@dataclass
class Rectangle:
    width: float
    height: float
    ll: list[float] # lower left corner on plane parallel to rectangle
    
    def scale_size(self, scale: float) -> None:
        self.width *= scale
        self.height *= scale
        return self
    
    def scale_position(self, scale: float):
        self.ll = [c * scale for c in self.ll]
        return self

@dataclass
class Annulus:
    r_inner: float
    r_outer: float
    angular_sector: tuple
    center: list[float] # [R, phi, Z]

    def scale_size(self, scale: float) -> None:
        self.r_inner *= scale
        self.r_outer *= scale
        return self

    def scale_position(self, scale: float):
        self.center = [c * scale for c in self.center]
        return self

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if (self.r_inner <= 0.):
            raise ValueError("Annulus r_inner must be greater than 0")
        if (self.r_outer <= 0.):
            raise ValueError("Annulus r_outer must be greater than 0")
        if (self.angular_sector[0] >= self.angular_sector[1]):
            raise ValueError("Annulus angular sector final angle must be smaller than the first angle.")


@dataclass
class AngledSample:
    center: list[float] = field(default_factory=lambda: [0., 0., 0.])
    r: float = 0.3
    height: float = 0.01
    angle: float = 10 * math.pi / 180
    z_cut: float = 0.015

    def scale_size(self, scale: float) -> None:
        self.r *= scale
        self.height *= scale
        self.z_cut *= scale
        return self

    def scale_position(self, scale: float):
        self.center = [c * scale for c in self.center]
        return self
    
    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if (self.r <= 0.):
            raise ValueError("Angled Sample radius must be greater than 0")
        if (self.height <= 0.):
            raise ValueError("Angled Sample height must be greater than 0")
        if (self.z_cut <= 0.):
            raise ValueError("Angled Sample z_cut must be greater than 0")
        if abs(self.angle) >= math.pi / 2:
            raise ValueError("Angle must be below pi/2")

# Type alias for a single shape or list of shapes
Shape = Union[Rectangle, Disk, Annulus, AngledSample, Cube]
ShapeOrList = Union[Shape, list[Shape]]

@dataclass
class Object2D:
    shape: ShapeOrList
    label: Union[str, list[str]]
    surface_tag: Union[int, list[int], None] = None

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
        
    def scale_parameters(self, scale):
        self.size_min *= scale
        self.size_max *= scale
        return self

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
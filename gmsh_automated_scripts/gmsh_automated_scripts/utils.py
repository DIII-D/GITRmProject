#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:06:17 2024

@author: cappellil
"""
from .helpers import z_on_tilted_surface
from .utils2 import Disk, Rectangle, Annulus, MeshConfig
import gmsh
import math
import numpy as np
import os
from dataclasses import dataclass, field
from typing import Union, Optional

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
        self.center = [c * scale for c in self.center]
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


#%%
""" miscellaneous functions """

def rectangle_loop(x, y, z, width, height):
    # Create points for the rectangle corners
    p1 = gmsh.model.occ.addPoint(x, y, z)          # Bottom-left corner
    p2 = gmsh.model.occ.addPoint(x + width, y, z)  # Bottom-right corner
    p3 = gmsh.model.occ.addPoint(x + width, y + height, z)  # Top-right corner
    p4 = gmsh.model.occ.addPoint(x, y + height, z)  # Top-left corner
    
    # Create lines for the rectangle edges
    l1 = gmsh.model.occ.addLine(p1, p2)  # Bottom edge
    l2 = gmsh.model.occ.addLine(p2, p3)  # Right edge
    l3 = gmsh.model.occ.addLine(p3, p4)  # Top edge
    l4 = gmsh.model.occ.addLine(p4, p1)  # Left edge# Define a rectangle in the XY plane
    
    loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
    
    return p1, p2, p3, p4, l1, l2, l3, l4, loop


def annulus_loop(x, y, z, r_inner, r_outer, phi_start=0., angle=360.):
    """
    Flat annulus / annular sector on the plane z. Returns (dimes_outer_loop, inner_loops).
    Built entirely from SHARED vertices so the wire has no duplicate/orphan curves
    that would be left behind (stranded at z=0) when the surfaces are rotated.
    """
    if not (0 < angle <= 360):
        raise ValueError(f"angle must be in (0, 360], got {angle}")
    if r_inner >= r_outer:
        raise ValueError(f"r_inner ({r_inner}) must be less than r_outer ({r_outer})")

    # --- full ring: two concentric circles, central hole punched in the ring surface ---
    if abs(angle - 360.) < 1e-9:
        dimes_outer_loop = gmsh.model.occ.addCurveLoop([gmsh.model.occ.addCircle(x, y, z, r_outer)])
        inner_loop = gmsh.model.occ.addCurveLoop([gmsh.model.occ.addCircle(x, y, z, r_inner)])
        return dimes_outer_loop, [inner_loop]

    # --- sector: arcs through shared points via addCircleArc (center=True) ---
    pc   = gmsh.model.occ.addPoint(x, y, z)                 # shared circle center
    phi1 = math.radians(phi_start)
    phi2 = phi1 + math.radians(angle)
    nseg = 1 if angle < 180. else 2                         # one OCC arc must span < 180 deg
    phis = [phi1 + (phi2 - phi1) * k / nseg for k in range(nseg + 1)]

    pt = lambda r, p: gmsh.model.occ.addPoint(x + r * math.cos(p), y + r * math.sin(p), z)
    o_pts = [pt(r_outer, p) for p in phis]
    i_pts = [pt(r_inner, p) for p in phis]

    outer_arcs = [gmsh.model.occ.addCircleArc(o_pts[k], pc, o_pts[k + 1]) for k in range(nseg)]
    inner_arcs = [gmsh.model.occ.addCircleArc(i_pts[k], pc, i_pts[k + 1]) for k in range(nseg)]
    line_end   = gmsh.model.occ.addLine(o_pts[-1], i_pts[-1])
    line_start = gmsh.model.occ.addLine(i_pts[0],  o_pts[0])

    sector_loop = gmsh.model.occ.addCurveLoop(
        outer_arcs + [line_end] + [-a for a in reversed(inner_arcs)] + [line_start]
    )
    gmsh.model.occ.remove([(0, pc)])     # center point is only a construction aid; drop it
    return sector_loop, []


  
def roto_Ztranslation(curve, coo: list[float], rot_axis: list[float], angle: float, dimtag=1, translate = True, height = 0.):
    """object is rotated around an axis of revolution 
       AND translated along the z axis to ensure zmin = 0..
       Where z_center is the central z-coo before rotation."""
    
    x, y, z = coo
    ax, ay, az = rot_axis

    # 1. rotate around the origin
    gmsh.model.occ.rotate([(dimtag, curve)], x, y, z, ax, ay, az, angle)

    # 2. synchronize to query the bounding box after rotation
    gmsh.model.occ.synchronize()

    # 3. find zmin of the rotated curve
    x_min, y_min, z_min, x_max, y_max, z_max = gmsh.model.getBoundingBox(dimtag, curve)

    # 4. translate along Z so that zmin = z
    if translate:
        dz = z - z_min
        gmsh.model.occ.translate([(dimtag, curve)], 0, 0, dz + height)
        gmsh.model.occ.synchronize()

def truncated_cylinder_Ellipse(x, y, z, r, ax, ay, angle):
    
    if angle != 0.:
        r1 = r / math.cos(angle) if abs(ay) > abs(ax) else r  # stretch x when rotating about y
        r2 = r / math.cos(angle) if abs(ax) > abs(ay) else r  # stretch y when rotating about x

        # if r1 != r2 disk perimeter is an Ellipse
        curve = gmsh.model.occ.addEllipse(x , y , z , r1, r2) # r1 along x, r2 along y by default from .addEllipse method
    else:
        curve = gmsh.model.occ.addCircle(x, y, z, r)

    return curve

def truncated_cylinder_loop(x, y, z, r,
              ax, ay, az,
              angle, height = 0.):
    
    """
    The top surface of a tilted DiMES head or a tilted button is approx. a slice of a cylinder (top surface of a truncated cylinder).

    r1 along x, r2 along y by default from .addEllipse method

    if rotation about y-axis -> r1 increases
    if rotation about x-axis -> r2 increases
    """

    curve = truncated_cylinder_Ellipse(x, y, z, r, ax, ay, angle)

    roto_Ztranslation(curve, [x, y, z], [ax, ay, az], angle, translate=True, height=height)

    # add loop after rotation. Remember: you can't rotate loops

    loop = gmsh.model.occ.addCurveLoop([curve]) 

    return curve, loop # -> tuple(curve_tag:int, loop_tag:int) 

def get_surface_center(surface_tag: int) -> tuple[float, float, float]:
    """Returns the (x, y, z) center of mass of a planar surface."""
    gmsh.model.occ.synchronize()
    x, y, z = gmsh.model.occ.getCenterOfMass(2, surface_tag)
    return x, y, z


def get_surface_rotation(surface_tag: int) -> tuple[float, float, float]:
    """
    Returns the rotation (rx, ry, rz) in radians of a planar surface
    by computing its normal vector and deriving Euler angles from it.
    """
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)

    # Get nodes on the surface and their parametric coordinates
    _, _, param = gmsh.model.mesh.getNodes(2, surface_tag, includeBoundary=True)

    # Get the normal at the first node (all normals are equal on a planar surface)
    normal = gmsh.model.getNormal(surface_tag, param[:2])
    nx, ny, nz = normal[0], normal[1], normal[2]
    n = np.array([nx, ny, nz])
    n = n / np.linalg.norm(n)  # ensure unit vector

    # Rotation around X axis (pitch): angle between n and the XZ plane
    rx = np.arctan2(ny, nz)

    # Rotation around Y axis (roll): angle between n and the YZ plane
    ry = np.arctan2(nx, nz)

    # Rotation around Z axis (yaw): angle of the projection onto the XY plane
    rz = np.arctan2(ny, nx)

    return rx, ry, rz


# 3. add samples loops to top surface
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
        from .export import save_component_meshes
        base_path = os.path.splitext(filename)[0]
        save_component_meshes(component_surfaces, base_path, save_npz=_save_npz,
                              save_ply=save_ply, ply_normals=ply_normals, ply_ascii=ply_ascii)

    # Finalize GMSH
    gmsh.finalize()

# 0. Create Tile around DiMES

gmsh.initialize()

L_tile = 10.

Tile = Rectangle(L_tile,L_tile,[-L_tile/2,-L_tile/2, 0.])

*_, tile_loop = rectangle_loop(*Tile.ll, Tile.width, Tile.height) 

# 1. create DiMES top surface as truncated cylinder

# 1.1 create Ellipse

DimesTop = Disk(4., [0.,0.,0.])
rot_axis = [0.,-1.,0.]
angle = np.deg2rad(5)

curve = truncated_cylinder_Ellipse(*DimesTop.center, DimesTop.r, rot_axis[0], rot_axis[1], angle)
dimes_outer_loop = gmsh.model.occ.addCurveLoop([curve])

# 1.3 create DiMES hole and tile surface

tile_surface = gmsh.model.occ.addPlaneSurface([tile_loop, dimes_outer_loop])

# 1.4 Create hole loops and surfaces

shapes = [Rectangle(1,1,[-0.5,-0.5, 0.]), Disk(0.5, [-2, 0., 0.]), Annulus(3.3, 3.8, (90, 270), [0., 0., 0.])]

inner_loops, inner_surfaces = [], []
for shape in shapes:
    loop, surface = add_sample(shape)
    inner_loops.append(loop)
    inner_surfaces.append(surface)

# 1.5 Create surface with holes
dimes_top_surface = gmsh.model.occ.addPlaneSurface([dimes_outer_loop] + inner_loops)

# 1.6 Rotate surfaces
if angle != 0.:
    ax, ay, az = rot_axis
    cx, cy, cz = DimesTop.center

    all_surfaces = [(2, dimes_top_surface)] + [(2, s) for s in inner_surfaces]

    gmsh.model.occ.rotate(all_surfaces, cx, cy, cz, ax, ay, az, angle)   # surfaces only
    gmsh.model.occ.synchronize()
    _, _, z_min, _, _, _ = gmsh.model.getBoundingBox(2, dimes_top_surface)
    gmsh.model.occ.translate(all_surfaces, 0, 0, cz - z_min)
    gmsh.model.occ.synchronize()

# 1.7 Add DiMES side surface
if angle != 0.:
    dimes_base_loop = gmsh.model.occ.addCurveLoop(
        [gmsh.model.occ.addCircle(*DimesTop.center, DimesTop.r)]
    )
    result = gmsh.model.occ.addThruSections(
        [dimes_base_loop, dimes_outer_loop], makeSolid=False
    )
    gmsh.model.occ.synchronize()
    # result is a list of (dim, tag) tuples — extract the surface tag
    dimes_side_surface = result[0][1]
    print(f"Side surface tag: {dimes_side_surface}")  # verify it's an int


# mesh = MeshConfig()
# Add Mesh Options

component_surfaces = {
    "tile":       [tile_surface],
    "dimes_top":  [dimes_top_surface],
    "dimes_side": [dimes_side_surface],   # an int, not a tuple
    **{f"sample_{i}": [s] for i, s in enumerate(inner_surfaces)}
}

make_dimes_mesh(
    mesh=None,
    component_surfaces=component_surfaces,
    GUI_msh=True,
    GUI_geo=False
)

gmsh.model.mesh.generate(2)

# Run gmsh
gmsh.fltk.run()

gmsh.finalize()






def create_loops(input_dict, z_dimes, volumes_surfaces, dot_loops, ax, ay, az, theta_dimes):
    
    # check if user is not tilting both DiMES head and dots
    
    if theta_dimes != 0:
    
        for elem_def in input_dict.values():
    # Check if theta_dot exists and is not zero
            if "theta_dot" in elem_def and elem_def["theta_dot"] != 0:
                print(f"{elem_def}: theta_dot is not zero (theta_dot = {elem_def['theta_dot']})")
                raise Exception("You are setting a tilted dot (theta_dot != 0) " \
                      "on a tilted DiMES head (theta_dimes !=0), this could be problematic")


    component_surfaces = {}

    for name, elem_def in input_dict.items():
    
        if elem_def["shape"] == "circle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            r = elem_def['radius']
            theta_dot = theta_dimes + math.pi / 180 * elem_def.get('theta_dot', 0)
            
            # z_dimes_tilt indicates the z position of dot over a tilted DiMES
            
            z = z_on_tilted_surface(z_dimes, x, y, ax, ay, theta_dimes)

            # if dot coating surface is not tilted, dot simulated as a Disk
            # coplanar with the DiMES head
            
            # otherwise dot simulated as a surface delimitated by an allipse at the top
            # and a circle at the base
            
            # create base circle coplanar with DiMES head
            r1 = r / math.cos(theta_dimes) # greater radius
            r2 = r # smaller radius
            
            dot_base_curve = gmsh.model.occ.addEllipse(x , y  , z , r1, r2)
            
            # if you want to rotate the dot around the x-axis, since the ellipse must
            # be turned before around the z-axis of 90 deg so that the major radius is along the y-axis
            # also the circle must be rotated otherwise the addThruSection function won't work properly
            # when you want to create the side surface
            
            if ax != 0:
                gmsh.model.occ.rotate([(1, dot_base_curve)], x,
                                      y, z, 0, 0, 1, math.pi / 2)
                
            
            gmsh.model.occ.rotate([(1 , dot_base_curve)], x, y, z, ax, ay, az, theta_dimes)
            
            dot_base_loop = gmsh.model.occ.addCurveLoop([dot_base_curve])
            
            # append base loop to list of holes to create DiMES head surface
            #dot_loops.append(dot_base_loop)
            
            
            if theta_dot > theta_dimes:

                if "deposits" in elem_def and elem_def["deposits"]:
                    raise Exception("Deposits on tilted dots are not supported")

                #---------------------

                # PLEASE NOTE: this code only works for rotations
                # - about the y-axis (ax=0, ay=±1, az=0)
                # - about the x-axis (ax=±1, ay=0, az=0)
                #
                # ** all other combinations might not work for rectangular dots
                #
                #
                # ** if theta_dimes !=0  dots should be centered symmetrical along axis perp.
                # to rotation axis. For instace (if rotation about y, dots should be 
                # located in a position where their central position (x_center==0) .
                # Otherwise you might encounter issues related to the shift along the z-axis.
                # For the moment it is better to avoid setting theta_dimes!=0 when 
                # also theta_dot!=0
                
                #---------------------
                
                # create top ellipse:
                
                # 1. shift ellipse along z to avoid intersection with DiMES head surface
                # after rotation
                
                delta_z = r * abs(math.tan(theta_dot - theta_dimes)) + 0.0001 # + 0.0001 to avoid intersecting facets
                
                z += delta_z
                
                # 2. make top elliptical curve
                dot_top_curve = gmsh.model.occ.addEllipse(x  , y , z, r / math.cos(theta_dot), r)
                
                
                # 3. rotate top elliptical curve (remember: you can't rotate loops) and make surface
                
                if ax != 0:
                    gmsh.model.occ.rotate([(1, dot_top_curve)], x,
                                          y, z, 0, 0, 1, math.pi / 2)
                
                gmsh.model.occ.rotate([(1 , dot_top_curve)], x, y, z, ax, ay, az, theta_dot)
                
                dot_top_loop = gmsh.model.occ.addCurveLoop([dot_top_curve])
                dot_top_surface = gmsh.model.occ.addPlaneSurface([dot_top_loop])
                
                # 4. create side surface
                dot_side_surface = gmsh.model.occ.addThruSections([dot_base_loop, dot_top_loop], makeSolid = False)
                
                # 5. append side and top surfaces to list of surface delimiting the plasma volume
                volumes_surfaces.append(dot_side_surface[0][1])
                volumes_surfaces.append(dot_top_surface)
                component_surfaces[name] = [dot_side_surface[0][1], dot_top_surface]

            else:

                # collect child deposit loops to perforate this dot surface
                child_loops = []
                child_surfaces = {}
                if "deposits" in elem_def and elem_def["deposits"]:
                    child_surfaces = create_loops(elem_def["deposits"], z_dimes, volumes_surfaces,
                                                  child_loops, ax, ay, az, theta_dimes)

                dot_surface = gmsh.model.occ.addPlaneSurface([dot_base_loop] + child_loops)

                # append dot_base_loop to list of holes composing parent (DiMES or dot) surface
                dot_loops.append(dot_base_loop)

                # append disk surface to list of surfaces delimiting the plasma volume
                volumes_surfaces.append(dot_surface)
                component_surfaces[name] = [dot_surface]
                component_surfaces.update(child_surfaces)
            
                    
        elif elem_def["shape"] == "rectangle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            width = elem_def['width']
            height = elem_def['height']
            theta_dot = theta_dimes + math.pi / 180 * elem_def.get('theta_dot', 0)
            
            # if dot coating surface is not tilted, dot simulated as a Rectangular surface
            # coplanar with the DiMES head
            
            # otherwise dot simulated as a wedge delimitated by a tilted plane at the top
            # and a plane coplanar to DiMES head at the base
            
            # z_dimes_tilt indicates the z position of dot over a tilted DiMES

            z = z_on_tilted_surface(z, x, y, ax, ay, theta_dimes)
            
            # translate base rectangle to be coplanar with DiMES head and to avoid intersections
            delta_z = 0.001 # to avoid overlapping between curves
            delta_z_dimes = 0

            # 1. Create base lines and curve loop shifted along z
            base_l1, base_l2, base_l3, base_l4, dot_base_loop = rectangle_loop(x, y, z, width , height)[-5:]

            # 2. rotate rectangle using base_l1 edge as pivotal point around y-axis
            #    and create loop

            gmsh.model.occ.rotate([(1, base_l1), (1, base_l2), (1, base_l3), (1, base_l4)], x, y, z, ax, ay, az, theta_dimes)
            dot_base_loop = gmsh.model.occ.addCurveLoop([base_l1, base_l2, base_l3, base_l4])
            
            # append dot_base_loop to list of holes composing DiMES head surface
            dot_loops.append(dot_base_loop)

            if theta_dot > theta_dimes:

                if "deposits" in elem_def and elem_def["deposits"]:
                    raise Exception("Deposits on tilted dots are not supported")

                #---------------------

                # PLEASE NOTE: this code only works for rotations
                # - about the y-axis (ax=0, ay=±1, az=0)
                # - about the x-axis (ax=±1, ay=0, az=0)
                #
                # ** all other combinations might not work
                #
                #
                # ** if theta_dimes !=0  dots should be centered symmetrical along axis perp.
                # to rotation axis. For instace (if rotation about y, dots should be
                # located in a position where their central position (x_center==0) .
                # Otherwise you might encounter issues related to the shift along the z-axis.
                # For the moment it is better to avoid setting theta_dimes!=0 when
                # also theta_dot!=0

                #---------------------

                #---------------------

                # PLEASE NOTE 2: even if theta_dimes==0 for rectangular shapes
                # this code only works for rotations about the y-axis (ax=0, ay=±1, az=0)
                # and the x-axis (ax=±1, ay=0, az=0)
                
                #---------------------
                
                # Create top plane:
                    
                # actualization of z axis at dot position
                
                x_c = x + width / 2
                
                y_c = y + height / 2

                z = z_on_tilted_surface(z_dimes, x_c, y_c, ax, ay, theta_dimes)
                
                # 1. Get the top plane individual lines (curves) lengths before rotation
                w = width * math.cos(theta_dimes) / math.cos(theta_dot) # width of rotated rectangle
                side_edge_height = w * math.sin(theta_dot - theta_dimes) # maximum height of wedge along z
            
                # 1a. add side_edge_height to z shift because rotation always happens about l1
                if ay == 1:
                    delta_z_dimes += side_edge_height 
                    
                if ax == -1:
                    h = height * math.cos(theta_dimes) / math.cos(theta_dot) # width of rotated rectangle
                    side_edge_height = h * math.sin(theta_dot - theta_dimes)
                    delta_z_dimes += side_edge_height 
                    
                
                # 2. create top plane
                l1, l2, l3, l4, dot_loop = rectangle_loop(x, y, z + delta_z_dimes + delta_z, w, height)[-5:]
        
                # 3. Rotate top plane by rotating its individual lines l1, l2, l3, l4
                gmsh.model.occ.rotate([(1, l1), (1, l2), (1, l3), (1, l4)], x, y, z + delta_z_dimes + delta_z, ax, ay, az, theta_dot)
        
                # Synchronize to apply the rotation
                gmsh.model.occ.synchronize()
        
                # Get the rotated curve IDs (Gmsh might assign new IDs after rotation, so we synchronize first)
                new_l1, new_l2, new_l3, new_l4 = [l1, l2, l3, l4]  # In Gmsh, the curve IDs should stay the same, but we reassign them for clarity
        
                # Rebuild the rotated curve loop from the rotated curves
                dot_rotated_loop = gmsh.model.occ.addCurveLoop([new_l1, new_l2, new_l3, new_l4])

                # 4. create top surface
                dot_surface = gmsh.model.occ.addPlaneSurface([dot_rotated_loop])
        
                # 5. Create the through section (side surface) between the original and rotated loops
                dot_side_surface = gmsh.model.occ.addThruSections([dot_base_loop, dot_rotated_loop], makeSolid=False)
        
                # 6. append side and top surfaces to list of surface delimiting the plasma volume
                for dim, tag in dot_side_surface:
                    volumes_surfaces.append(tag)

                volumes_surfaces.append(dot_surface)
                component_surfaces[name] = [tag for _, tag in dot_side_surface] + [dot_surface]

                # Synchronize the model to update geometry
                gmsh.model.occ.synchronize()
            else:
                # collect child deposit loops to perforate this dot surface
                child_loops = []
                child_surfaces = {}
                if "deposits" in elem_def and elem_def["deposits"]:
                    child_surfaces = create_loops(elem_def["deposits"], z_dimes, volumes_surfaces,
                                                  child_loops, ax, ay, az, theta_dimes)

                dot_base_surface = gmsh.model.occ.addPlaneSurface([dot_base_loop] + child_loops)

                # append rectangle surface to list of surfaces delimiting the plasma volume
                volumes_surfaces.append(dot_base_surface)
                component_surfaces[name] = [dot_base_surface]
                component_surfaces.update(child_surfaces)

        elif elem_def["shape"] == "annulus":

            x         = elem_def['x']
            y         = elem_def['y']
            r_inner   = elem_def['r_inner']
            r_outer   = elem_def['r_outer']
            phi_start = elem_def.get('phi_start', 0.)
            if 'phi_end' in elem_def:
                angle = elem_def['phi_end'] - phi_start
            else:
                angle = elem_def.get('angle', 360.)
            theta_dot = theta_dimes + math.pi / 180 * elem_def.get('theta_dot', 0)

            if not (0 < angle <= 360):
                raise Exception(f"angle must be in (0, 360], got {angle}")
            if r_inner >= r_outer:
                raise Exception(f"r_inner ({r_inner}) must be less than r_outer ({r_outer})")
            if theta_dot > theta_dimes:
                raise Exception("Tilted annuli (theta_dot != 0) are not supported")

            z        = z_on_tilted_surface(z_dimes, x, y, ax, ay, theta_dimes)
            r1_outer = r_outer / math.cos(theta_dimes)
            r1_inner = r_inner / math.cos(theta_dimes)

            child_loops    = []
            child_surfaces = {}
            if "deposits" in elem_def and elem_def["deposits"]:
                child_surfaces = create_loops(elem_def["deposits"], z_dimes, volumes_surfaces,
                                              child_loops, ax, ay, az, theta_dimes)

            if abs(angle - 360.) < 1e-9:
                outer_curve = gmsh.model.occ.addEllipse(x, y, z, r1_outer, r_outer)
                inner_curve = gmsh.model.occ.addEllipse(x, y, z, r1_inner, r_inner)
                if ax != 0:
                    gmsh.model.occ.rotate([(1, outer_curve), (1, inner_curve)],
                                          x, y, z, 0, 0, 1, math.pi / 2)
                gmsh.model.occ.rotate([(1, outer_curve), (1, inner_curve)],
                                      x, y, z, ax, ay, az, theta_dimes)
                dimes_outer_loop = gmsh.model.occ.addCurveLoop([outer_curve])
                inner_loop = gmsh.model.occ.addCurveLoop([inner_curve])
                annulus_surface = gmsh.model.occ.addPlaneSurface(
                    [dimes_outer_loop, inner_loop] + child_loops)
                dot_loops.append(dimes_outer_loop)

            else:
                phi1 = math.radians(phi_start)
                phi2 = phi1 + math.radians(angle)
                c1, s1 = math.cos(phi1), math.sin(phi1)
                c2, s2 = math.cos(phi2), math.sin(phi2)

                p_o1 = gmsh.model.occ.addPoint(x + r1_outer * c1, y + r_outer * s1, z)
                p_o2 = gmsh.model.occ.addPoint(x + r1_outer * c2, y + r_outer * s2, z)
                p_i1 = gmsh.model.occ.addPoint(x + r1_inner * c1, y + r_inner * s1, z)
                p_i2 = gmsh.model.occ.addPoint(x + r1_inner * c2, y + r_inner * s2, z)

                outer_arc = gmsh.model.occ.addEllipse(x, y, z, r1_outer, r_outer,
                                                      angle1=phi1, angle2=phi2)
                inner_arc = gmsh.model.occ.addEllipse(x, y, z, r1_inner, r_inner,
                                                      angle1=phi1, angle2=phi2)
                line_end   = gmsh.model.occ.addLine(p_o2, p_i2)
                line_start = gmsh.model.occ.addLine(p_i1, p_o1)

                all_ents = [(1, outer_arc), (1, inner_arc), (1, line_end), (1, line_start),
                            (0, p_o1), (0, p_o2), (0, p_i1), (0, p_i2)]
                if ax != 0:
                    gmsh.model.occ.rotate(all_ents, x, y, z, 0, 0, 1, math.pi / 2)
                gmsh.model.occ.rotate(all_ents, x, y, z, ax, ay, az, theta_dimes)
                gmsh.model.occ.synchronize()

                sector_loop = gmsh.model.occ.addCurveLoop(
                    [outer_arc, line_end, -inner_arc, line_start])
                annulus_surface = gmsh.model.occ.addPlaneSurface(
                    [sector_loop] + child_loops)
                dot_loops.append(sector_loop)

            volumes_surfaces.append(annulus_surface)
            component_surfaces[name] = [annulus_surface]
            component_surfaces.update(child_surfaces)

    return component_surfaces


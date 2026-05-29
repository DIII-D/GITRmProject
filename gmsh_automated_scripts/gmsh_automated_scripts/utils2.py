#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:07:27 2024

@author: cappellil
"""
import os
import gmsh
import math
from .utils import rectangle_def, create_loops
from dataclasses import dataclass
from typing import Union

@dataclass
class Cube:
    L_R: float      # radial length
    L_phi: float    # toroidal length
    L_Z: float      # z-axis length
    center: list[float]  # [R, phi, Z]

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def scale_size(self, scale: float) -> None:
        self.L_R *= scale
        self.L_phi *= scale
        self.L_Z *= scale
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

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if self.r <= 0.:
            raise ValueError("Disk radius must be greater than 0")

@dataclass
class Annulus:
    r_minor: float
    r_major: float
    angular_sector: tuple
    center: list[float] # [R, phi, Z]

    def scale_size(self, scale: float) -> None:
        self.r_minor *= scale
        self.r_major *= scale
        return self

    def __post_init__(self): # dataclass calls __post_init__ automatically right after construction
        self.check_input()

    def check_input(self) -> None:
        if len(self.center) != 3:
            raise ValueError("len(self.center) must be equal to 3")
        if (self.r_minor <= 0.):
            raise ValueError("Annulus r_minor must be greater than 0")
        if (self.r_major <= 0.):
            raise ValueError("Annulus r_major must be greater than 0")
        if (self.angular_sector[0] >= self.angular_sector[1]):
            raise ValueError("Annulus angular sector final angle must be smaller than the first angle.")

@dataclass
class Samples:
    items: list[Union[Cube, Disk, Annulus]]

def make_dimes_geom(input_dict, L_R=8, L_phi=8, L_Z=8, box_center = [0.,0.,0.],
                    r_dimes = 2.5, dimes_center = [0., 0., 0.],
                    ax=0, ay=-1, az=0, theta_dimes=0, no_dots=False, scale = 1.):
    
    boxCube = Cube(L_R, L_phi, L_Z, box_center).scale_size(scale)
    DiMESTop = Disk(r_dimes, dimes_center).scale_size(scale)

    print(f"type(boxCube) = {type(boxCube)}")

    # convert deg to rad
    theta_dimes = math.pi / 180 * theta_dimes
    
    # avoid not allowed rotation angles to be used
    
    ax_rot = (ax, ay, az)
        
    ax_rot_allowed = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
        
    if ax_rot in ax_rot_allowed:
        pass
    else:
        raise Exception(f"You entered a combination of {ax, ay, az} which " \
                        "is not allowed. Please inspect function to see " \
                            "directions you can use")

    # %%

    """ Plasma volume geometry """
    # Cartesian coordinates of bottom plasma volume surface lower left corner
    
    # boxCube: coordinates of lower left vertex
    R_ll = -boxCube.L_R / 2 + boxCube.center[0]  # radial
    phi_ll = -boxCube.L_phi / 2 + boxCube.center[1]  # toroidal
    z_ll = boxCube.center[2]  # vertical

    # Create a recangular curve loop for box base
    p1, p2, p3, p4, l1, l2, l3, l4, base_rectangle_loop = \
        rectangle_def(R_ll, phi_ll, z_ll, boxCube.L_R, boxCube.L_phi)

    # Create a rectangular curve loop for box top
    p5, p6, p7, p8, l5, l6, l7, l8, top_rectangle_loop = \
        rectangle_def(R_ll, phi_ll, z_ll + boxCube.L_Z, boxCube.L_R, boxCube.L_phi)

    # Create vertical lines connecting bottom and top
    l9 = gmsh.model.occ.addLine(p1, p5)
    l10 = gmsh.model.occ.addLine(p2, p6)
    l11 = gmsh.model.occ.addLine(p3, p7)
    l12 = gmsh.model.occ.addLine(p4, p8)

    # list of surfaces IDs enclosing the plasma volume
    volumes_surfaces = []

    # Create surfaces for the sides (lateral surfaces)
    # minus when instead of going for extremity A to B you go from B to A close a loop
    box_side1 = gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([l1, l10, -l5, -l9])])
    box_side2 = gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([l2, l11, -l6, -l10])])
    box_side3 = gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([l3, l12, -l7, -l11])])
    box_side4 = gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([l4, l9, -l8, -l12])])

    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(box_side1)
    volumes_surfaces.append(box_side2)
    volumes_surfaces.append(box_side3)
    volumes_surfaces.append(box_side4)

    # Create the top surface
    box_top_surface = gmsh.model.occ.addPlaneSurface(
        [gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])])

    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(box_top_surface)

    # Synchronize the GMSH model
    gmsh.model.occ.synchronize()

    # %%
    """ DiMES geometry """
        
    DiMES_base_circle = gmsh.model.occ.addCircle(*DiMESTop.center, DiMESTop.r)

    # if you want to rotate the DiMES head around the x-axis, the ellipse must
    # be turned before around the z-axis of 90 deg so that the major radius is along the y-axis.
    # The circle must be rotated as well otherwise the addThruSection function won't work properly
    # when you want to create the side surface of DiMES
    
    if ax != 0:
        gmsh.model.occ.rotate([(1, DiMES_base_circle)], *DiMESTop.center, 0, 0, 1, math.pi / 2)

    DiMES_base_circle_loop = gmsh.model.occ.addCurveLoop([DiMES_base_circle])

    box_base_surface = gmsh.model.occ.addPlaneSurface(
        [base_rectangle_loop, DiMES_base_circle_loop])

    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(box_base_surface)

    if theta_dimes != 0:
        
        # translate along z to avoid overlapping with base surface
        delta_z_dimes = DiMESTop.r * math.tan(theta_dimes)         
        DiMESTop.center[2] += delta_z_dimes + 0.001*scale # + 0.001*scale to avoid overlapping facets from different surfaces
        
        # once tilted, the disk perimeter turns from a circle into an ellipse
        r_major = r_dimes / math.cos(theta_dimes)  # major radius
        r_minor = r_dimes  # minor radius (r1 >= r2)
        
        top_ellipse = gmsh.model.occ.addEllipse(*DiMESTop.center, r_major, r_minor)
        
        # addEllipse only creates ellipses with major radius along x-axis
        
        # if you want to rotate the DiMES head around the x-axis, the ellipse must
        # be turned before around the z-axis of 90 deg so that the major radius is along the y-axis
        
        if ax != 0:
            gmsh.model.occ.rotate([(1, top_ellipse)], *DiMESTop.center, 0, 0, 1, math.pi / 2)
            

        # rotate about rotation axis about an angle equal to theta_dimes 
        gmsh.model.occ.rotate([(1, top_ellipse)], *DiMESTop.center, ax, ay, az, theta_dimes)
        
        top_ellipse_loop = gmsh.model.occ.addCurveLoop([top_ellipse])

        # Volumes and surfaces can be constructed from (closed) curve loops thanks to the
        # `addThruSections()' function

        DiMES_side_surface = gmsh.model.occ.addThruSections(
            [DiMES_base_circle_loop, top_ellipse_loop], makeSolid=False)

        # store surfaces enclosing volume in a variable
        for tup in DiMES_side_surface:
            if tup[0] == 2:
                DiMES_side_surface_id = tup[1]
                break

        volumes_surfaces.append(DiMES_side_surface_id)

        # identify top_circle loop
        DiMES_top_curve_loop = top_ellipse_loop

    else:
        DiMES_top_curve_loop = DiMES_base_circle_loop

    # %%

    """ geometry Dots (coatings)"""

    dot_loops = []
    
    # if theta_dimes != 0:
    #     pass
    # else:
        
    dot_component_surfaces = {}
    if no_dots:
        pass
    else:
        dot_component_surfaces = create_loops(input_dict, DiMESTop.center[2], volumes_surfaces,
                                              dot_loops, ax, ay, az, theta_dimes)

    # Generate DiMES top surface
    gmsh.model.occ.synchronize()
    DiMES_top_surface = gmsh.model.occ.addPlaneSurface([DiMES_top_curve_loop] + dot_loops)

    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(DiMES_top_surface)

    # %% generate the volume
    plasma_volume = gmsh.model.occ.addVolume(
        [gmsh.model.occ.addSurfaceLoop(volumes_surfaces)])

    # build component surface map (name → list of OCC surface tags)
    component_surfaces = {
        "box_side1": [box_side1],
        "box_side2": [box_side2],
        "box_side3": [box_side3],
        "box_side4": [box_side4],
        "box_top": [box_top_surface],
        "box_base": [box_base_surface],
        "DiMES_top": [DiMES_top_surface]
    }
    if theta_dimes != 0:
        component_surfaces["DiMES_side"] = [DiMES_side_surface_id]
    component_surfaces.update(dot_component_surfaces)

    return plasma_volume, component_surfaces


# %%
""" function """


def make_dimes_mesh(filename="test.msh", save_msh=False, GUI_geo=False, GUI_msh=True,
                    msh_dim=3, component_surfaces=None, save_output = False, save_npz=None, save_ply=False,
                    ply_normals=True, ply_ascii=True):
    # %% Generate the mesh and visualize the result

    # Remove duplicates (coherence)
    gmsh.model.occ.removeAllDuplicates()

    # Final synchronization of the CAD model
    gmsh.model.occ.synchronize()

    # Label each component as a named Physical Surface (2-D meshes only)
    if msh_dim == 2 and component_surfaces is not None:
        for comp_name, surf_tags in component_surfaces.items():
            pg = gmsh.model.addPhysicalGroup(2, surf_tags)
            gmsh.model.setPhysicalName(2, pg, comp_name)

    if GUI_geo:
        gmsh.fltk.run()
    
    # meshing options
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)

    # Prevent very small elements in small dots
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
    # Set maximum mesh characteristic length for the whole model
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)
    gmsh.model.mesh.generate(msh_dim)

    # run GUI
    if GUI_msh:
        gmsh.fltk.run()

    # save .msh file 
    if save_msh:
        gmsh.write(filename)

    # save_npz=None means "follow save_msh"; explicit True/False overrides
    _save_npz = save_msh if save_npz is None else save_npz

    # Save per-component files for 2-D meshes (independent of save_msh)
    if msh_dim == 2 and (_save_npz or save_ply) and save_output and component_surfaces is not None:
        from .export import save_component_meshes
        base_path = os.path.splitext(filename)[0]
        save_component_meshes(component_surfaces, base_path, save_npz=_save_npz, save_ply=save_ply,
                              ply_normals=ply_normals, ply_ascii=ply_ascii)

    # Finalize GMSH
    gmsh.finalize()


# %%
""" function """

def generate_dimes_mesh(Box: Cube, DiMESTop: Disk, samples : Samples):
    return


def generate_dimes_mesh_legacy(input_dict, **kwargs):
    """
This function automatically generates the DiMES mesh for a given set of dot geometries.

### PARAMETERS: (KEYWORDS DESCRIPTION IS AFTER kwargs)

- input_dict: A dictionary containing the geometries of the dots to be simulated, where each dot is defined by a unique label followed by its properties.

- kwargs: Additional keyword arguments. Supported keywords include:

  - x_center_dimes: The X-coordinate of the center of the DiMES base (float, default: 0).
  - y_center_dimes: The Y-coordinate of the center of the DiMES base (float, default: 0).
  - z_top_dimes: The Z-coordinate of the center of the DiMES base (float, default: 0).
  - r_dimes: The radius of the DiMES base (float, default: 2.5).
  - x_center: The X-coordinate of the plasma volume base surface center (float, default: 0).
  - y_center: The Y-coordinate of the plasma volume base surface center (float, default: 0).
  - z_center: The Z-coordinate of the plasma volume base surface center (float, default: 0).
  - L_R: The extent of plasma volume in the radial direction (half-length, float, default: 8).
  - L_phi: The extent of plasma volume in the toroidal direction (half-length, float, default: 8).
  - L_Z: The extentof plasma volume  in the vertical direction (float, default: 8).
  - ax: The X-component of the rotation axis (float, default: 0).
  - ay: The Y-component of the rotation axis (float, default: -1).
  - az: The Z-component of the rotation axis (float, default: 0).
  - msh_dim: The dimension of the mesh, which can be 1D, 2D, or 3D (int, default: 3).
  - filename: The name of the mesh file to be generated (string, default: "test.msh").
  - save_msh: A flag indicating whether to save the mesh (bool, default: False).
  - GUI_geo: A flag to use a GUI for visualizing the geometry before meshing (bool, default: False).
  - GUI_msh: A flag to use a GUI for visualizing both geometry and mesh (bool, default: True).
  - no_dots: A flag to remove all coatings (bool, default: False).

### GEOMETRY KEYWORDS:

The DiMES geometry can take one of the following shapes:

- disk (DiMES head flushed with lower divertor)
- cylinder (DiMES head protruding in the plasma but top face still parallel to divertor)
- truncated cylinder (DiMES head tilted)

On top of the DiMES head different coatings can be added as discussed later.

The DiMES geometry is defined by its center and radius:

- `x_center_dimes`: The X-coordinate of the center of the DiMES base (float).
- `y_center_dimes`: The Y-coordinate of the center of the DiMES base (float).
- `z_top_dimes`: The Z-coordinate of the center of the DiMES base (float).
- `r_dimes`: The radius of the DiMES base (float).
- `theta_dimes`: Tilting angle
-  `ax`, `ay` ,`az`: components of rotation' direction (gmsh rotates objects using Rodrigues formula)'

Right now only the following 4 directions are allowed:
    
    ax = ± 1, ay = 0, az = 0
    ax = 0, ay = ± 1, az = 0

#### Plasma Volume:

The DiMES geometry is enclosed within a volume known as `plasma_volume`, representing the plasma background. This volume is a symmetrical box extending in three directions:

- `x_center`: The X-coordinate of the plasma volume center (float).
- `y_center`: The Y-coordinate of the plasma volume center (float).
- `z_center`: The Z-coordinate of the plasma volume center (float).

The plasma volume extends symmetrically in these directions:

- `L_R`: Twice the radial extent of the volume (float).
- `L_phi`: Twice the toroidal extent of the volume (float).
- `L_Z`: The vertical extent of the volume (float).

#### Dots (Material Coatings):

On the DiMES top surface, dots represent the material coatings. Two geometries are currently supported:

- circle
- rectangle
- annulus

Each dot's position and dimensions are specified within an `input_dict` dictionary. The user must manually create this dictionary, where each dot is defined by a unique label followed by values that describe its geometry, shape, and position.

Each dot can be tilted _around the same rotation axis as DiMES_, identified by its
components (`ax`, `ay`, `az`) and the angle `theta_dot` (deg). 

`theta_dot` can be set in the `input_dict`, as shown later.

**NOTE: tilting both DiMES head and dots can result in wrong geometry and mesh generation! Try to avoid it.**

Right now only the following 4 directions are allowed:
    
    ax = ± 1, ay = 0, az = 0
    ax = 0, ay = ± 1, az = 0

#### Example Input Dictionary:

This example illustrates the structure for two dots: one circle and one rectangle.

```

Names = ["Dot_1", "Dot_2"]

input_dict = {
    Names[0]: {
        "shape": "circle",
        "x": 0,
        "y": 0.75,
        "radius": 0.05  # Radius for the circle
        "theta_dot": 10
    },
    Names[1]: {
        "shape": "rectangle",
        "x": -0.5,
        "y": -0.25,
        "width": 1,    # Width of the rectangle
        "height": 0.5  # Height of the rectangle
        "theta_dot": 5
    }
    # Additional dots can be added as needed
}

```

Each dot has the following properties:

- `x`: The X-coordinate of the dot’s position on the DiMES top surface (float).
- `y`: The Y-coordinate of the dot’s position on the DiMES top surface (float).

For "Dot_1" (shape: "circle"), the `x` and `y` coordinates represent the center of the circle, and the `radius` defines its size.

For "Dot_2" (shape: "rectangle"), the `x` and `y` coordinates denote the lower-left corner of the rectangle, while the following values define its size:

- `width`: The width of the rectangle (float).
- `height`: The height of the rectangle (float).

For an annulus (shape: "annulus"), all parameters are:

- `x`: X-coordinate of the annulus centre (float).
- `y`: Y-coordinate of the annulus centre (float).
- `r_inner`: Inner radius (float, must be < `r_outer`).
- `r_outer`: Outer radius (float).
- `phi_start`: Start angle in degrees measured from the +X axis (float, default 0).
- `phi_end`: End angle in degrees measured from the +X axis (float, optional). When provided, takes priority over `angle`.
- `angle`: Angular extent of the sector in degrees (float, default 360 = full annulus). Ignored when `phi_end` is set.

#### Deposits (Dots on Dots):

Any flat (untilted) dot can itself act as a substrate by adding a `deposits` key whose value is a nested `input_dict` with the same structure. Deposits are coplanar with their parent dot; the parent surface is automatically perforated where each deposit sits, exactly as the DiMES head is perforated for top-level dots.

Deposits can be nested to arbitrary depth. Tilted dots (`theta_dot != 0`) cannot carry deposits; attempting this raises an exception.

```

Names = ["Dot_1", "Dot_2", "Dot_1a"]

input_dict = {
    Names[0]: {
        "shape": "circle",
        "x": 0,
        "y": 0,
        "radius": 0.5,
        "deposits": {
            Names[2]: {
                "shape": "circle",
                "x": 0,
                "y": 0,
                "radius": 0.1   # sits inside Dot_1; Dot_1 surface gets a hole here
            }
        }
    },
    Names[1]: {
        "shape": "rectangle",
        "x": -0.5,
        "y": -0.25,
        "width": 1,
        "height": 0.5
    }
}

```

#### Rotation:

The entire DiMES top surface and dots can be rotated by an angle `theta_dimes` (in degrees) relative to a direction defined by the unit vector components:

- `ax`: The X-component of the rotation axis (float).
- `ay`: The Y-component of the rotation axis (float).
- `az`: The Z-component of the rotation axis (float).

Right now only the following 4 directions are allowed:
    
    ax = ± 1, ay = 0, az = 0
    ax = 0, ay = ± 1, az = 0

if you pick other directions the code will produce an error. If requested, the feature
of rotating around any direction will be added in future releases.

To make a DiMES head with no coatings you can simply set the `no_dots` flag:
    
- `no_dots`: A flag to remove all coatings (bool, default: False).
---

### MESH KEYWORDS:

Before generating the mesh, users can configure several options:

- `msh_dim`: The dimension of the mesh. It can be 1D, 2D, or 3D (int).

- `filename`: The name of the generated mesh file (e.g., "test.msh"; string).

- `save_msh`: A flag indicating whether to save the mesh file.
  - `True`: Save the mesh file.
  - `False`: Do not save the mesh file.

- `GUI_geo`: A flag to visualize the geometry before meshing using a graphical interface.
  - `True`: Use the GUI for geometry visualization.
  - `False`: Do not use the GUI.

- `GUI_msh`: A flag to visualize both the geometry and mesh using a graphical interface.
  - `True`: Use the GUI for both geometry and mesh visualization.
  - `False`: Do not use the GUI.
   """
   
    try:
        gmsh.finalize()
    except:
        pass
    finally:
        gmsh.initialize()

    # Defining keys specific to geometry and mesh
    geo_specific_keys = ['input_dict', 'L_R', 'L_phi', 'L_Z', 'x_center_dimes',
                         'y_center_dimes', 'z_top_dimes', 'r_dimes', 'ax', 'ay', 'az', 'theta_dimes', 'no_dots']
    mesh_specific_keys = ['filename', 'save_msh', 'save_output', 'scale'
                          'GUI_geo', 'GUI_msh', 'msh_dim', 'save_npz', 'save_ply',
                          'ply_normals', 'ply_ascii']

    kw_geo = {key: value for key,
              value in kwargs.items() if key in geo_specific_keys}
    kw_msh = {key: value for key,
              value in kwargs.items() if key in mesh_specific_keys}

    _, component_surfaces = make_dimes_geom(input_dict, **kw_geo)
    make_dimes_mesh(component_surfaces=component_surfaces, **kw_msh)

    try:
        gmsh.finalize()
    except:
        pass

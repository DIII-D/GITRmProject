#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:07:27 2024

@author: cappellil
"""
import gmsh
import math
from .utils import rectangle_def , create_loops

def make_dimes_geom(input_dict, l_radial=4, l_toroidal=4, l_vertical=3,
                    x_center=0, y_center=0, z_center=0,
                    x_center_dimes = 0, y_center_dimes =0 , z_center_dimes=1, r_dimes = 1, 
                    ax=0, ay=1, az=0, theta_dimes=0):
    
    # convert deg to rad
    theta_dimes = math.pi / 180 * theta_dimes 
    
    
    #%% 
    
    """ Plasma volume geometry """
    # Cartesian coordinates of bottom plasma volume surface lower left corner
   
    x_plasma_volume_ll = -l_radial / 2  + x_center# toroidal
    y_plasma_volume_ll = -l_toroidal / 2+ y_center # radial
    z_plasma_volume_ll = z_center# vertical
    
    # Width and height of plasma volume base surface 
    width = l_radial
    height = l_toroidal
    dz_plasma_volume = l_vertical  # This variable is defined but not used in this surface creation
    
    # Create a curve loop for plasma volume base
    p1, p2, p3, p4, l1, l2, l3, l4, plasma_base_rectangle_loop = \
        rectangle_def(x_plasma_volume_ll, y_plasma_volume_ll, z_plasma_volume_ll, width, height)
    
    # Create a curve loop for plasma_volume top
    
    p5, p6, p7, p8, l5, l6, l7, l8, plasma_top_rectangle_loop = rectangle_def(x_plasma_volume_ll, y_plasma_volume_ll, \
                                                              z_plasma_volume_ll + dz_plasma_volume, width, height)


    # Create vertical lines connecting bottom and top
    l9 = gmsh.model.occ.addLine(p1, p5)
    l10 = gmsh.model.occ.addLine(p2, p6)
    l11 = gmsh.model.occ.addLine(p3, p7)
    l12 = gmsh.model.occ.addLine(p4, p8)
    
    # list of surfaces IDs enclosing the plasma volume
    volumes_surfaces = []
    
    # Create surfaces for the sides (lateral surfaces)
    plasma_side1 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l1, l10, -l5, -l9])]) # minus when instead of going for extremity A to B you go from B to A close a loop
    plasma_side2 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l2, l11, -l6, -l10])])
    plasma_side3 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l3, l12, -l7, -l11])])
    plasma_side4 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l4, l9, -l8, -l12])])
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_side1)
    volumes_surfaces.append(plasma_side2)
    volumes_surfaces.append(plasma_side3)
    volumes_surfaces.append(plasma_side4)
    
    # Create the top surface
    plasma_top = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop([l5, l6, l7, l8])])
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_top)
    
    # Synchronize the GMSH model
    gmsh.model.occ.synchronize()
    
    #%% 
    """ DiMES geometry """
    
    DiMES_circle = gmsh.model.occ.addCircle(x_center_dimes , y_center_dimes , z_center , r_dimes)
    
    DiMES_circle_loop = gmsh.model.occ.addCurveLoop([DiMES_circle])
    
    plasma_base = gmsh.model.occ.addPlaneSurface([plasma_base_rectangle_loop, DiMES_circle_loop]) 
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(plasma_base)
    
    if z_center_dimes > z_center: 
        
        r1 = r_dimes / math.cos(theta_dimes) # major radius
        r2 = r_dimes  #minor radius (r1 >= r2)
        ellipse = gmsh.model.occ.addEllipse(x_center_dimes, y_center_dimes, z_center_dimes, r1, r2)
        
        gmsh.model.occ.rotate([(1 , ellipse)], x_center_dimes, y_center_dimes, z_center_dimes, ax, ay, az, theta_dimes)
        ellipse_loop = gmsh.model.occ.addCurveLoop([ellipse])               
    
        # Volumes and surfaces can be constructed from (closed) curve loops thanks to the
        # `addThruSections()' function
    
        DiMES_side_surface = gmsh.model.occ.addThruSections([DiMES_circle_loop, ellipse_loop], makeSolid = False)
      
            
        # store surfaces enclosing volume in a variable
        for tup in DiMES_side_surface:
            if tup[0] == 2:
                DiMES_side_surface_id = tup[1]
                break
        
        volumes_surfaces.append(DiMES_side_surface_id)
        
        #identify top_circle loop
        DiMES_top_curve_loop = gmsh.model.occ.addCurveLoop([ellipse])
        
    else:
        DiMES_top_curve_loop = DiMES_circle_loop
        
    #%%
    
    """ geometry Dots (coatings)"""
        
    dot_loops = []
    
    create_loops(input_dict, z_center_dimes, volumes_surfaces, dot_loops, ax, ay, az, theta_dimes)

    # Generate DiMES top surface
    gmsh.model.occ.synchronize()
    DiMES_top_surface = gmsh.model.occ.addPlaneSurface([DiMES_top_curve_loop] + dot_loops)
    
    # store surfaces enclosing volume in a variable
    volumes_surfaces.append(DiMES_top_surface)
    
    #%% generate the volume
    plasma_volume = gmsh.model.occ.addVolume([gmsh.model.occ.addSurfaceLoop(volumes_surfaces)])
    
    return plasma_volume

#%% 
""" function """
    
def make_dimes_mesh(filename = "test.msh" , save_msh=False, GUI_geo=False, GUI_msh=True , msh_dim=3):
    #%% Generate the mesh and visualize the result
    
    if GUI_geo:
    # Launch the GUI to see the results:
    # Optionally, run the GUI to visualize
        gmsh.fltk.run()
    # Final synchronization of the CAD model
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 1)
    gmsh.option.setNumber("Mesh.MinimumElementsPerTwoPi", 20)
    
    # # Prevent very small elements in small dots
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.05)
    # Set maximum mesh characteristic length for the whole model
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2) 
    gmsh.model.mesh.generate(msh_dim)
    
    if GUI_msh:
    # Launch the GUI to see the results:
    # Optionally, run the GUI to visualize
        gmsh.fltk.run()
        
    if save_msh:
        gmsh.write(filename)
    
    # Finalize GMSH
    gmsh.finalize()
    

#%%
""" function """

def generate_dimes_mesh(input_dict, **kwargs): 
    
    """
This function automatically generates the DiMES mesh for a given set of dot geometries.

### PARAMETERS: (KEYWORDS DESCRIPTION IS AFTER kwargs)
    
- input_dict: A dictionary containing the geometries of the dots to be simulated, where each dot is defined by a unique label followed by its properties.
  
- kwargs: Additional keyword arguments. Supported keywords include:
    
  - x_center_dimes: The X-coordinate of the center of the DiMES base (float, default: 0).
  - y_center_dimes: The Y-coordinate of the center of the DiMES base (float, default: 0).
  - z_center_dimes: The Z-coordinate of the center of the DiMES base (float, default: 0).
  - r_dimes: The radius of the DiMES base (float, default: 1).
  - x_center: The X-coordinate of the plasma volume base surface center (float, default: 0).
  - y_center: The Y-coordinate of the plasma volume base surface center (float, default: 0).
  - z_center: The Z-coordinate of the plasma volume base surface center (float, default: 0).
  - l_radial: The extent of plasma volume in the radial direction (half-length, float, default: 4).
  - l_toroidal: The extent of plasma volume in the toroidal direction (half-length, float, default: 4).
  - l_vertical: The extentof plasma volume  in the vertical direction (float, default: 3).
  - ax: The X-component of the rotation axis (float, default: 0).
  - ay: The Y-component of the rotation axis (float, default: 0).
  - az: The Z-component of the rotation axis (float, default: 0).
  - msh_dim: The dimension of the mesh, which can be 1D, 2D, or 3D (int, default: 3).
  - filename: The name of the mesh file to be generated (string, default: "test.msh").
  - save_msh: A flag indicating whether to save the mesh (bool, default: False).
  - GUI_geo: A flag to use a GUI for visualizing the geometry before meshing (bool, default: False).
  - GUI_msh: A flag to use a GUI for visualizing both geometry and mesh (bool, default: True).
  
### GEOMETRY KEYWORDS:

The DiMES geometry can take one of the following shapes:
    
- disk
- cylinder
- truncated cylinder

The geometry is defined by its center and radius:
    
- `x_center_dimes`: The X-coordinate of the center of the DiMES base (float).
- `y_center_dimes`: The Y-coordinate of the center of the DiMES base (float).
- `z_center_dimes`: The Z-coordinate of the center of the DiMES base (float).
- `r_dimes`: The radius of the DiMES base (float).

#### Plasma Volume:
    
The DiMES geometry is enclosed within a volume known as `plasma_volume`, representing the plasma background. This volume is a symmetrical box extending in three directions:

- `x_center`: The X-coordinate of the plasma volume center (float).
- `y_center`: The Y-coordinate of the plasma volume center (float).
- `z_center`: The Z-coordinate of the plasma volume center (float).

The plasma volume extends symmetrically in these directions:
    
- `l_radial`: Twice the radial extent of the volume (float).
- `l_toroidal`: Twice the toroidal extent of the volume (float).
- `l_vertical`: The vertical extent of the volume (float).

#### Dots (Material Coatings):
    
On the DiMES top surface, dots represent the material coatings. Two geometries are currently supported:
    
- circle
- rectangle

Each dot's position and dimensions are specified within an `input_dict` dictionary. The user must manually create this dictionary, where each dot is defined by a unique label followed by values that describe its geometry, shape, and position.

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
    },
    Names[1]: {
        "shape": "rectangle",
        "x": -0.5,
        "y": -0.25,
        "width": 1,    # Width of the rectangle
        "height": 0.5  # Height of the rectangle
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

#### Rotation:
    
The entire DiMES top surface and dots can be rotated by an angle `theta_dimes` (in degrees) relative to a direction defined by the unit vector components:

- `ax`: The X-component of the rotation axis (float).
- `ay`: The Y-component of the rotation axis (float).
- `az`: The Z-component of the rotation axis (float).

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
    geo_specific_keys = ['input_dict', 'l_radial', 'l_toroidal','l_vertical','x_center_dimes',
                         'y_center_dimes','z_center_dimes','r_dimes' ,'ax','ay', 'az', 'theta_dimes']
    mesh_specific_keys = ['filename', 'save_msh', 'GUI_geo', 'GUI_msh', 'msh_dim']    
    
    kw_geo = {key: value for key, value in kwargs.items() if key in geo_specific_keys}
    kw_msh = {key: value for key, value in kwargs.items() if key in mesh_specific_keys}
    
    make_dimes_geom(input_dict, **kw_geo)
    make_dimes_mesh(**kw_msh)
    
    try:
        gmsh.finalize()
    except:
        pass
    
help(generate_dimes_mesh)
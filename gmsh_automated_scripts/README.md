# GITRmProject - gmsh_automated_scripts

Hi User!
As part of the effort to automate the workflow of DiMES local impurity transport simulations at DIII-D using the GITRm transport code,
we propose a tool that generates meshes for simple DiMES geometries using the gmsh Python API. 
This user-friendly tool allows anyone to easily simulate their preferred DiMES setup without 
the need to manually draw and mesh the geometry each time.

# Dependecies

This package was tested using pyhton version >= 3.10, gmsh version 4.13.1 and requires the math package

# Build

This package was built using poetry , if you don't know what poetry is please check out: https://python-poetry.org/. 
Also, a short tutorial on how to use poetry with conda or mamba is provided in the wiki: https://github.com/DIII-D/GITRmProject/wiki/Use-gmsh-scripts--with-poetry

# Documentation 

The package consists of two files encoding all functions needed to run **generate_dimes_mesh(input_dict,**kwargs)**

**generate_dimes_mesh(input_dict,**kwargs)** automatically generates the DiMES mesh for a given input.

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

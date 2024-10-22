#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 16:06:17 2024

@author: cappellil
"""
import gmsh
import math

#%%
""" miscellaneous functions """

def rectangle_def(x, y, z, width, height):
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

#%%
def create_loops(input_dict, z_dimes, volumes_surfaces, dot_loops, ax, ay, az, theta_dimes):
    
    # check if user is not tilting both DiMES head and dots
    
    if theta_dimes != 0:
    
        for elem_def in input_dict.values():
    # Check if theta_dot exists and is not zero
            if "theta_dot" in elem_def and elem_def["theta_dot"] != 0:
                print(f"{elem_def}: theta_dot is not zero (theta_dot = {elem_def['theta_dot']})")
                raise Exception("You are setting a tilted dot (theta_dot != 0) " \
                      "on a tilted DiMES head (theta_dimes !=0), this could be problematic")
                    
                    
    for elem_def in input_dict.values():
    
        if elem_def["shape"] == "circle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            r = elem_def['radius']
            theta_dot =  theta_dimes + math.pi / 180 * elem_def['theta_dot']
            
            # z_dimes_tilt indicates the z position of dot over a tilted DiMES
            
            z_dimes_tilt = math.tan(theta_dimes) * (ay * x + ax * y)
            
            z += z_dimes_tilt
            
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
            dot_loops.append(dot_base_loop)
            
            
            if theta_dot > theta_dimes:
                
                
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
                
            else:
                
                # if dot is not tilted, directly use the base to make a disk
                dot_surface = gmsh.model.occ.addPlaneSurface([dot_base_loop])
                
                # plane has to be coplanar to DiMES head
                gmsh.model.occ.rotate([(2 , dot_surface)], x, y, z, ax, ay, az, theta_dimes)
                
                # append dot_base_loop to list of holes composing DiMES head surface
                dot_loops.append(dot_base_loop)
                
                # append disk surfaces to list of surface delimiting the plasma volume
                volumes_surfaces.append(dot_surface)
            
                    
        elif elem_def["shape"] == "rectangle":
            
            x = elem_def['x']
            y = elem_def['y']
            z = z_dimes
            width = elem_def['width']
            height = elem_def['height']
            theta_dot = theta_dimes +  math.pi / 180 * elem_def['theta_dot']
            
            # if dot coating surface is not tilted, dot simulated as a Rectangular surface
            # coplanar with the DiMES head
            
            # otherwise dot simulated as a wedge delimitated by a tilted plane at the top
            # and a plane coplanar to DiMES head at the base
            
            # z_dimes_tilt indicates the z position of dot over a tilted DiMES
            
            z_dimes_tilt = math.tan(theta_dimes) * (ay * x + ax * y) 
            
            z += z_dimes_tilt
            
            # translate base rectangle to be coplanar with DiMES head and to avoid intersections
            delta_z = 0.001 # to avoid overlapping between curves
            delta_z_dimes = ay * width / 2 * math.tan(theta_dimes)
            
            # 1. Create base lines and curve loop shifted along z
            base_l1, base_l2, base_l3, base_l4, dot_base_loop = rectangle_def(x, y, z + delta_z_dimes, width , height)[-5:]

            # 2. rotate rectangle using base_l1 edge as pivotal point around y-axis
            #    and create loop
            
            gmsh.model.occ.rotate([(1, base_l1), (1, base_l2), (1, base_l3), (1, base_l4)], x, y, z + delta_z_dimes, ax, ay, az, theta_dimes)
            dot_base_loop = gmsh.model.occ.addCurveLoop([base_l1, base_l2, base_l3, base_l4])
            
            # append dot_base_loop to list of holes composing DiMES head surface
            dot_loops.append(dot_base_loop)

            if theta_dot != 0:
                
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
                
                z_dimes_tilt = math.tan(theta_dimes) * (ay * x_c + ax * y_c)
                z += z_dimes_tilt
                
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
                l1, l2, l3, l4, dot_loop = rectangle_def(x, y, z + delta_z_dimes + delta_z, w, height)[-5:]
        
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
                    
                # Synchronize the model to update geometry
                gmsh.model.occ.synchronize()
            else:
                # create base surface
                dot_base_surface = gmsh.model.occ.addPlaneSurface([dot_base_loop])
                
                # append disk surfaces to list of surface delimiting the plasma volume
                volumes_surfaces.append(dot_base_surface)
            
            
        

import numpy as np
import gmsh
import math
from .data_structures import Object2D, AngledSample

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


def make_angled_sample(sample: Object2D):

    assert isinstance(sample.shape, AngledSample), \
        TypeError(f"Sample shape is of type {type(sample.shape)} while it must be of type AngledSample.")

    occ = gmsh.model.occ
    shape = sample.shape

    dz     = shape.r * np.tan(shape.angle)
    z_flat = shape.height + shape.z_cut

    assert shape.z_cut < 2 * dz - 1e-10, \
        f"z_cut={shape.z_cut} too large! max allowed: {2*dz:.6f}"

    margin = 0.1 * shape.r

    # 1. Cylinder
    cyl_bottom = shape.height - dz - margin
    cyl_top    = shape.height + dz + shape.z_cut + margin
    cyl = occ.addCylinder(0, 0, cyl_bottom, 0, 0, cyl_top - cyl_bottom, shape.r)

    # 2. Tilted box cut
    box_tilt = occ.addBox(-2*shape.r, -2*shape.r, shape.height,
                           4*shape.r,  4*shape.r, dz + shape.z_cut + margin)
    occ.rotate([(3, box_tilt)], 0, 0, shape.height, 0, 1, 0, -shape.angle) # rotate around Y-axis

    # 3. Flat box cut
    box_flat = occ.addBox(-2*shape.r, -2*shape.r, z_flat,
                           4*shape.r,  4*shape.r, dz + margin)

    # 4. Cut
    truncated, _ = occ.cut(
        [(3, cyl)],
        [(3, box_tilt), (3, box_flat)],
        removeObject=True, removeTool=True
    )
    occ.synchronize()

# 5. Identify bottom face BEFORE translation:
    #    must be FLAT (zmin ≈ zmax) AND lowest — the side wall also touches
    #    the bottom z, so filtering by min(zmin) alone is ambiguous.
    faces = gmsh.model.getBoundary(truncated, oriented=False)

    z_tol = min(dz, shape.z_cut) * 1e-2   # separates flat faces from spanning ones

    face_info = {}
    for dim, tag in faces:
        _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(dim, tag)
        face_info[tag] = (zmin, zmax, dim)

    # flat faces only (near-zero z-span)
    flat_faces = {t: v for t, v in face_info.items() if abs(v[1] - v[0]) < z_tol}
    bottom_tag = min(flat_faces, key=lambda t: flat_faces[t][0])
    print(f"  Bottom face (flat, lowest): tag={bottom_tag}, "
          f"zmin={face_info[bottom_tag][0]:.6f}, zmax={face_info[bottom_tag][1]:.6f}")

    faces_to_keep = [(dim, tag) for dim, tag in faces if tag != bottom_tag]
    for dim, tag in faces_to_keep:
        print(f"  Keeping face {tag}: zmin={face_info[tag][0]:.6f}, zmax={face_info[tag][1]:.6f}")

    # 6. Translate so zmin of solid = 0
    vol_tag = truncated[0][1]
    _, _, zmin_solid, _, _, _ = gmsh.model.getBoundingBox(3, vol_tag)
    occ.translate(truncated, 0, 0, -zmin_solid)
    occ.synchronize()

    # 7. Copy kept surfaces and remove solid
    surfaces = occ.copy(faces_to_keep)
    occ.remove(truncated, recursive=True)
    occ.synchronize()

    # 8. Translate to final (x, y) location
    occ.translate(surfaces, *shape.center)
    occ.synchronize()

    # 9. Build the hole loop explicitly — we already KNOW it's a circle of
    #    radius shape.r at the sample center, sitting at z=0. No fishing.
    cx, cy, _ = shape.center
    bottom_circle = occ.addCircle(cx, cy, 0.0, shape.r)
    bottom_loop   = occ.addCurveLoop([bottom_circle])
    occ.synchronize()

    surface_tags = [tag for _, tag in surfaces]
    return bottom_loop, surface_tags
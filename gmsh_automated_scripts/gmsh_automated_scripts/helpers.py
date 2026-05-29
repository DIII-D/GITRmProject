import math

def z_on_tilted_surface(x, y, z, ax, ay, theta):
    delta_z = - math.tan(theta) * (ay * x + ax * y)
    return z + delta_z

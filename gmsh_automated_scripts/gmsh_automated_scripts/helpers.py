import math

def z_on_tilted_surface(z0, x, y, ax, ay, theta):
    delta_z = - math.tan(theta) * (ay * x + ax * y)
    return z0 + delta_z

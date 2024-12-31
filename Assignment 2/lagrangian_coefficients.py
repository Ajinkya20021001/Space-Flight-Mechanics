# Step 6 of Algorithm 5.1 (Lambert's problem) from Curtis

from stumpff_functions import C, S
import math

def lagrangianCoefficients(r1, r2, z, A):
    y = r1 + r2 + A*(z*S(z) - 1)/math.sqrt(C(z))

    mu = 398600.4418

    f = 1 - y/r1
    g = A * math.sqrt(y / mu)
    f_dot = (math.sqrt(mu) / (r1 * r2)) * math.sqrt(y / C(z)) * (z*S(z) - 1)
    g_dot = 1 - y/r2

    return f, g, f_dot, g_dot
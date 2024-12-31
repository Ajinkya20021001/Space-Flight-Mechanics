import numpy as np


def lambertProblem(R1, R2, t):
    # Step 1: Calculate r1 and r2 using Equation 5.24.
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    # Step 2:
    # a. Choose either prograde or retrograde motion.
    isPrograde = True

    # b. Calculate deltaTheta using Equation 5.26.
    from delta_theta import deltaTheta
    deltaTheta = deltaTheta(R1, R2, isPrograde)

    # Step 3: Calculate A using Equation 5.35.
    A = np.sin(deltaTheta) * np.sqrt((r1 * r2) / (1 - np.cos(deltaTheta)))

    # Step 4: Solve for z using Equation 5.40, 5.43, and 5.45.
    from newton_raphson import newtonRaphson
    z = newtonRaphson(t, A, r1, r2)
    print(f"Z: {z}")

    # # (Redundant) Step 5: Calculate y using Equation 5.38.
    # from y import y
    # y = y(R1, R2, z, A)

    # Step 6: Calculate f, g, f_dot, and g_dot using Equation 5.46.
    from lagrangian_coefficients import lagrangianCoefficients
    f, g, f_dot, g_dot = lagrangianCoefficients(r1, r2, z, A)

    # Step 7: Calculate the velocity vectors at R1 and R2 using Equation 5.28 and 5.29.
    V1 = (R2 - f*R1)/g
    V2 = ((g_dot/g) * R2) - ((((f*g_dot) - (f_dot*g))/g) * R1)

    return V1, V2





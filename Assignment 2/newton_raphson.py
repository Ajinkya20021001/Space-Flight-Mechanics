# 4. By iteration, using Equations 5.40, 5.43 and 5.45, solve Equation 5.39 for z.
# The sign of z tells us whether the orbit is a hyperbola (z < 0), parabola (z = 0) or ellipse (z > 0).

import numpy as np
import math

from stumpff_functions import C, S

mu = 398600.4418

def newtonRaphson(delta_t, A, r1, r2):
    def F(z):
        F_z = (((y(z)/C(z))**1.5) * S(z)) + (A * math.sqrt(y(z))) - (math.sqrt(mu) * delta_t)
        return F_z

    def F_prime(z): 
        if (z!= 0):
            a = ((y(z))/(C(z))) ** 1.5
            b1 = C(z)
            b2 = 1.5 * (S(z)/C(z))  
            b = (0.5/z) * (b1 - b2)
            c = 0.75 * (S(z)**2)/(C(z))

            term1 = a * (b + c)
            term2 = (A/8)*((3*(S(z)/C(z))*math.sqrt(y(z))) + (A * math.sqrt(C(z)/y(z))))
            F_prime_z = term1 + term2
        else:
            F_prime_z = ((math.sqrt(2)/40) * ((y(0))**1.5)) + (A/8) * (math.sqrt(y(0)) + A*math.sqrt(1/(2*y(0))))
        return F_prime_z

    def y(z):
        y_z = r1 + r2 + A * ((z * S(z)) - 1) / math.sqrt(C(z))
        return y_z  
        
    # for the initial guess of z, find where F(z) changes sign
    z = -2
    while (F(z) < 0):
        z += 0.1     

    # Newton Raphson iteration:
    
    maxIterations = 1000
    convergenceCriterion = 1e-8
    iteration = 1
    ratio = 1
    
    while (iteration < maxIterations) and (abs(ratio) > convergenceCriterion):
        ratio = F(z)/F_prime(z)
        z -= ratio
        iteration += 1
    
    return z

if __name__ == '__main__':
    # Test the function
    R1 = np.array([5000, 10000, 2100])      # km
    R2 = np.array([-14600, 2500, 7000])     # km
    deltaT = 60 * 60                        # 60 minutes, but in seconds
    A = 12372
    z = newtonRaphson(deltaT, A, np.linalg.norm(R1), np.linalg.norm(R2))
    print(z)
    # print(y(z))
    # print(C(z))
    # print(S(z))
    # print(F(z))
    # print(F_prime(z))
    # print(num_iter)

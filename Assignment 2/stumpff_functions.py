from math import sqrt, cos, sin, cosh, sinh

def C(z):
        if z > 0:
            return (1-cos(sqrt(z)))/(z)
        elif z < 0:
            return (cosh(sqrt(-z)) - 1)/(-z)
        else:
            return 1/2
        
def S(z):
        if z > 0:
            return (sqrt(z)-sin(sqrt(z)))/(sqrt(z)**3)
        elif z < 0:
            return (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)**3)
        else:
            return 1/6
# Step 2 of Algorithm 5.1 (Lambert's problem) from Curtis

import numpy as np

def deltaTheta(R1, R2, isPrograde):

    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    if isPrograde:
        if np.cross(R1, R2)[2] >= 0:
            dTheta = np.arccos(np.dot(R1, R2) / (r1 * r2))
        else:
            dTheta = 2 * np.pi - np.arccos(np.dot(R1, R2) / (r1 * r2))
    else:
        if np.cross(R1, R2)[2] >= 0:
            dTheta = 2 * np.pi - np.arccos(np.dot(R1, R2) / (r1 * r2))
        else:
            dTheta = np.arccos(np.dot(R1, R2) / (r1 * r2))

    return dTheta


if __name__ == '__main__':
    # Test the function

    R1 = np.array([5000, 10000, 2100])      # km
    R2 = np.array([-14600, 2500, 7000])     # km

    deltaTheta = deltaTheta(R1, R2, True)
    print(np.degrees(deltaTheta))

    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    A = np.sin(deltaTheta) * np.sqrt((r1 * r2) / (1 - np.cos(deltaTheta)))
    print(A)
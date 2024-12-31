import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_orbit(a, e, i, RAAN, arg_periapsis, true_anomaly):
    # Constants
    mu = 398600  # Earth's gravitational parameter, km^3/s^2
    earth_radius = 6371  # Earth's radius in km

    # Calculate the semi-latus rectum
    p = abs(a) * (1 - e**2)

    # Generate true anomaly values
    if e < 1:
        theta = np.linspace(0, 2 * np.pi, 1000)
    else:  # Hyperbolic orbit
        theta = np.linspace(-np.radians(60), np.radians(60), 1000)

    # Calculate the orbit in the perifocal coordinate system
    r = p / (1 + e * np.cos(theta))
    x_perifocal = r * np.cos(theta)
    y_perifocal = r * np.sin(theta)
    z_perifocal = np.zeros_like(theta)

    # Rotation matrices
    R3_W = np.array([[np.cos(RAAN), -np.sin(RAAN), 0],
                     [np.sin(RAAN),  np.cos(RAAN), 0],
                     [0,             0,            1]])

    R1_i = np.array([[1, 0,             0],
                     [0, np.cos(i), -np.sin(i)],
                     [0, np.sin(i),  np.cos(i)]])

    R3_w = np.array([[np.cos(arg_periapsis), -np.sin(arg_periapsis), 0],
                     [np.sin(arg_periapsis),  np.cos(arg_periapsis), 0],
                     [0,                      0,                     1]])

    # Combined rotation matrix
    Q = R3_W @ R1_i @ R3_w

    # Rotate the perifocal coordinates to the geocentric equatorial frame
    r_geocentric = Q @ np.vstack((x_perifocal, y_perifocal, z_perifocal))

    # Plot the orbit
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(r_geocentric[0], r_geocentric[1], r_geocentric[2], label='Orbit', color='r')

    # Plot the Earth
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:40j]
    x = earth_radius * np.cos(u) * np.sin(v)
    y = earth_radius * np.sin(u) * np.sin(v)
    z = earth_radius * np.cos(v)
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.6)

    # Set the x, y, z limits
    ax.set_xlim([-40000, 40000])
    ax.set_ylim([-40000, 40000])
    ax.set_zlim([-40000, 40000])

    # Set the default view to front view
    ax.view_init(elev=30, azim=135)

    # Set the aspect ratio to be equal
    ax.set_aspect('auto')

    # Labels and legend
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()

    plt.show()

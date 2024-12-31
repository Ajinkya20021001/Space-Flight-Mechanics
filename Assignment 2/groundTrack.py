import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ...existing code...

def plot_ground_track(orbital_elements):
    a = orbital_elements['semi_major_axis']     # Semi-major axis, km
    e = orbital_elements['eccentricity']        # Eccentricity
    i = orbital_elements['inclination']         # Inclination, radians
    RAAN = orbital_elements['RAAN']             # Right Ascension of Ascending Node, radians
    w = orbital_elements['argument_of_perigee'] # Argument of perigee, radians
    M0 = orbital_elements['mean_anomaly']       # Mean anomaly at epoch, radians
    t0 = orbital_elements['epoch']              # Epoch time, seconds

    # Constants
    mu = 398600.4418  # Earth's gravitational parameter, km^3/s^2
    Re = 6378.137  # Earth's radius, km

    # Time array
    t = np.linspace(0, 86400, 1000)  # 1 day period

    if a > 0:  # Elliptical orbit
        # Mean motion for elliptical orbit
        n = np.sqrt(mu / a**3)
        M = M0 + n * (t - t0)
        E = M
        for _ in range(10):  # Solve Kepler's equation
            E = M + e * np.sin(E)
        r = a * (1 - e * np.cos(E))
        nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    elif a < 0:  # Hyperbolic orbit
        p = abs(a) * (1 - e**2)  # Semi-latus rectum
        F = np.linspace(-3, 3, 1000)  # Hyperbolic anomaly
        r = p / (1 + e * np.cosh(F))  # Distance in orbit plane
        nu = 2 * np.arctan(np.sqrt((e + 1) / (e - 1)) * np.tanh(F / 2))
    else:
        raise ValueError("Semi-major axis (a) cannot be zero for a valid orbit.")

    # Position in orbital plane
    x_orb = r * np.cos(nu)
    y_orb = r * np.sin(nu)

    # Rotation matrices
    R3_W = np.array([[np.cos(RAAN), -np.sin(RAAN), 0],
                     [np.sin(RAAN), np.cos(RAAN), 0],
                     [0, 0, 1]])

    R1_i = np.array([[1, 0, 0],
                     [0, np.cos(i), -np.sin(i)],
                     [0, np.sin(i), np.cos(i)]])

    R3_w = np.array([[np.cos(w), -np.sin(w), 0],
                     [np.sin(w), np.cos(w), 0],
                     [0, 0, 1]])

    # Position in inertial frame
    r_orb = np.vstack((x_orb, y_orb, np.zeros_like(x_orb)))
    r_eci = R3_W @ R1_i @ R3_w @ r_orb

    # Convert to latitude and longitude
    lon = np.arctan2(r_eci[1], r_eci[0])
    lat = np.arcsin(r_eci[2] / np.linalg.norm(r_eci, axis=0))

    # Convert to degrees
    lon = np.degrees(lon)
    lat = np.degrees(lat)

    # Plotting
    plt.figure(figsize=(8, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = True
    gl.right_labels = True

    ax.plot(lon, lat, 'b', transform=ccrs.Geodetic())

    # Show the full map irrespective of the ground track
    ax.set_global()

    plt.show()


import numpy as np
import matplotlib.pyplot as plt
from lambert_problem import lambertProblem
from coe_from_sv import coe_from_sv, print_orbital_elements
from traceOrbit import plot_orbit
from groundTrack import plot_ground_track

# Input data
R1 = np.array([5644, -2830, 4170])      # km
R2 = np.array([-2240, 7320, -4980])     # km
deltaT = 20 * 60                        # 20 minutes, but in seconds

# Solve Lambert's problem for velocity vectors
V1, V2 = lambertProblem(R1, R2, deltaT)

# Convert state vectors to orbital elements
orbitalElements = coe_from_sv(R1, V1)

# Debugging: Print orbital elements structure
print("Orbital Elements:", orbitalElements)

# Extract all orbital elements including angular momentum and semi-major axis
h, e, i, RAAN, arg_periapsis, true_anomaly, a = orbitalElements

# Convert angles from radians to degrees for clarity (optional)
i_deg = np.degrees(i)
RAAN_deg = np.degrees(RAAN)
arg_periapsis_deg = np.degrees(arg_periapsis)
true_anomaly_deg = np.degrees(true_anomaly)

# Display elements for clarity
print(f"Angular momentum (h): {h}")
print(f"Eccentricity (e): {e}")
print(f"Inclination (i): {i_deg} degrees")
print(f"RAAN: {RAAN_deg} degrees")
print(f"Argument of Periapsis: {arg_periapsis_deg} degrees")
print(f"True Anomaly: {true_anomaly_deg} degrees")
print(f"Semi-major Axis (a): {a} km")

# Plot the orbit (even for hyperbolic trajectories)
plot_orbit(a, e, i, RAAN, arg_periapsis, true_anomaly)

# Prepare orbital elements dictionary for ground track
orbital_elements_dict = {
    'semi_major_axis': a,
    'eccentricity': e,
    'inclination': i,  # Already in radians
    'RAAN': RAAN,  # Already in radians
    'argument_of_perigee': arg_periapsis,  # Already in radians
    'mean_anomaly': true_anomaly,  # Use true anomaly as an approximation for mean anomaly
    'epoch': 0  # Assuming epoch is at time zero
}

# Plot the ground track
plot_ground_track(orbital_elements_dict)

import numpy as np

# Constants
G = 6.67430e-11  # Gravitational constant, m^3/kg/s^2
M_Earth = 5.972e24  # Mass of Earth, kg 
M_Moon = 7.34767309e22  # Mass of Moon, kg
R_Earth = 6371e3  # Radius of Earth, m
R_Moon = 1737e3  # Radius of Moon, m
d_Earth_to_Moon = 384400e3  # Average distance from Earth to Moon, m
mu_Earth = G * M_Earth  # Gravitational parameter of Earth, m^3/s^2
mu_Moon = G * M_Moon  # Gravitational parameter of Moon, m^3/s^2

# Inputs
r_parking_earth = R_Earth + 200e3  # Earth parking orbit altitude, m
r_transfer_apogee = d_Earth_to_Moon  # Apogee of transfer orbit, m
r_parking_moon = R_Moon + 100e3  # Moon parking orbit altitude, m

# Step 1: Velocities in Earth's Sphere of Influence
v_circular_earth = np.sqrt(mu_Earth / r_parking_earth)  # Circular orbit velocity around Earth
v_transfer_perigee = np.sqrt(mu_Earth * (2 / r_parking_earth - 1 / ((r_parking_earth + r_transfer_apogee) / 2)))  # Velocity at perigee of transfer orbit
delta_v1 = v_transfer_perigee - v_circular_earth  # Burn to enter transfer orbit

# Step 2: Velocities at Lunar Sphere of Influence
v_transfer_apogee = np.sqrt(mu_Earth * (2 / r_transfer_apogee - 1 / ((r_parking_earth + r_transfer_apogee) / 2)))  # Velocity at apogee of transfer orbit (Earth SOI)
v_circular_moon = np.sqrt(mu_Moon / r_parking_moon)  # Circular orbit velocity around Moon
v_moon_influence_entry = v_transfer_apogee  # Velocity entering Moon's SOI (approximated)
delta_v2 = np.abs(v_circular_moon - v_moon_influence_entry)  # Burn to circularize around Moon

# Total Delta-V
delta_v_total = delta_v1 + delta_v2

# Orbital elements of transfer orbit
semi_major_axis = (r_parking_earth + r_transfer_apogee) / 2  # Semi-major axis, m
eccentricity = (r_transfer_apogee - r_parking_earth) / (r_parking_earth + r_transfer_apogee)  # Eccentricity

# Time of transfer
orbital_period = 2 * np.pi * np.sqrt(semi_major_axis**3 / mu_Earth)  # Full orbital period, s
transfer_time = orbital_period / 2  # Time for half the orbit, s
transfer_time_hours = transfer_time / 3600  # Convert time to hours

# Position vectors
position_vector_first_burn = [r_parking_earth, 0, 0]  # Position vector at first burn (Earth parking orbit)
position_vector_second_burn = [r_transfer_apogee, 0, 0]  # Position vector at second burn (apogee of transfer orbit)

# Results
print(f"Delta-V for the mission:")
print(f"  First burn (Earth parking orbit to transfer orbit): {delta_v1:.2f} m/s")
print(f"  Second burn (Moon SOI to Moon parking orbit): {delta_v2:.2f} m/s")
print(f"  Total Delta-V: {delta_v_total:.2f} m/s\n")

print("Orbital elements of the transfer orbit:")
print(f"  Semi-major axis (a): {semi_major_axis:.2f} m")
print(f"  Eccentricity (e): {eccentricity:.4f}")
print(f"  Perigee radius (r_perigee): {r_parking_earth:.2f} m")
print(f"  Apogee radius (r_apogee): {r_transfer_apogee:.2f} m\n")

print(f"Time required to complete the transfer: {transfer_time:.2f} seconds ({transfer_time_hours:.2f} hours)\n")

print("Position Vectors:")
print(f"  Position vector at first burn (Earth parking orbit): {position_vector_first_burn} m")
print(f"  Position vector at second burn (transfer orbit apogee): {position_vector_second_burn} m")






# Solution
# Delta-V for the mission:
#   First burn (Earth parking orbit to transfer orbit): 3133.10 m/s
#   Second burn (Moon SOI to Moon parking orbit): 1447.20 m/s
#   Total Delta-V: 4580.30 m/s

# Orbital elements of the transfer orbit:
#   Semi-major axis (a): 195485500.00 m
#   Eccentricity (e): 0.9664
#   Perigee radius (r_perigee): 6571000.00 m
#   Apogee radius (r_apogee): 384400000.00 m

# Time required to complete the transfer: 430089.59 seconds (119.47 hours)

# Position Vectors:
#   Position vector at first burn (Earth parking orbit): [6571000.0, 0, 0] m
#   Position vector at second burn (transfer orbit apogee): [384400000.0, 0, 0] m
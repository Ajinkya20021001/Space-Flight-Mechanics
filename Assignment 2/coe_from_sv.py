# Step 8 of Algorithm 5.1 (Lambert's problem) from Curtis

import numpy as np

def coe_from_sv(position_vector, velocity_vector):
    """Computes the classical orbital elements (coe) from the state vector (position and velocity).
    
    Parameters:
        position_vector (numpy array): Position vector in the geocentric equatorial frame (km)
        velocity_vector (numpy array): Velocity vector in the geocentric equatorial frame (km/s)
    
    Returns:
        list: Orbital elements [angular momentum, eccentricity, right ascension, inclination, argument of perigee, true anomaly, semi-major axis]

    Example:
        >>> position_vector = np.array([-6045, -3490, 2500])
        >>> velocity_vector = np.array([-3.457, 6.618, 2.533])
        >>> orbital_elements = coe_from_sv(position_vector, velocity_vector)
    """
    gravitational_parameter = 398600  # gravitational parameter (km^3/s^2)
    epsilon = 1.e-10
    position_magnitude = np.linalg.norm(position_vector)
    velocity_magnitude = np.linalg.norm(velocity_vector)
    radial_velocity = np.dot(position_vector, velocity_vector) / position_magnitude
    angular_momentum_vector = np.cross(position_vector, velocity_vector)
    angular_momentum_magnitude = np.linalg.norm(angular_momentum_vector)
    
    # Using Equation 4.7 from Curtis:
    inclination = np.arccos(angular_momentum_vector[2] / angular_momentum_magnitude)
    
    # Using Equation 4.8 from Curtis:
    node_line_vector = np.cross([0, 0, 1], angular_momentum_vector)
    node_line_magnitude = np.linalg.norm(node_line_vector)
    
    # Using Equation 4.9 from Curtis:
    if node_line_magnitude != 0:
        right_ascension = np.arccos(node_line_vector[0] / node_line_magnitude)
        if node_line_vector[1] < 0:
            right_ascension = 2 * np.pi - right_ascension
    else:
        right_ascension = 0
    
    # Using Equation 4.10 from Curtis:
    eccentricity_vector = (1 / gravitational_parameter) * ((velocity_magnitude**2 - gravitational_parameter / position_magnitude) * position_vector - position_magnitude * radial_velocity * velocity_vector)
    eccentricity = np.linalg.norm(eccentricity_vector)
    
    # Equation 4.12 (incorporating the case e=0):
    if node_line_magnitude != 0:
        if eccentricity > epsilon:
            argument_of_perigee = np.arccos(np.dot(node_line_vector, eccentricity_vector) / (node_line_magnitude * eccentricity))
            if eccentricity_vector[2] < 0:
                argument_of_perigee = 2 * np.pi - argument_of_perigee
        else:
            argument_of_perigee = 0
    else:
        argument_of_perigee = 0
    
    # Equation 4.13a (incorporating the case e=0):
    if eccentricity > epsilon:
        true_anomaly = np.arccos(np.dot(eccentricity_vector, position_vector) / (eccentricity * position_magnitude))
        if radial_velocity < 0:
            true_anomaly = 2 * np.pi - true_anomaly
    else:
        cross_product = np.cross(node_line_vector, position_vector)
        if cross_product[2] >= 0:
            true_anomaly = np.arccos(np.dot(node_line_vector, position_vector) / (node_line_magnitude * position_magnitude))
        else:
            true_anomaly = 2 * np.pi - np.arccos(np.dot(node_line_vector, position_vector) / (node_line_magnitude * position_magnitude))
    
    # Equation 2.61 (a < 0 for a hyperbola):
    semi_major_axis = angular_momentum_magnitude**2 / gravitational_parameter / (1 - eccentricity**2)
    orbital_elements = [angular_momentum_magnitude, eccentricity, right_ascension, inclination, argument_of_perigee, true_anomaly, semi_major_axis]
    
    return orbital_elements

def print_orbital_elements(orbital_elements):
    """Prints the orbital elements.
    
    Parameters:
        orbital_elements (list): Orbital elements [angular momentum, eccentricity, right ascension, inclination, argument of perigee, true anomaly, semi-major axis]
    """
    degrees_conversion_factor = np.pi / 180
    print('Orbital Elements:')
    print(f'Angular momentum (km^2/s) = {orbital_elements[0]}')
    print(f'Eccentricity = {orbital_elements[1]}')
    print(f'Right ascension (deg) = {orbital_elements[2] / degrees_conversion_factor}')
    print(f'Inclination (deg) = {orbital_elements[3] / degrees_conversion_factor}')
    print(f'Argument of perigee (deg) = {orbital_elements[4] / degrees_conversion_factor}')
    print(f'True anomaly (deg) = {orbital_elements[5] / degrees_conversion_factor}')
    print(f'Semimajor axis (km) = {orbital_elements[6]}')
    if orbital_elements[1] < 1:
        gravitational_parameter = 398600  # gravitational parameter (km^3/s^2)
        period = 2 * np.pi / np.sqrt(gravitational_parameter) * orbital_elements[6]**1.5  # Equation 2.73
        print('Period:')
        print(f'Seconds = {period}')
        print(f'Minutes = {period / 60}')
        print(f'Hours = {period / 3600}')
        print(f'Days = {period / 24 / 3600}')

def main():
    gravitational_parameter = 398600  # gravitational parameter (km^3/s^2)
    degrees_conversion_factor = np.pi / 180
    
    # Input data:
    position_vector = np.array([-6045, -3490, 2500])
    velocity_vector = np.array([-3.457, 6.618, 2.533])
    
    # Algorithm 4.1:
    orbital_elements = coe_from_sv(position_vector, velocity_vector)
    
    # Echo the input data and output results to the command window:
    print('---------------------------------------------------')
    print('Example 4.3')
    print(f'Gravitational parameter (km^3/s^2) = {gravitational_parameter}')
    print('State vector:')
    print(f'Position vector (km) = [{position_vector[0]} {position_vector[1]} {position_vector[2]}]')
    print(f'Velocity vector (km/s) = [{velocity_vector[0]} {velocity_vector[1]} {velocity_vector[2]}]')
    print_orbital_elements(orbital_elements)

if __name__ == "__main__":
    main()
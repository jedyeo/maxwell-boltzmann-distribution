# Maxwell-Boltzmann Distribution for air and helium

import math
import numpy as np
import matplotlib.pyplot as plt

# Parameters (to change or add)
molar_mass_air = 0.029  # kg/mol for air
molar_mass_he2 = 0.004  # kg/mol for He2

# Constants
GAMMA_AIR = 1.4
GAMMA_HE2 = 1.66
MOLECULES_IN_MOL = 6.022E23
R = 8.314
KB = 1.38064852E-23  # Boltzmann constant
MASS_AIR = molar_mass_air / MOLECULES_IN_MOL  # mass per molecule of air
MASS_HE2 = molar_mass_he2 / MOLECULES_IN_MOL  # mass per molecule of He2


def mb_dist(m, v, T):
    return (math.sqrt(pow(m / (2 * math.pi * KB * T), 3))) * (4 * math.pi * pow(v, 2)) * np.exp((-m * pow(v, 2)) / (2 * KB * T))


def get_rms_velocity(m, T):
    return math.sqrt(3 * KB * T / m)


def get_speed_of_sound(gamma, m, T):
    return math.sqrt(gamma * R * T / m)


if __name__ == '__main__':
    # Analysis of air
    T = 273.15 + 20

    v = np.linspace(1, 5000, 5000)
    plt.plot(v, mb_dist(MASS_AIR, v, T))
    plt.title('Maxwell-Boltzmann Distribution of Air at 20C')
    plt.xlabel('Speed (m/s)')
    plt.ylabel('Probablity Density (s/m)')
    plt.show()

    print('v_rms of air is', get_rms_velocity(MASS_AIR, T))
    print('speed of sound in air at', T, 'K is', get_speed_of_sound(GAMMA_AIR, molar_mass_air, T))

    # Analysis of He2
    T = 273.15 + 20

    v = np.linspace(1, 5000, 5000)
    plt.plot(v, mb_dist(MASS_HE2, v, T))
    plt.title('Maxwell-Boltzmann Distribution of He2 at 20C')
    plt.xlabel('Speed (m/s)')
    plt.ylabel('Probablity Density (s/m)')
    plt.show()

    print('v_rms of He2 is', get_rms_velocity(MASS_HE2, T))
    print('speed of sound in He2 at', T, 'K is', get_speed_of_sound(GAMMA_HE2, molar_mass_he2, T))

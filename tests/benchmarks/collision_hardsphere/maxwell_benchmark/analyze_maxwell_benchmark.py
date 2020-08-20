"""
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 Simple script for the analysis of maxwell boltzmann benchmark of a field free ion 
 with hard sphere collision model

 @author Walter Wissdorf
 """

import numpy as np
import matplotlib.pyplot as plt
import sys

K_BOLTZMANN = 1.3806505e-23
AMU_TO_KG = 1.66048e-27


def maxwellBoltzmannDensity(v, temp_K, mass_amu):
    """
    Three dimensional Maxwell-Boltzmann density distribution
    """
    mass_kg = mass_amu * AMU_TO_KG    
    dens = 4* np.pi * (mass_kg / (2*np.pi*K_BOLTZMANN*temp_K) )**(3.0/2.0)* v**2.0 * np.exp( -(mass_kg*v**2.0)/(2*K_BOLTZMANN*temp_K)) 
    return(dens)


if __name__ == '__main__':
    # read file with samples from hs_maxwell_benchmark.cpp and plot them in comparison to MB distribution

    samples_fn = sys.argv[1]
    samples = np.genfromtxt(samples_fn)

    samples_mag = [np.linalg.norm(s) for s in samples]

    velocity_bins = np.linspace(0, 1800, 100)
    mb_dens = np.array([maxwellBoltzmannDensity(v, 298, 28) for v in velocity_bins])

    plt.hist(samples_mag, bins = velocity_bins, density = True)
    plt.plot(velocity_bins,mb_dens)

    plt.show()



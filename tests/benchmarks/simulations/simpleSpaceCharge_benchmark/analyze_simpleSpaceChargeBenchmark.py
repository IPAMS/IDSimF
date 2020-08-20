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
import IDSimF_analysis as ia

K_BOLTZMANN = 1.3806505e-23
AMU_TO_KG = 1.66048e-27


if __name__ == '__main__':
    # read file with samples from hs_maxwell_benchmark.cpp and plot them in comparison to MB distribution

    tra_basename = sys.argv[1]
    result_name = sys.argv[2]
    #dat = ia.read_hdf5_trajectory_file(tra_file_name)
    ia.render_scatter_animation(tra_basename, result_name)


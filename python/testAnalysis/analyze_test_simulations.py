# analyze results of test_integrationQuality.cpp

import numpy as np
import matplotlib.pyplot as plt
import sys
import os.path


def analyze(result_path, result_filename):
    data = np.genfromtxt(os.path.join(result_path, result_filename+'.txt'))
    fig,ax = plt.subplots(1,1, figsize=(10,10))
    ax.scatter(data[:,1], data[:,2], s=2, c=data[:,0], cmap='plasma')
    ax.set_xlabel('position (m)')
    ax.set_ylabel('velocity (m/s)')
    fig.savefig(result_filename+'.png', dpi=200)

if __name__ == '__main__':
    if (len(sys.argv) <=1):
        result_path = "."
    else:
        result_path = sys.argv[1]

    analyze(result_path, "integration_test_verlet_serial")





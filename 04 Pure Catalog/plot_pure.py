import os
import numpy as np
import matplotlib.pyplot as plt

# Pick a few parameters for a sample galaxy.

main_dir = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# The user must set this to the directory and file that needs to be read in.
inputfile = main_dir + '/' + 'pure_CAS_results.csv'

d = np.genfromtxt(inputfile, \
    names = 'dr7objid,dr8objid,C_array,A_array,S_array', \
    dtype = 'int64,int64,float64,float64,float64', \
    delimiter=',')
    

plt.scatter(d['C_array'],d['A_array'])

# Add labels to the plot.
plt.xlabel(r'Concentration (C)', fontsize=16)
plt.ylabel(r'Asymmetry (A)', fontsize=16)

plt.title('CAS Values for Galaxies', fontsize=16)

# Set the range for the plot.
plt.axis([0,5, -4,2])

plt.show()
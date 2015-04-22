import os
import numpy as np
import matplotlib.pyplot as plt

main_dir = os.path.dirname(os.path.abspath(__file__)) + '/' # directory of the script being run
#FITS_dir = '/Users/acrider/FITS/DR10/'
FITS_dir = '/Volumes/My Book for Mac/FITS/DR10/'

# The user must set this to the directory and file that needs to be read in.
#inputfile = main_dir  + 'NGC_4713.csv'
inputfile = main_dir + 'SF_2015330.csv'

d = np.genfromtxt(inputfile, \
    names = 'specObjID,bestobjid,oiii_5007_flux,h_beta_flux,nii_6584_flux,h_alpha_flux,elliptical,spiral,other,u,g,r,i,z,petroRadErr_r', \
    dtype = 'int64,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64', \
    delimiter=',', skip_header=1)
# Get a list of all of the unique DR7 and DR8 galaxy IDs.
oldObjID                   = d['bestobjid']
oldObjID, unique_indices   = np.unique(oldObjID, return_index=True)
nUnique                    = len(oldObjID)

color = np.transpose((d['elliptical'][unique_indices], d['other'][unique_indices], d['spiral'][unique_indices])) # Use all three colors...

color = color.clip(0.0, 1.0)
    
u = d['u'][unique_indices]
g = d['g'][unique_indices]
r = d['r'][unique_indices]
i = d['i'][unique_indices]
z = d['z'][unique_indices]

plt.clf()
plt.title(inputfile)
plt.xlabel('r')
plt.ylabel('g-r')
plt.scatter(r,g-r, c=color, s=30)
plt.ylim(-0.5,2)
plt.xlim(20,15)
plt.show()
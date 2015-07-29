import os
import numpy as np
import matplotlib.pyplot as plt

# Pick a few parameters for a sample galaxy.

main_dir = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# The user must set this to the directory and file that needs to be read in.
inputfile = main_dir + '/' + 'pure_CAS_results_545.csv'

d = np.genfromtxt(inputfile, \
    names = 'dr7objid,dr8objid,C_array,A_array,S_array', \
    dtype = 'int64,int64,float64,float64,float64', \
    delimiter = ',')
    
inputfile2 = main_dir + '/' + 'SF_2015330.csv'

d2 = np.genfromtxt(inputfile2, \
    names = 'specObjID,bestobjid,oiii_5007_flux,h_beta_flux,nii_6584_flux,h_alpha_flux,elliptical,spiral,other,u,g,r,i,z,petroRadErr_r', \
    dtype = 'int64,int64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64,float64', \
    delimiter = ',')

#y = d['C_array']
#y = d['A_array']
y = d['S_array']
x = np.zeros_like(y)

for i in np.arange(len(y)):
    i2 = np.where(d2['bestobjid']==d['dr7objid'][i])
    if len(i2[0]) > 0:
        print d['dr7objid'][i], d2['bestobjid'][i2[0][0]]
        x[i] = np.log10(d2['oiii_5007_flux'][i2[0][0]] / d2['h_beta_flux'][i2[0][0]] )
    else:
        print 'NO MATCH'

plt.clf()
plt.scatter(x,y)

# Add labels to the plot.
plt.xlabel(r'$\rm{log}_{10}$([OIII] / $\rm{H}_{\beta}$)',  fontsize=16)
#plt.ylabel(r'Concentration (C)', fontsize=16)
#plt.ylabel(r'Asymmetry (A)', fontsize=16)
plt.ylabel(r'Clumpiness (S)', fontsize=16)

plt.title('CAS Values for Star-Forming Galaxies', fontsize=16)

# Set the range for the plot.
#plt.axis([-1.5,1.5, 0.0, 4.5]) # for C_array
#plt.axis([-1.5,1.5, -0.52,2.0]) # for A_array
plt.axis([-1.5,1.5, -0.7,0.5]) # for S_array

plt.show()
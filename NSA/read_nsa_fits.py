import pyfits
import numpy as np
import matplotlib.pyplot as plt

# These data were pulled manually from various SDSS websites (see above).
nsa_file = '/Users/acrider/FITS/NSA/nsa_v0_1_2.fits' # from http://nsatlas.org/data

# Open the FITS file and read in the table data.
hdulist = pyfits.open(nsa_file)
tbdata  = hdulist[1].data

# Calculate the x and y values using line rations for the BPT diagram.
x = np.log10(tbdata['n2flux'] / tbdata['haflux'])
y = np.log10(tbdata['o3flux'] / tbdata['hbflux'])

sersic_n      = tbdata['SERSIC_N'] 
asymmetry     = tbdata['ASYMMETRY'][:,2] # g=2 
clumpiness    = tbdata['CLUMPY'][:,2] # g=2 

# Use an increasingly dark greyscale for each dot to make the dense areas stand out.
r = np.linspace(0.5,0.0,len(x)) 
color = np.transpose((r,r,r))

# Make a plot of x vs. y with very thing crosses.
p0 = plt.scatter(x, y, s=2, marker='+', facecolor=color, lw=1)

# Label the plot
plt.xlabel(r'$\log_{10}( NII / H_{\alpha})$', fontsize=16)
plt.ylabel(r'$\log_{10}(OIII / H_{\beta})$', fontsize=16)
plt.title('BPT from NASA Sloan Atlas', fontsize=16)

# Set the range for the plot.
#plt.axis([-1.4, 0.3, -1.5,1.3])
plt.axis([-2, 1, -1.25,1.5])

# Now that everything is set, show the plot!
plt.show()
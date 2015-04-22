class SDSS_power_law:
    def __init__(self,a, b, redshift):
        self.redshift = redshift
        self.a = a
        self.b = b
        
    def value(self,wavelength):
        return np.exp(np.log(self.a) + self.b * np.log(wavelength/6500.))
  
def fetch_filter_2(filt):
    
    main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

    # Downloaded from https://www.sdss3.org/instruments/camera.php
    ff = pyfits.open(main_dir + '/filter_curves.fits')
    tbdata = ff[filt].data
    X = tbdata['wavelength']
    Y = tbdata['respt']          
    return np.transpose(np.array([X,Y]))  
    
def SDSS_spec(x,a,b):
    
    wavelength = np.arange(2000.0,9000.0,10.0)
    s = SDSS_power_law(a, b, 0.0)
    flux = s.value(wavelength)
    
    filter_set = np.array(['u','g','r','i','z'])
    counts = [0,0,0,0,0]
    
    for i in np.arange(0,5):
        ugriz = filter_set[i]
        response = fetch_filter_2(ugriz)
        filter_response = interpolate.interp1d(response[:,0],response[:,1], \
            bounds_error=False,fill_value=0.0)
    
        # See http://en.wikipedia.org/wiki/AB_magnitude
        f_nu = 3.34e-13 * (wavelength * wavelength) * (filter_response(wavelength) * flux) # convert from 1e-17 to Jy
        f_nu = f_nu / 3.631e-6 # convert to Jansky then to nanomaggies
        
        counts[i] = np.sum(filter_response(wavelength) * f_nu) 
         
    return counts                            
                                                                                                      
import os
import pyfits                                                      
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# Approximate the continuum spectrum with a power law.
wavelength = np.arange(4000.0,9000.0,10.0)
s = SDSS_power_law(60.0, -0.85, 0.11)
flux = s.value(wavelength)

plt.figure(1)
plt.clf()
plt.plot(wavelength, flux)
plt.xlabel(r'Wavelength', fontsize=16)
plt.ylabel(r'Flux', fontsize=16)
plt.title('Power-Law Spectrum', fontsize=16)
plt.show()

print 'ugriz for flux of 60 at 6000 Angstroms and alpha = 00.85'
print SDSS_spec([0,1,2,3,4],1.5,-3)
    

    
    

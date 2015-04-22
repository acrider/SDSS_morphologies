class SDSS_power_law:
    def __init__(self,a, b, redshift):
        self.redshift = redshift
        self.a = a
        self.b = b
        
    def value(self,wavelength):
        return np.exp(np.log(self.a) + self.b * np.log(wavelength/6500.))
  
def fetch_SDSS_spectrum(filename):
    
    ff = pyfits.open(filename)
    tbdata = ff[1].data
    X = 10.0**tbdata['loglam']
    Y = tbdata['flux']          
    return np.transpose(np.array([X,Y]))
    
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
s = SDSS_power_law(60.0, -0.85, 0.0)
flux = s.value(wavelength)

# Read in an actual spectrum.
Xref = fetch_SDSS_spectrum(main_dir + '/NGC 4713/spec-0847-52426-0556.fits')

plt.plot(wavelength, flux)
plt.xlabel(r'Wavelength', fontsize=16)
plt.ylabel(r'Flux', fontsize=16)
plt.title('Power-Law Spectrum', fontsize=16)
plt.plot(Xref[:, 0], Xref[:, 1], '-k', lw=0.5)
plt.show()


print SDSS_spec([0,1,2,3,4],60,-0.85)
    
#plt.plot(wavelength, filter_response(wavelength) * flux)
#plt.plot(Xref[:,0],  Xref[:,1] * filter_response(Xref[:,0]))
    
    

def make_stamp(large_image, stamp_X, stamp_Y, stampsize):
    """Given a large image this returns a stamp image."""
    xmin = int(stamp_X-stampsize/2)  
    xmax = xmin+stampsize-1
    ymin = np.max([int(stamp_Y-stampsize/2),0])
    ymax = ymin+stampsize-1
    stamp_image = large_image[xmin:xmax+1,ymin:ymax+1]
    return stamp_image

#---BEGIN---#
import os

import img_scale

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from SDSS_dark_gain import *
from power_law import *
from SDSS_objid_to_values import SDSS_objid_to_values
from get_DR10_FITS_library import *
from extract_galaxy import galaxy_layers, galaxy_XY, create_FITS_stamp

from matplotlib.colors import hsv_to_rgb

main_dir = os.path.dirname(os.path.abspath(__file__)) + '/' # directory of the script being run
#FITS_dir = '/Users/acrider/FITS/DR10/'
FITS_dir = '/Volumes/My Book for Mac/FITS/DR10/'

# The user must set this to the directory and file that needs to be read in.
#inputfile = main_dir  + 'NGC_4713.csv'
inputfile = main_dir + 'dr7todr8_SF2.csv'

d = np.genfromtxt(inputfile, \
    names = 'dr7objid,dr8objid,distance,RA,DEC,petroR50_r,petroR90_r', \
    dtype = 'int64,int64,float64,float64,float64, float64,float64', \
    delimiter=',', skip_header=1)

# Get a list of all of the unique DR7 and DR8 galaxy IDs.
oldObjID                   = d['dr7objid']
oldObjID, unique_indices   = np.unique(oldObjID, return_index=True)
nUnique                    = len(oldObjID)

newObjID = d['dr8objid'][unique_indices]
distance = d['distance'][unique_indices]
RA       = d['RA'][unique_indices]
DEC      = d['DEC'][unique_indices]
petroR50 = d['petroR50_r'][unique_indices]
petroR90 = d['petroR90_r'][unique_indices]
boxsize  = np.zeros_like(petroR90, dtype=int)
boxsize[:] = 100
boxsize[:] = 3.0 * petroR90 / 0.396 # SDSS pixels are 0.396" each # For NGC 4713
nGalaxies = len(newObjID)

# Get information about line to sort based on line strength
inputfile2 = main_dir + '/' + 'BPT_SF_data.csv'
d2 = np.genfromtxt(inputfile2, \
    names = 'specObjID,bestobjid,oiii_5007_flux,h_beta_flux,nii_6584_flux,h_alpha_flux,elliptical,spiral,other', \
    dtype = 'int64,int64,float64,float64,float64,float64,float64,float64,', \
    delimiter=',', skip_header=1)

oiii_5007_flux = np.zeros(nGalaxies, dtype=float)
h_beta_flux    = np.zeros(nGalaxies, dtype=float)
nii_6584_flux  = np.zeros(nGalaxies, dtype=float)
h_alpha_flux   = np.zeros(nGalaxies, dtype=float)
       
for k in xrange(0,nGalaxies):
    index             = np.where(d2['bestobjid']==oldObjID[k])[0]
    if len(index) > 0:
        index = index[0]
        oiii_5007_flux[k] = d2['oiii_5007_flux'][index]
        h_beta_flux[k]    = d2['h_beta_flux'][index]
        nii_6584_flux[k]  = d2['nii_6584_flux'][index]
        h_alpha_flux[k]   = d2['h_alpha_flux'][index]
    else:
        print 'ERROR: Unable to find' + str(oldObjID[k])

isort = np.argsort(h_alpha_flux) # Sort (Low to High)
irevsort = isort[::-1] # Reverse sort (High to Low)       

oldObjID = oldObjID[irevsort]
newObjID = newObjID[irevsort]         
RA       = RA[irevsort]
DEC      = DEC[irevsort]
petroR50 = petroR50[irevsort]
petroR90 = petroR90[irevsort]
boxsize  = boxsize[irevsort]
oiii_5007_flux = oiii_5007_flux[irevsort]
h_beta_flux    = h_beta_flux[irevsort]
nii_6584_flux  = nii_6584_flux[irevsort]
h_alpha_flux   = h_alpha_flux[irevsort]                         
                                                      
filters = ['u','g','r','i','z']

#for k in xrange(0,nGalaxies): # Use this to run ALL of the galaxies.
for k in xrange(20,30): # Use this to run ALL of the galaxies.
    
    print oldObjID[k], newObjID[k]

    skyVersion, rerun, run, camcol, ff, field, object_num = SDSS_objid_to_values(newObjID[k])
    gain, dark_var = SDSS_gain_dark(camcol, ugriz, run)

    plt.figure(1)
    plt.clf()
    n_plot_counter = 1
    
    stamp_stack = np.zeros((boxsize[k], boxsize[k],5))
    
    for ugriz in filters: # Read in all five filters with the first loop.
        
        # ...so you can name the FIT file...    
        frame_FITS = FITS_dir + ugriz + '-' + str(run) + '-' + str(camcol) + '-' + str(field) + '.fits'
    
        # ...and then get it from the SDSS database.    
        get_SDSS_FITS(run, rerun, camcol, field, ugriz, frame_FITS)
        
        # Now get the image data itself out of the FITS file.
        frame_image, sky_image, calib_image, dn_image = galaxy_layers(frame_FITS)
        galaxy_X, galaxy_Y = galaxy_XY(RA[k],DEC[k],frame_FITS)
        
        # Calculate the error in the frame image for use fitting algorithms later.
        dn_image        = frame_image / calib_image + sky_image # counts
        dn_err_image    = np.sqrt(dn_image / gain + dark_var)
        frame_image_err = dn_err_image * calib_image

        #stamp = make_stamp(frame_image / frame_image_err, galaxy_X, galaxy_Y, boxsize[k])
        stamp = make_stamp(frame_image, galaxy_X, galaxy_Y, boxsize[k])
        stamp_stack[:,:,n_plot_counter-1] = stamp

        plt.subplot(2,3,n_plot_counter) # 2 down, 3 across, nth plot
        n_plot_counter = n_plot_counter + 1
        plt.imshow(stamp, cmap='gnuplot2', interpolation='none')
        plt.title('SDSS filter = ' + ugriz)
        plt.show()
    
    plt.subplot(2,3,n_plot_counter) # 2 down, 3 across, nth plot
    
    u = stamp_stack[:,:,0]
    g = stamp_stack[:,:,1]
    r = stamp_stack[:,:,2]
    i = stamp_stack[:,:,3]
    z = stamp_stack[:,:,4]
    
    #----

    #plt.imshow(r-i, interpolation='none', cmap='ocean')
    #plt.colorbar()
    #plt.show()

    #----

    #hue = (i - r)
    #hue = hue - np.min(hue)
    #hue = hue / np.max(hue) * 0.69  # The 0.68 stops the spectrum at blue rather than purple/red.
    
    #value = stamp_stack[:,:,2]
    #value = value - np.min(value)
    #value = value / np.max(value)
    
    #saturation = np.ones_like(hue)
    #saturation = value
    
    #HSV = np.dstack((hue, saturation, value))
    #RGB = hsv_to_rgb(HSV)
    #plt.imshow(RGB, interpolation='none')
    #plt.show()
    
    #----

    b_img = u + g + z
    g_img = r 
    r_img = i 
    
    b_img = b_img - np.min(b_img)
    g_img = g_img - np.min(g_img)
    r_img = r_img - np.min(r_img)

    RGB = np.zeros((r_img.shape[0], r_img.shape[1], 3), dtype=float)

    #i_img = (b_img + g_img + r_img) / 3.0
    #RGB[:,:,0] = r_img * img_scale.linear(i_img) / i_img
    #RGB[:,:,1] = g_img * img_scale.linear(i_img) / i_img
    #RGB[:,:,2] = b_img * img_scale.linear(i_img) / i_img
  
    RGB[:,:,0] = img_scale.linear(r_img) 
    RGB[:,:,1] = img_scale.linear(g_img) 
    RGB[:,:,2] = img_scale.asinh(b_img)
      
    #plt.imshow(r_img, interpolation='none', cmap='ocean')
    plt.title(str(oldObjID[k]))
    plt.imshow(RGB, interpolation='none')
    plt.show()
    
    #plt.figure(2)
    #plt.clf()
    
    #----

    #x = stamp_stack[:,:,2] / stamp_stack[:,:,1]
    #y = stamp_stack[:,:,3] / stamp_stack[:,:,2]
    #plt.scatter(x,y)
    #plt.show()
        
    #for j in xrange(0,boxsize):
    #     for i in xrange(0,boxsize):
    #         plt.scatter([0,1,2,3,4],stamp_stack[i,j,:])
     #        plt.show()
              
        

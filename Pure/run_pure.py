#-------------------------------------- 

# IMPORT LIBRARIES   

import sys
import os.path

import urllib
import urllib2
import csv
import time

import pylab as py
import numpy as np
import matplotlib.pyplot as plt

import math as mth
import subprocess
import pyfits
import matplotlib.cm as cm
import scipy.ndimage as ndimage

from scipy.optimize import curve_fit

from get_DR10_FITS_library import *
from SDSS_objid_to_values import SDSS_objid_to_values
from extract_galaxy import galaxy_layers, galaxy_XY, create_FITS_stamp
from conselice_functions import *

# Pick a few parameters for a sample galaxy.

main_dir       = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# The user must set this to the directory and file that needs to be read in.
inputfile = main_dir + '/' + 'dr7todr8_AGN2.csv'

d = np.genfromtxt(inputfile, \
    names = 'dr7objid,dr8objid,distance,RA,DEC,petroR50_r,petroR90_r', \
    dtype = 'int64,int64,float64,float64,float64, float64,float64', \
    delimiter=',', skip_header=1)

# Get a list of all of the unique DR7 and DR8 galaxy IDs.
#oldObjID       = np.unique(DR7_ObjID)
#newObjID       = np.unique(DR8_ObjID)
newObjID       = d['dr8objid']
nUnique        = len(newObjID)

# Create an empty array to store the Sersic fit parameters.
C_array = np.zeros(nUnique)
A_array = np.zeros(nUnique)

# Set a few other things...
FITS_directory = '/Users/acrider/FITS/DR10/'
ugriz = 'g'
#ugriz = 'r'
#ugriz = 'i' # GUNZIP FIRST!
obj = ''

# MAIN LOOP

for i in xrange(0,nUnique): # Use this to run ALL of the galaxies.
#for i in xrange(0,10): # Use this to test that things work.

    objid = newObjID[i]
    print 'objid = ', objid
    
    # Run retrieve_SDSS_params...
    # run, rerun, camcol, field, obj, rowc, colc = retrieve_SDSS_params(objid)
    skyVersion, rerun, run, camcol, ff, field, object_num = SDSS_objid_to_values(objid)
    print 'skyVersion, rerun, run, camcol, ff, field, object_num = ', skyVersion, rerun, run, camcol, ff, field, object_num 

    # ...so you can name the FIT file...    
    frame_FITS = FITS_directory + ugriz + '-' + str(run) + '-' + str(camcol) + '-' + str(field) + '.fits'
    print frame_FITS
    
    # ...and then get it from the SDSS database.    
    get_SDSS_FITS(run, rerun, camcol, field, ugriz, frame_FITS)

    # Now get the image data itself out of the FITS file.
    frame_image, sky_image, calib_image, dn_image = galaxy_layers(frame_FITS)
    print 'frame_image.shape = ', frame_image.shape
    plt.figure(1) # first window
    plt.imshow(frame_image)
    plt.show()
    
    RA  = d['RA'][i]
    DEC = d['DEC'][i]
    print 'RA,DEC = ',RA, DEC
    
    galaxy_X, galaxy_Y = galaxy_XY(RA,DEC,frame_FITS)
    print  'galaxy_X, galaxy_Y = ', galaxy_X, galaxy_Y

    # Extract a small box around the galaxy.
    boxsize = 100
    xmin = int(galaxy_X-boxsize/2)  
    xmax = xmin+boxsize-1
    ymin = np.max([int(galaxy_Y-boxsize/2),0])
    ymax = ymin+boxsize-1
    galaxy_image = frame_image[xmin:xmax,ymin:ymax]
    
    plt.figure(2) # second window
    plt.imshow(galaxy_image, cmap='gray', interpolation='none')
    plt.show()
    
    #sys.exit("Error message")
    
    C_array[i] = concentration(galaxy_image)
    A_array[i] = asymmetry(galaxy_image)
    print objid,  C_array[i],  A_array[i]
    
    print
    print '======'
    print
    
    

    
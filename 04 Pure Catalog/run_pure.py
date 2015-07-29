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

import img_scale

from scipy.optimize import curve_fit

from get_DR10_FITS_library import *
from SDSS_objid_to_values import SDSS_objid_to_values
from extract_galaxy import galaxy_layers, galaxy_XY, create_FITS_stamp
from conselice_functions import *

# Pick a few parameters for a sample galaxy.

main_dir = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# The user must set this to the directory and file that needs to be read in.
inputfile = main_dir + '/' + 'dr7todr8_BOTH2.csv'

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

nUnique        = len(oldObjID)

# Create an empty array to store the Sersic fit parameters.
C_array = np.zeros(nUnique)
A_array = np.zeros(nUnique)
S_array = np.zeros(nUnique)

# Set a few other things...
FITS_directory = '/Volumes/My Book for Mac/FITS/DR10/'
ugriz = 'g'
#ugriz = 'r'
#ugriz = 'i' # GUNZIP FIRST!
obj = ''
filters = ['u','g','r','i','z']

# MAIN LOOP

#for i in xrange(0,nUnique): # Use this to run ALL of the galaxies.
for i in xrange(0,20): # Use this to test that things work.

    objid = newObjID[i]
    print 'objid = ', objid
    
    # Run retrieve_SDSS_params...
    # run, rerun, camcol, field, obj, rowc, colc = retrieve_SDSS_params(objid)
    skyVersion, rerun, run, camcol, ff, field, object_num = SDSS_objid_to_values(objid)
    print 'skyVersion, rerun, run, camcol, ff, field, object_num = ', skyVersion, rerun, run, camcol, ff, field, object_num 

    plt.figure(1) # second window
    plt.clf()
    n_plot_counter = 1
    #plt.suptitle('SDSS ' + str(objid), fontsize=24)

    for ugriz in filters:
        
        # ...so you can name the FIT file...    
        frame_FITS = FITS_directory + ugriz + '-' + str(run) + '-' + str(camcol) + '-' + str(field) + '.fits'
        print frame_FITS
    
        # ...and then get it from the SDSS database.    
        get_SDSS_FITS(run, rerun, camcol, field, ugriz, frame_FITS)

        # Now get the image data itself out of the FITS file.
        frame_image, sky_image, calib_image, dn_image = galaxy_layers(frame_FITS)
        print 'frame_image.shape = ', frame_image.shape
    
        print 'RA,DEC = ',RA[i], DEC[i]
    
        galaxy_X, galaxy_Y = galaxy_XY(RA[i],DEC[i],frame_FITS)
        print  'galaxy_X, galaxy_Y = ', galaxy_X, galaxy_Y

        # Extract a small box around the galaxy.
        boxsize = 100
        xmin = int(galaxy_X-boxsize/2)  
        xmax = xmin+boxsize-1
        ymin = np.max([int(galaxy_Y-boxsize/2),0])
        ymax = ymin+boxsize-1
        galaxy_image = frame_image[xmin:xmax+1,ymin:ymax+1]
    
        if ugriz == 'u':
            SDSS_u_img = np.copy(galaxy_image)
            print 'u'
        elif ugriz == 'g':
            SDSS_g_img = np.copy(galaxy_image)
        elif ugriz == 'r':
            SDSS_r_img = np.copy(galaxy_image)
            print 'r'
        elif ugriz == 'i':
            SDSS_i_img = np.copy(galaxy_image)
        elif ugriz == 'z':
            SDSS_z_img = np.copy(galaxy_image)
        else:
            print 'Error! Filter does not exist!'

        plt.subplot(2,3,n_plot_counter) # 2 down, 3 across, nth plot
        n_plot_counter = n_plot_counter + 1
        plt.imshow(galaxy_image, cmap='gray', interpolation='none')
        plt.title('SDSS filter = ' + ugriz)
        plt.show()
        
    #b_img = np.copy(SDSS_g_img) + np.copy(SDSS_r_img) + np.copy(SDSS_i_img) + np.copy(SDSS_z_img)
    #g_img = np.copy(SDSS_r_img) + np.copy(SDSS_i_img) + np.copy(SDSS_z_img)
    #r_img = np.copy(SDSS_i_img) + np.copy(SDSS_z_img)
    #b_img = np.clip(b_img, 0.0, np.max(b_img))
    #g_img = np.clip(g_img, 0.0, np.max(g_img))
    #r_img = np.clip(r_img, 0.0, np.max(r_img))
    #i_img = (b_img + g_img + r_img) / 3.0
    #img = np.zeros((r_img.shape[0], r_img.shape[1], 3), dtype=float)
    #img[:,:,0] = r_img * img_scale.linear(i_img) / i_img
    #img[:,:,1] = g_img * img_scale.linear(i_img) / i_img
    #img[:,:,2] = b_img * img_scale.linear(i_img) / i_img
    #plt.subplot(236) # 2 down, 3 across, nth plot
    #py.imshow(img, aspect='equal',interpolation='none')
    #plt.title('Color Composite')
    #plt.show()
    
    
    HII_image = np.copy(SDSS_i_img) / (np.copy(SDSS_r_img) + np.copy(SDSS_z_img))
    HII_image = np.clip(HII_image, 0.0, 5)
    plt.subplot(236)
    plt.imshow(HII_image, aspect='equal',interpolation='none')
    plt.title('HII regions')
    plt.show()
    
    # Create an image where each value is the radius from the center.
    lx, ly = SDSS_r_img.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    radius_image = np.sqrt((X - lx / 2) ** 2 + (Y - ly / 2) ** 2)
    
    plt.figure(3)
    plt.clf()
    plt.xlim([0,20])
    plt.scatter(radius_image, HII_image )
    plt.show()

    #plt.figure(2) # first window
    #plt.clf()
    #plt.subplot(131) # 1 down, 3 across, 1st plot
    #C_array[i] = concentration(SDSS_r_img)
    # print 'Concentration (C) = ', C_array[i]
    #plt.subplot(132) # 1 down, 3 across, 1st plot
    #A_array[i] = asymmetry_min(SDSS_r_img)
    #print 'Asymmetry     (A) = ', A_array[i]
    #plt.subplot(133) # 1 down, 3 across, 1st plot
    #S_array[i] = clumpiness(SDSS_r_img)
    #print 'Clumpiness    (S) = ', S_array[i]  
    #print
    #print '======'
    #print

f = open(main_dir + '/pure_CAS_results.csv', 'w')
for i in np.arange(len(oldObjID)):
    tmp_string = str(oldObjID[i]) + ',' + str(newObjID[i]) + ',' + str(C_array[i]) + ',' + str(A_array[i]) + ',' + str(S_array[i])
    print tmp_string
    f.write(tmp_string + '\n')
f.close()
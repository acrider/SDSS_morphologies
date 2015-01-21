# The data below was retrieved from multiple SDSS sites including:
#  http://dr10.sdss3.org/fields/name?name=1237671264962150496
#  http://skyserver.sdss.org/public/en/tools/quicklook/quicksummary.aspx?id=0x112d177540c70060&spec=0x0d3c8b097a006800
#  http://skyserver.sdss.org/public/en/tools/explore/Summary.aspx?plate=847&mjd=52426&fiber=556

import pyfits
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import interpolate

from SDSS_dark_gain import *

get_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run
get_file = get_dir + '/frame-g-006005-2-0199.fits' # NGC 4713 Spiral 
#get_file = get_dir + '/frame-g-005318-5-0041.fits' # NGC 4526 Elliptical 

# These data were pulled manually from various SDSS websites (see above).
galaxy_RA  = 192.491112046059 # degrees
galaxy_DEC = 5.31141005442379 # degrees
galaxy_PetroRad = 53.92 # arcsec

# Open the FITS file.
FITS_file = pyfits.open(get_file)

# These values allow translation from RA,DEC to X,Y and vice versa.
crpix1 = FITS_file[0].header['CRPIX1'] #   X of reference pixel
crpix2 = FITS_file[0].header['CRPIX2'] #   Y of reference pixel
crval1 = FITS_file[0].header['CRVAL1'] #  RA of reference pixel
crval2 = FITS_file[0].header['CRVAL2'] # DEC of reference pixel
cd1_1  = FITS_file[0].header['CD1_1']  #  RA deg per column pixel
cd1_2  = FITS_file[0].header['CD1_2']  #  RA deg per row pixel
cd2_1  = FITS_file[0].header['CD2_1']  # DEC deg per column pixel
cd2_2  = FITS_file[0].header['CD2_2']  # DEC deg per row pixel

# Determine the gain (electrons/count) and dark variance from the CAMCOL, FILTER, and RUN
camcol     = FITS_file[0].header['CAMCOL']  # camcol
ugriz      = FITS_file[0].header['FILTER']  # ugriz filter
run        = FITS_file[0].header['RUN']     # run
gain, dark_var = SDSS_gain_dark(camcol, ugriz, run)

frame_image = FITS_file[0].data.transpose()

# Create SKY and CALIBRATION images.
# http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html

allsky     = FITS_file[2].data['ALLSKY'].transpose()
allsky     = allsky[:,:,0]
xinterp    = FITS_file[2].data['XINTERP'].transpose()
xinterp    = xinterp[:,0]
yinterp    = FITS_file[2].data['YINTERP'].transpose()
yinterp    = yinterp[:,0]

sky_function = interpolate.interp2d(np.arange(allsky.shape[1]), np.arange(allsky.shape[0]), allsky, kind='linear')
sky_image    = sky_function(yinterp, xinterp) # in counts

calib     = FITS_file[1].data #  nanomaggies per count
calib_image = np.empty_like(frame_image)
for i in np.arange(calib_image.shape[1]):
    calib_image[:,i] = calib
    
# Close the FITS file.
FITS_file.close()

# Calculate the error in the frame image for use fitting algorithms later.
dn_image        = frame_image / calib_image + sky_image # counts

dn_err_image    = np.sqrt(dn_image / gain + dark_var)
frame_image_err = dn_err_image * calib_image

# Find the X,Y values of the galaxy's RA and DEC.
galaxy_X = crpix1 + (galaxy_RA - crval1) / cd1_1 # I've temporarily left out the cd2_1 and cd1_2 components since they're small.
galaxy_Y = crpix2 + (galaxy_DEC - crval2) / cd2_2
galaxy_R = galaxy_PetroRad / 3600.0 / abs(cd1_1)

# Extract a small box around the galaxy.
boxsize = int(galaxy_R * 2.0 * 1.5)
xmin = int(galaxy_X-boxsize/2)  
xmax = xmin+boxsize-1
ymin = int(galaxy_Y-boxsize/2)
ymax = ymin+boxsize-1
galaxy_image = frame_image[xmin:xmax,ymin:ymax] # Why did I have to switch these from image_full[xmin:xmax,ymin:ymax]?

# Make a plot to show the processed, background, calibration, and original images. And also the galaxy.
plt.figure(1) # first window
plt.subplot(141) # 1 down, 4 across, 1st plot
plt.title('Final Image (in nmgy)')
plt.imshow(frame_image, norm=LogNorm())
plt.subplot(142) # 1 down, 4 across,  2nd plot
plt.title('Sky Background (in counts)')
plt.imshow(sky_image, norm=LogNorm())
plt.subplot(143) # 1 down, 4 across,  3rd plot
plt.title('Calibration (in nmgy/count)')
plt.imshow(calib_image, norm=LogNorm())
plt.subplot(144) # 1 down, 4 across,  3rd plot
plt.title('Original Image (in counts)')
plt.imshow(dn_image, norm=LogNorm())


plt.figure(2) # second window
plt.imshow(galaxy_image, cmap='gray', norm=LogNorm())

old_header = FITS_file[0].header
new_header      = old_header
new_header['RA']   = galaxy_RA
new_header['DEC']  = galaxy_DEC
# STILL NEED TO CHANGE FITS HEADER COMMENTS FOR THESE

hdu1 = pyfits.PrimaryHDU(galaxy_image,header=new_header)
hdu2 = pyfits.ImageHDU(sky_image[xmin:xmax,ymin:ymax],name='SKY_IMAGE')
hdu3 = pyfits.ImageHDU(calib_image[xmin:xmax,ymin:ymax],name='CALIB_IMAGE')

hdulist = pyfits.HDUList([hdu1,hdu2,hdu3])
hdulist.writeto(get_dir + '/galaxy_stamp.fits',clobber=True)

plt.show()


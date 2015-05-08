# FROM http://www.astrobetter.com/making-rgb-images-from-fits-files-with-pythonmatplotlib/

import os
import pyfits
import numpy as np
import pylab as py
import img_scale
from conselice_functions import make_mask_image
import matplotlib.pyplot as plt
 
main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run
main_dir  = main_dir + '/Haro_29'

SDSS_u_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_u.fits')
SDSS_g_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_g.fits')
SDSS_r_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_r.fits')
SDSS_i_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_i.fits')
SDSS_z_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_z.fits')

b_img = make_mask_image(SDSS_g_img,90)
g_img = make_mask_image(SDSS_r_img,90)
r_img = make_mask_image(SDSS_i_img,90)

b_img = np.clip(b_img, 0.0, np.max(b_img))
g_img = np.clip(g_img, 0.0, np.max(g_img))
r_img = np.clip(r_img, 0.0, np.max(r_img))

i_img = (b_img + g_img + r_img) / 3.0

img = np.zeros((r_img.shape[0], r_img.shape[1], 3), dtype=float)
img[:,:,0] = r_img * img_scale.linear(i_img) / i_img
img[:,:,1] = g_img * img_scale.linear(i_img) / i_img
img[:,:,2] = b_img * img_scale.linear(i_img) / i_img

plt.figure(1)
py.clf()
py.imshow(img,interpolation='none')
#py.imshow(make_mask_image(SDSS_g_img,90)[300:500,300:500], aspect='equal',interpolation='none',cmap='cubehelix')
py.title('SDSS - RGB band')
py.savefig(main_dir + '/' + 'my_rgb_image.png')

plt.figure(2)
py.clf()
py.imshow(SDSS_g_img,interpolation='none',cmap='cubehelix')
#py.imshow(make_mask_image(SDSS_g_img,90)[300:500,300:500], aspect='equal',interpolation='none',cmap='cubehelix')
py.title('SDSS - g band')

# Create an image where each value is the radius from the center.
lx, ly = SDSS_r_img.shape
X, Y = np.ogrid[0:lx, 0:ly]
radius_image = np.sqrt((X - lx / 2) ** 2 + (Y - ly / 2) ** 2)

plt.figure(3)
plt.clf()
plt.xlabel(r'radius from center (in pixels)', fontsize=16)
plt.ylabel(r'r-band / i-band', fontsize=16)
#plt.title(main_dir, fontsize=16)
plt.scatter(make_mask_image(radius_image,90), make_mask_image(SDSS_i_img,90) / make_mask_image(SDSS_z_img,90), s=0.1)
plt.show()

plt.figure(4)
plt.clf()
r = make_mask_image(SDSS_r_img,90)
i = make_mask_image(SDSS_i_img,90)
plt.title(r'r-band - i-band ($H_\alpha$)', fontsize=16)
plt.imshow((SDSS_g_img + SDSS_r_img) - (SDSS_i_img + SDSS_z_img), aspect='equal',interpolation='none', cmap='bwr')
plt.show()
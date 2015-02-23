# FROM http://www.astrobetter.com/making-rgb-images-from-fits-files-with-pythonmatplotlib/

import os
import pyfits
import numpy as np
import pylab as py
import img_scale
from conselice_functions import make_mask_image
 
main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

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

py.clf()
py.imshow(img, aspect='equal',interpolation='none')
py.title('SDSS - g band')
py.savefig(main_dir + '/' + 'my_rgb_image.png')
# FROM http://www.astrobetter.com/making-rgb-images-from-fits-files-with-pythonmatplotlib/

import os
import pyfits
import numpy as np
import pylab as py
import img_scale
 
main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

SDSS_u_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_u.fits')
SDSS_g_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_g.fits')
SDSS_r_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_r.fits')
SDSS_i_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_i.fits')
SDSS_z_img =  pyfits.getdata(main_dir + '/' + 'galaxy_stamp_z.fits')

b_img = SDSS_g_img 
g_img = SDSS_r_img 
r_img = SDSS_i_img

img = np.zeros((r_img.shape[0], r_img.shape[1], 3), dtype=float)
img[:,:,0] = img_scale.linear(r_img)
img[:,:,1] = img_scale.linear(g_img) 
img[:,:,2] = img_scale.linear(b_img) 

py.clf()
py.imshow(img, aspect='equal',interpolation='none')
py.title('Blue = g, Green = r, Red = i')
py.savefig(main_dir + '/' + 'my_rgb_image.png')
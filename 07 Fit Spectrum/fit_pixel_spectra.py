def make_mask_image(image, pr):
    
    lx, ly = image.shape
    
    mask_image = np.copy(image)
    lx, ly = mask_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]

    mask = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (pr)**2
    mask_image[mask] = 0
    
    return mask_image
    
# FROM http://www.astrobetter.com/making-rgb-images-from-fits-files-with-pythonmatplotlib/

import os
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import img_scale
import scipy.optimize as optimization

import power_law

main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# Read in the DN_IMAGE for each filter.
SDSS_u_img =  pyfits.getdata(main_dir + '/NGC 4713/' + 'galaxy_stamp_u.fits', 'DN_IMAGE')
SDSS_g_img =  pyfits.getdata(main_dir + '/NGC 4713/' + 'galaxy_stamp_g.fits', 'DN_IMAGE')
SDSS_r_img =  pyfits.getdata(main_dir + '/NGC 4713/' + 'galaxy_stamp_r.fits', 'DN_IMAGE')
SDSS_i_img =  pyfits.getdata(main_dir + '/NGC 4713/' + 'galaxy_stamp_i.fits', 'DN_IMAGE')
SDSS_z_img =  pyfits.getdata(main_dir + '/NGC 4713/' + 'galaxy_stamp_z.fits', 'DN_IMAGE')

plt.clf()
plt.subplot(1,5,1)
plt.imshow(make_mask_image(SDSS_u_img,150) )
plt.subplot(1,5,2)
plt.imshow(make_mask_image(SDSS_g_img,150) )
plt.subplot(1,5,3)
plt.imshow(make_mask_image(SDSS_r_img,150) )
plt.subplot(1,5,4)
plt.imshow(make_mask_image(SDSS_i_img,150) )
plt.subplot(1,5,5)
plt.imshow(make_mask_image(SDSS_z_img,150) )

N = SDSS_z_img.shape[0]

for j in np.arange(3):
    for i in np.arange(3):
        x = [0,1,2,3,4]
        observed = [SDSS_u_img[i,j], SDSS_g_img[i,j], SDSS_r_img[i,j], SDSS_i_img[i,j], SDSS_z_img[i,j]]
        a = SDSS_r_img[i,j] / 100.
        b = (SDSS_i_img[i,j] - SDSS_g_img[i,j]) / (7480.0 - 4686.0)
        p0 = [a,b]
        print observed
        print power_law.SDSS_spec(x,a,b)
        popt,pcov = optimization.curve_fit(power_law.SDSS_spec, x, observed)
        
#print optimization.curve_fit(func, xdata, ydata, x0, sigma)



plt.show()
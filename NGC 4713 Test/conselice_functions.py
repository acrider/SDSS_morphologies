def aperture_crider(image, radius):
    mask_image = np.copy(image)
    lx, ly = mask_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    mask = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > radius**2
    mask_image[mask] = 0
    n_pixels = float(lx*ly - mask.sum())
    return float(mask_image.sum()), n_pixels
    
def petrosian_radius(image):
    # http://spiff.rit.edu/classes/phys443/lectures/gal_1/petro/petro.html
    lx, ly = image.shape
    n = min([lx/2, ly/2])

    r = np.arange(n)
    counts = np.arange(n, dtype=np.float)
    pixels = np.arange(n)

    for i in r:
        counts[i], pixels[i] = aperture_crider(image, i)
            
    local_intensity = (counts[1:n-1] - counts[0:n-2])/(pixels[1:n-1] - pixels[0:n-2])
    eta = local_intensity / (counts[1:n-1] / pixels[1:n-1])
    
    # Use ETA cutoff of 0.2 as SDSS does.
    pr = r[np.argmax(eta<=0.2)]
    return pr, counts
    
def compactness(image):
    pr, counts = petrosian_radius(image)
    r = np.arange(len(counts))
    
    c100 = counts[min([1.5 * pr, len(counts)-1])]

    r80 = r[np.argmin(counts <= 0.8 * c100)]
    r20 = r[np.argmin(counts <= 0.2 * c100)]
    
    print 'Petrosian radius = ', petrosian_radius, 'pixels'
    print 'r80 radius = ', r80, 'pixels'
    print 'r20 radius = ', r20, 'pixels'

    #plt.plot(r, counts, 'b-')
    #plt.xlabel('radius (in pixels)')
    #plt.ylabel('curve of growth')
    #plt.plot([r20,r20],[min(counts),max(counts)], 'r--')
    #plt.plot([r80,r80],[min(counts),max(counts)], 'r--')
    #plt.show()
    
    # Compactness = 5.0 * np.log10(r80/r20) from Bershady< Jangren, & Conselice 2000
    print c100, r80, r20
    return 5 * np.log10(float(r80)/float(r20))
    
def assymetry(image):
    pr, counts = petrosian_radius(image)
    mask_image = np.copy(image)
    lx, ly = mask_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    mask = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > pr**2
    mask_image[mask] = 0
    mask_image_180 = np.rot90(mask_image,2) # Rotate 180 degrees.
    
    #plt.imshow(np.abs(mask_image - mask_image_180))
    #plt.show()

    assymetry = 0.5 * np.sum(np.abs(mask_image - mask_image_180)) / np.sum(mask_image)
    return assymetry

#---

import os    
import numpy as np
import matplotlib.pyplot as plt
import pyfits

main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run
sample_file = main_dir + '/galaxy_stamp.fits'

image = pyfits.getdata(sample_file,0)

c = compactness(image)
print 'Compactness = ', c
print 'Assymetry   = ', assymetry(image)

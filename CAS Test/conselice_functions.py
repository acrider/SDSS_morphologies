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
    
def concentration(image):
    pr, counts = petrosian_radius(image)
    r = np.arange(len(counts))
    
    c100 = counts[min([1.5 * pr, len(counts)-1])]

    r80 = r[np.argmin(counts <= 0.8 * c100)]
    r20 = r[np.argmin(counts <= 0.2 * c100)]
    
    #print 'Petrosian radius = ', pr, 'pixels'
    #print 'r80 radius = ', r80, 'pixels'
    #print 'r20 radius = ', r20, 'pixels'

    #plt.plot(r, counts, 'b-')
    #plt.xlabel('radius (in pixels)')
    #plt.ylabel('curve of growth')
    #plt.plot([r20,r20],[min(counts),max(counts)], 'r--')
    #plt.plot([r80,r80],[min(counts),max(counts)], 'r--')
    #plt.show()
    
    lx, ly = image.shape
    
    # Show plot of data used to calculate assymetry
    plt.imshow(image, cmap='gray_r', norm=LogNorm(), interpolation='none')
    plt.axhline(lx/2., color='g')
    plt.axvline(lx/2., color='g')
    circle1 = plt.Circle((lx/2.,ly/2.), radius=r20, color='g', fill=False)
    circle2 = plt.Circle((lx/2.,ly/2.), radius=r80, color='g', fill=False)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)

    plt.show()
    
    # Compactness = 5.0 * np.log10(r80/r20) from Bershady< Jangren, & Conselice 2000
    # print c100, r80, r20
    return 5 * np.log10(float(r80)/float(r20))
    
def assymetry(image):
    pr, counts = petrosian_radius(image)
    mask_image = np.copy(image)
    lx, ly = mask_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    mask = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (1.5*pr)**2
    mask_image[mask] = 0
    mask_image_180 = np.rot90(mask_image,2) # Rotate 180 degrees.
    
    # Show plot of data used to calculate assymetry
    plt.imshow(mask_image - mask_image_180, cmap='gray_r', interpolation='none')
    plt.axhline(lx/2., color='g')
    plt.axvline(lx/2., color='g')
    circle1 = plt.Circle((lx/2.,ly/2.), radius=1.50*pr, color='g', fill=False)
    plt.gca().add_patch(circle1)
    plt.show()

    assymetry = 0.5 * np.sum(np.abs(mask_image - mask_image_180)) / np.sum(mask_image)
    return assymetry

def clumpiness(image):
    
    pr, counts = petrosian_radius(image)
    
    copy_image = np.copy(image)
    lx, ly = copy_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]

    #Convolve image with a 2D boxcar the size of the 1/4 Petrosian radius.
    boxcar_kernel  = Box2DKernel(0.25 * pr, mode='center')
    smoothed_image = convolve(copy_image, boxcar_kernel)

    # See Equation 11 for clumpiness from Lotz, Primack, and Madua (2004)
    clump_image = copy_image - smoothed_image
    
    #background_clump_image = np.copy(clump_image)
   
    mask1 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (1.50*pr)**2
    mask2 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 < (0.25*pr)**2
    
    clump_image[mask1] = 0
    clump_image[mask2] = 0
    copy_image[mask1] = 0
    copy_image[mask2] = 0
    
    # Show plot of data used to calculate assymetry
    plt.imshow(clump_image, cmap='gray_r',interpolation='none')
    plt.axhline(lx/2., color='g')
    plt.axvline(lx/2., color='g')
    circle1 = plt.Circle((lx/2.,ly/2.), radius=1.50*pr, color='g', fill=False)
    circle2 = plt.Circle((lx/2.,ly/2.), radius=0.25*pr, color='g', fill=False)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)


    plt.show()
    
    #galaxy_clumpiness =  clump_image[clump_image>0].sum() / image[image>0].sum()
    galaxy_clumpiness = np.sum(abs(clump_image)) / np.sum(abs(copy_image))
 
    # Use background annulus of 2 to 2.5 Petrosian radii
    #background_image       = np.copy(image)
    #mask3 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (2.5*pr)**2
    #mask4 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 < (2.0*pr)**2
    #background_clump_image[mask3] = 0
    #background_clump_image[mask4] = 0
    #background_image[mask3] = 0
    #background_image[mask4] = 0
    #background_clumpiness = np.sum(background_clump_image)/np.sum(background_image)
    #print background_clumpiness
    
    background_clumpiness = 0
    
    return galaxy_clumpiness - background_clumpiness
    
#---
import os    
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from astropy.convolution import convolve, Tophat2DKernel, Box2DKernel
from matplotlib.colors import LogNorm

main_dir  = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run
sample_file = main_dir + '/galaxy_stamp.fits'

image = pyfits.getdata(sample_file,0)

plt.figure(3) # first window
plt.subplot(131) # 1 down, 3 across, 1st plot
print 'Concentration (C) = ', concentration(image)
plt.subplot(132) # 1 down, 3 across, 1st plot
print 'Assymetry     (A) = ', assymetry(image)
plt.subplot(133) # 1 down, 3 across, 1st plot
print 'Clumpiness    (S) = ', clumpiness(image)

def rebin( a, newshape ):
    '''Rebin an array to a new shape.
    '''
    assert len(a.shape) == len(newshape)
    
    slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')   #choose the biggest smaller integer index
    return a[tuple(indices)]
    
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
    minimum_pr = 5 # set minimum to 5 pixels
    pr = max([r[np.argmax(eta<=0.2)], minimum_pr]) 
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
    
    mask_image = np.copy(image)
    lx, ly = mask_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]

    mask = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (1.5*pr)**2
    mask_image[mask] = 0
    
    # Show plot of data used to calculate assymetry
    plt.imshow(mask_image, cmap='gray_r', norm=LogNorm(vmin=0.01), interpolation='none')
    plt.axhline(lx/2., color='g')
    plt.axvline(lx/2., color='g')
    circle1 = plt.Circle((lx/2.,ly/2.), radius=r20, color='g', fill=False)
    circle2 = plt.Circle((lx/2.,ly/2.), radius=r80, color='g', fill=False)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    
    # Compactness = 5.0 * np.log10(r80/r20) from Bershady< Jangren, & Conselice 2000
    # print c100, r80, r20
    return 5 * np.log10(float(r80)/float(r20))
    
def asymmetry_min(image):
    
    Nx, Ny = image.shape
    
    xoff = 3
    yoff = 3
    
    A_matrix = np.zeros((xoff,yoff))
    B_matrix = np.zeros((xoff,yoff))
    
    for dy in np.arange(yoff):
        for dx in np.arange(xoff):
            A_matrix[dx,dy], B_matrix[dx,dy] = asymmetry(image[dx:Nx-xoff+dx+1,dy:Ny-yoff+dy+1])
            #A_matrix[dx,dy] = asymmetry(image[dx:Nx+1,dy:Ny+1])

            #print dx, dy, A_matrix[dx,dy]
    
    asymmetry_corrected = np.min(A_matrix) - np.min(B_matrix) # See Conselice, Bershady, and Jangren (2000)
             
    return asymmetry_corrected 
                
def asymmetry(image):
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
    circle2 = plt.Circle((lx/2.,ly/2.), radius=2.00*pr, color='r', fill=False)
    circle3 = plt.Circle((lx/2.,ly/2.), radius=2.50*pr, color='r', fill=False)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    plt.gca().add_patch(circle3)

    plt.show()
    
    # Use background annulus of 2 to 2.5 Petrosian radii (same area as r=1.5)
    background_image       = np.copy(image)
    mask3 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (2.5*pr)**2
    mask4 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 < (2.0*pr)**2
    background_image[mask3] = 0
    background_image[mask4] = 0
    background_image_180 = np.rot90(background_image,2)
    
    npixels_galaxy     = image.size - np.sum(mask)
    npixels_background = image.size - np.sum(mask3) - np.sum(mask4)
    
    bg = np.sum(np.abs(background_image - background_image_180))
    bg = bg * npixels_galaxy / npixels_background
    bg = bg / np.sum(mask_image)
        
    # Old version had 0.5 in it. Why?!?    
    # asymmetry = 0.5 * np.sum(np.abs(mask_image - mask_image_180)) / np.sum(mask_image)
    asymmetry = (np.sum(np.abs(mask_image - mask_image_180))) / np.sum(mask_image)

    return asymmetry, bg
    
def clumpiness(image, **kwargs):
    
    pr, counts = petrosian_radius(image)
        
    copy_image = np.copy(image)
    lx, ly = copy_image.shape
    X, Y = np.ogrid[0:lx, 0:ly]

    #Convolve image with a 2D boxcar the size of the 1/4 Petrosian radius.
    smoothing_kernel  = Box2DKernel(0.25 * pr)
    #smoothing_kernel  = Gaussian2DKernel(0.25 * pr)

    smoothed_image = convolve(copy_image, smoothing_kernel)

    # See Equation 11 for clumpiness from Lotz, Primack, and Madua (2004)
    clump_image = copy_image - smoothed_image
    #clump_image = clump_image * (np.sign(clump_image)+1)/2# From p. 11 of Conselice (?)
    
    # Show plot of data used to calculate assymetry
    plt.imshow(clump_image, interpolation='none')
    plt.axhline(lx/2., color='g')
    plt.axvline(lx/2., color='g')
    circle1 = plt.Circle((lx/2.,ly/2.), radius=1.50*pr, color='g', fill=False)
    circle2 = plt.Circle((lx/2.,ly/2.), radius=0.25*pr, color='g', fill=False)
    circle3 = plt.Circle((lx/2.,ly/2.), radius=2.00*pr, color='r', fill=False)
    circle4 = plt.Circle((lx/2.,ly/2.), radius=2.50*pr, color='r', fill=False)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    plt.gca().add_patch(circle3)
    plt.gca().add_patch(circle4)

    plt.show()

    mask1 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (1.50*pr)**2
    mask2 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 < (0.25*pr)**2
    
    clump_image[mask1] = 0
    clump_image[mask2] = 0
    copy_image[mask1] = 0
    copy_image[mask2] = 0
    
    # print np.sum(abs(clump_image)), np.sum(abs(copy_image)) 

    galaxy_clumpiness = np.sum(abs(clump_image)) / np.sum(abs(copy_image)) # LPM04
    #galaxy_clumpiness = np.sum((clump_image-copy_image)/copy_image) # Conselice 2003
    #galaxy_clumpiness = np.sum(np.clip(clump_image,0.0,np.max(clump_image))) / np.sum(np.clip(copy_image,0.0,np.max(copy_image))) # Conselice 2003

    # print 'Galaxy Clumpiness   = ', galaxy_clumpiness
    
    # Use background annulus of 2.0 to 2.5 Petrosian radii
    bg_copy_image       = np.copy(image)
        
    if 'bg' in kwargs:
        bg_smoothed_image = kwargs['bg']
        bg_smoothed_image = bg_smoothed_image * np.sum(bg_copy_image) / np.sum(bg_smoothed_image)
    else:
        bg_smoothed_image = np.copy(smoothed_image)
        
    bg_clump_image      = bg_copy_image - bg_smoothed_image # OPTION A

    plt.figure(3)
    plt.subplot(131)
    plt.imshow(image, interpolation='none')
    plt.subplot(132)
    plt.imshow(bg_smoothed_image, interpolation='none')
    plt.subplot(133)
    plt.imshow(bg_clump_image, interpolation='none')
    plt.show()
            
    mask3 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > (2.5*pr)**2
    mask4 = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 < (2.0*pr)**2
    bg_clump_image[mask3] = 0
    bg_clump_image[mask4] = 0
    bg_copy_image[mask3] = 0
    bg_copy_image[mask4] = 0
    
    nz12 = np.ones((lx,ly))
    nz12[mask1] = 0
    nz12[mask2] = 0
    nz34 = np.ones((lx,ly))
    nz34[mask3] = 0
    nz34[mask4] = 0
    
    # print 'Number of pixels in Mask1 + Mask 2', np.sum(nz12)
    # print 'Number of pixels in Mask3 + Mask 4', np.sum(nz34)
    
    # print np.sum(abs(bg_clump_image)), np.sum(abs(copy_image)) 
    
    bg_clumpiness = np.sum(abs(bg_clump_image))/np.sum(abs(copy_image)) * (np.sum(nz12)/np.sum(nz34))  # Denominator is the IMAGE, not the background.
    #bg_clumpiness = 0
    
    # print 'Background Clumpiness = ', bg_clumpiness
    
    print 'EXCEL',lx, galaxy_clumpiness, bg_clumpiness, galaxy_clumpiness - bg_clumpiness
        
    return galaxy_clumpiness - bg_clumpiness
    
#---

import pdb
import os    
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import sys
from astropy.convolution import convolve, Tophat2DKernel, Box2DKernel, Gaussian2DKernel
from matplotlib.colors import LogNorm, Normalize

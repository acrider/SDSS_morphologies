import csv
import pyfits
import numpy as np
import matplotlib.pyplot as plt
from SDSS_objid_to_values import SDSS_values_to_objid

# Set this to your own FITS data path!
FITS_dir  = '/Users/acrider/FITS/NSA/'

# These data were pulled manually from various SDSS websites (see above).
sdss_file = FITS_dir + 'sdss_atlas.fits' # from http://nsatlas.org/data

# Open the FITS file and read in the table data.
hdulist = pyfits.open(sdss_file)
tbdata  = hdulist[1].data

# http://www.sdss.org/dr12/help/glossary/#S
SKYVERSION = 2   # for DR8
RERUN      = 301 # for DR8
FF         = 0   # just a guess!!!

N = len(tbdata['RUN'])

DR8_OBJIDS = np.empty([N],dtype=int)

for i in range(N):
    DR8_OBJIDS[i] = SDSS_values_to_objid(SKYVERSION, RERUN, \
    tbdata['RUN'][i], \
    tbdata['CAMCOL'][i], \
    FF, \
    tbdata['FIELD'][i], \
    tbdata['ID'][i])
                                   
# from https://pythonhosted.org/pyfits/users_guide/users_tutorial.html
col1 = pyfits.Column(name='DR8_OBJID',format='K',array=DR8_OBJIDS)
cols = pyfits.ColDefs([col1])

tbhdu = pyfits.BinTableHDU.from_columns([col1])
print 
tbhdu.writeto( FITS_dir + 'sdss_atlas_objids.fits', clobber=True) 

print ' COMPLETE!'
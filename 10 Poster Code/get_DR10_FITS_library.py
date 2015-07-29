def get_SDSS_FITS(run, rerun, camcol, field, ugriz, put_file):
    
    # This function downloads an SDSS image FITS file given four parameters and the color and places it in put_file.
    
    image = urllib.URLopener()

    #get_file = 'http://das.sdss.org/www/cgi-bin/drC?RUN=' + run + '&RERUN=' + rerun + '&CAMCOL=' + camcol + '&FIELD=' + field + '&FILTER=' + ugriz
    
    #http://dr10.sdss3.org/sas/dr10/boss/photoObj/frames/301/94/3/frame-g-000094-3-0159.fits.bz2
    get_file = 'http://data.sdss3.org/sas/dr10/boss/photoObj/frames/301/' + str(run) + '/' + str(camcol) + '/frame-' + ugriz  + '-' + str(int(run)).zfill(6) + '-' + str(camcol) + '-' + str(int(field)).zfill(4) + '.fits.bz2'

    if not(os.path.isfile(put_file)):
        print 'Retrieving ' + get_file
        
        if get_file.endswith('.bz2'):
            put_file = put_file + '.bz2'
        
        print 'Putting ' + put_file
        image.retrieve(get_file, put_file)
    else:
        print 'File exists: ' + put_file
        
    # If you ended up downloading a BZ2 file, then unzip it.    
    if put_file.endswith('.bz2'):
        print 'Unzipping ' + put_file 
        p = subprocess.Popen(['bunzip2 -f "' + put_file + '"'],
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line,
        retval = p.wait() 
        
#--------------------------------------        

def retrieve_SDSS_params(objid):
    
    # This function returns seven parameters (run, rerun, camcol, field, obj, rowc, colc)
    # when given a DR7 objid number (e.g. from Allen et al. 2014).

    temp = 'http://cas.sdss.org/dr7/en/tools/explore/OETOC.asp?id='
    #temp = 'http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id='
    
    response = urllib2.urlopen(temp+str(objid))
    html     = response.read()
    start    = html.find('onLoad="loadSummary(')
    newId    = html[start+21:start+39] # This should look like 0x08290efc61020051
    newSpec  = html[start+49:start+67] # This should look like 0x0000000000000000
    newURL   = 'http://cas.sdss.org/dr7/en/tools/explore/summary.asp?id=' + newId + '&spec=' + newSpec

    response = urllib2.urlopen(newURL)
    html     = response.read()
    start    = html.find('colc')
    
    num_end = start + 1
    
    # Get the RUN
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    run =  html[num_start:num_end]
    
    # Get the RERUN
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    rerun =  html[num_start:num_end]
    
    # Get the CAMCOL
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    camcol =  html[num_start:num_end]
    
    # Get the FIELD
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    field =  html[num_start:num_end]
    
    # Get the OBJ
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    obj =  html[num_start:num_end]
    
    # Get the ROWC
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    rowc =  html[num_start:num_end]
    rowc = float(rowc)
  
    # Get the COLC
    class_start = html.find('class=',num_end)
    num_start   = html.find('>',class_start) + 1
    num_end     = html.find('<',num_start)
    colc =  html[num_start:num_end]
    colc = float(colc)  
        
    return run, rerun, camcol, field, obj, rowc, colc
    
#--------------------------------------    

def get_SDSS_img(put_file, rowc, colc):
    
    # Sample Code from - https://pythonhosted.org/pyfits/users_guide/users_tutorial.html
    FITS_file= pyfits.open(put_file)
    
    # These are things stored in the FITS file. r_img is the most important one.
    #r_img    = pyfits.getdata(put_file) # OLD VERSION
    r_img    = FITS_file[0].data # FAST VERSION

    #gain     = pyfits.getval(put_file,'GAIN') #GAIN = 4.59999990463257 / Gain averaged over all amplifiers (e/DN)  
    #dark_var = pyfits.getval(put_file,'DARK_VAR') # DARK_VAR = 1. / combined variance from dark cur. and read noise

# I've temporarily hard-wired the settings based on data from the following web page instead of reading them in from :
# http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    #camcol = FITS_file[0].header['CAMCOL'] # DARK_VAR = 1. / combined variance from dark cur. and read noise
    camcol = 1 # NOT TRUE!!!!

    if camcol == 1:
        gain     = 3.845
        dark_var = 15.6025
    elif camcol == 2:
        gain     = 3.855
        dark_var = 1.44
    elif camcol == 3:
        gain     = 3.845
        dark_var = 1.3225
    elif camcol == 4:
        gain     = 3.995
        dark_var = 1.96
    elif camcol == 5:
        gain     = 4.05
        dark_var = 1.1025
    elif camcol == 6:
        gain     = 4.035
        dark_var = 1.8225
    else:
        print 'ERROR: CAMCOL is out of range.'

    #sky_image = pyfits.getdata(put_file,2) # NOT SURE ABOUT THIS
    #sky_image = FITS_file[2].data

    #print sky_image
    
    #sys.exit('STOP: What is in the FITS_file[0].header?')     
    #sky      = FITS_file[0].header['SKY'] # SKY = 132.366674648241 / sky level (DN/pixel) 
    #sky_err  = FITS_file[0].header['SKYERR'] # SKYERR  = 3.76900984380701E-12 / The error of average sky value in the frame    
    sky = 0
    sky_err = 0
    
    #http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html

    FITS_file.close()

    # Extract a 100x100 piece centered on the galaxy.
    boxsize = 100   
    xmin = int(rowc-boxsize/2)  
    xmax = xmin+boxsize-1
    ymin = int(colc-boxsize/2)
    ymax = ymin+boxsize-1
    r_img = r_img[xmin:xmax,ymin:ymax]
    r_img = r_img - 1000 # This is the offset of 1000 built-in to SDSS images.
   
    return r_img, gain, sky, dark_var, sky_err
    
import urllib
import urllib2
import os.path
import subprocess

print 'RUNNING GET_DR10_FITS_library.py'
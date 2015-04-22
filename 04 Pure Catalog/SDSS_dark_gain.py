# +++
# 
# FUNCTION: SDSS_gain_dark
#
# This function return the SDSS gain and dark variance for a given camcol, filter, and run. 
# The values are pulled from the following website:
#
# http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
#
# A. Crider (2015-01-15)
#
# --- 

def SDSS_gain_dark(camcol, ugriz, run):

    if camcol == 1:
        if ugriz == 'u':
            gain = 1.62
            dark = 9.61
        elif ugriz == 'g':
            gain = 3.32
            dark = 15.6025
        elif ugriz == 'r':
            gain = 4.71
            dark = 1.8225
        elif ugriz == 'i':
            gain = 5.165
            dark = 7.84
        elif ugriz == 'z':
            gain = 4.745
            dark = 0.81
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 2:
        if ugriz == 'u':
            if run < 1100:
                gain = 1.595
            elif run > 1100:
                gain = 1.825
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'
            dark = 12.6025
        elif ugriz == 'g':
                gain = 3.855
                dark = 1.44
        elif ugriz == 'r':
                gain = 4.6
                dark = 1.00
        elif ugriz == 'i':
            gain = 6.565
            if run < 1500:
                dark = 5.76
            elif run > 1500:
                dark = 6.25
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'  
        elif ugriz == 'z':
            gain = 5.155
            dark = 1.0
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'
    elif camcol == 3:
        if ugriz == 'u':
            gain = 1.59
            dark = 8.7025
        elif ugriz == 'g':
            gain = 3.845
            dark = 1.3225
        elif ugriz == 'r':
            gain =  4.72
            dark = 1.3225
        elif ugriz == 'i':
            gain = 4.86
            dark = 4.6225
        elif ugriz == 'z':
            gain = 4.885
            dark = 1.0
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'  
    elif camcol == 4:
        if ugriz == 'u':
            gain = 1.6
            dark = 12.6025
        elif ugriz == 'g':
            gain = 3.995
            dark = 1.96
        elif ugriz == 'r':
            gain =  4.76
            dark = 1.3225
        elif ugriz == 'i':
            gain = 4.885
            if run < 1500:
                dark = 6.25
            elif run > 1500:
                dark = 7.5625
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'  
        elif ugriz == 'z':
            gain = 4.775
            if run < 1500:
                dark = 9.61
            elif run > 1500:
                dark = 12.6025
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'  
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'  
    elif camcol == 5:
        if ugriz == 'u':
            gain = 1.47
            dark = 9.3025
        elif ugriz == 'g':
            gain = 4.05
            dark = 1.1025
        elif ugriz == 'r':
            gain = 4.725
            dark = 0.81
        elif ugriz == 'i':
            gain = 4.64
            dark = 7.84
        elif ugriz == 'z':
            gain = 3.48
            if run < 1500:
                dark = 1.8225
            elif run > 1500:
                dark = 2.1025
            else:
                print 'ERROR in SDSS_dark_gain: RUN not set!'  
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'  
    elif camcol == 6:
        if ugriz == 'u':
            gain = 2.17
            dark = 7.0225
        elif ugriz == 'g':
            gain = 4.035
            dark = 1.8225
        elif ugriz == 'r':
            gain = 4.895
            dark = 0.9025
        elif ugriz == 'i':
            gain = 4.76
            dark = 5.0625
        elif ugriz == 'z':
            gain = 4.69
            dark = 1.21
        else:
            print 'ERROR in SDSS_dark_gain: UGRIZ not set!'  
    else:
        print 'ERROR in SDSS_dark_gain: CAMCOL is not set!'
    
    return gain, dark

#-----

import numpy as np

run = 0
#run = 2000

camcol_values = [1,2,3,4,5,6]
ugriz_values  = ['u', 'g', 'r', 'i', 'z']

print
print 'TEST'
print

for camcol in camcol_values:
    line = ''
    for ugriz in ugriz_values:
        gain, dark = SDSS_gain_dark(camcol, ugriz, run)
        line = line + ugriz + str(camcol) + ' '
    line = line + '\n'
    print line
    
print
print 'GAINS'
print

for camcol in camcol_values:
    line = ''
    for ugriz in ugriz_values:
        gain, dark = SDSS_gain_dark(camcol, ugriz, run)
        line = line + str(gain) + ' '
    line = line + '\n'
    print line
    
print  
print 'DARK VARIANCES'
print
   
for camcol in camcol_values:
    line = ''
    for ugriz in ugriz_values:
        gain, dark = SDSS_gain_dark(camcol, ugriz, run)
        line = line + str(dark) + ' '
    line = line + '\n'
    print line
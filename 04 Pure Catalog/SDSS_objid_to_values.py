def SDSS_objid_to_values(objid):
    
    # Determined from http://skyserver.sdss.org/dr7/en/help/docs/algorithm.asp?key=objID
    
    bin_objid = bin(objid)
    bin_objid = bin_objid[2:len(bin_objid)]
    bin_objid = bin_objid.zfill(64)

    empty      = int( '0b' + bin_objid[0],     base=0)
    skyVersion = int( '0b' + bin_objid[1:4+1],   base=0)
    rerun      = int( '0b' + bin_objid[5:15+1],  base=0)
    run        = int( '0b' + bin_objid[16:31+1], base=0)
    camcol     = int( '0b' + bin_objid[32:34+1], base=0)
    ff         = int( '0b' + bin_objid[35+1],    base=0)
    field      = int( '0b' + bin_objid[36:47+1], base=0)
    object_num = int( '0b' + bin_objid[48:63+1], base=0)

    return skyVersion, rerun, run, camcol, ff, field, object_num

def SDSS_values_to_objid(skyVersion, rerun, run, camcol, ff, field, object_num):
    
    # Determined from http://skyserver.sdss.org/dr7/en/help/docs/algorithm.asp?key=objID

    empty_bin      = '0'
    skyVersion_bin = '{0:04b}'.format(skyVersion)
    rerun_bin      = '{0:011b}'.format(rerun)
    run_bin        = '{0:016b}'.format(run)
    camcol_bin     = '{0:03b}'.format(camcol)
    ff_bin         = '{0:01b}'.format(ff)
    field_bin      = '{0:012b}'.format(field)
    object_num_bin = '{0:016b}'.format(object_num)

    objid_bin = empty_bin + skyVersion_bin + rerun_bin + run_bin + camcol_bin + ff_bin + field_bin + object_num_bin

    return int(objid_bin,2)
    
#============   
    
from bitarray import bitarray

DR7_objid  = 587722983372161650

skyVersion, rerun, run, camcol, ff, field, object_num  = SDSS_objid_to_values(DR7_objid)

print 'Values extracted from DR7 Objid include:'
print 'RERUN', 'RUN', 'CAMCOL', 'FIELD', 'OBJ'
print rerun, run, camcol, field, object_num
print

DR12_objid = 1237651191899619530

skyVersion, rerun, run, camcol, ff, field, object_num  = SDSS_objid_to_values(DR12_objid)

print 'Values extracted from DR12 Objid include:'
print 'RERUN', 'RUN', 'CAMCOL', 'FIELD', 'OBJ'
print rerun, run, camcol, field, object_num
print

print 'Reconstructed SDSS objid'
print SDSS_values_to_objid(skyVersion, rerun, run, camcol, ff, field, object_num)
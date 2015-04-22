# Import various Python libraries that you'll need for plotting your Google Docs data.
import os
import csv
import matplotlib.pyplot as plt
import numpy as np

#http://www.peterbe.com/plog/uniqifiers-benchmark
def f2(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked
   
main_dir = os.path.dirname(os.path.abspath(__file__)) # directory of the script being run

# The user must set this to the directory and file that needs to be read in.
inputfile = main_dir + '/' + 'BPT_Both_data.csv'

# Read in the CSV file from your own hard drive.
data = [];
with open(inputfile, 'rb') as f:
    csvReader = csv.reader(f, delimiter=',',skipinitialspace=True)
    headers = csvReader.next()
    for row in csvReader:
        data.append(row);
    data = np.asarray(data)

# Extract x and y from the first two columns of the array.
specObjID       = data[:,0]
objid           = data[:,1]
oiii_5007_flux  = np.double(data[:,2])
h_beta_flux     = np.double(data[:,3])
nii_6584_flux   = np.double(data[:,4])
h_alpha_flux    = np.double(data[:,5])


elliptical      = np.double(data[:,6])
spiral          = np.double(data[:,7])
other           = np.double(data[:,8])

x = np.log10(nii_6584_flux / h_alpha_flux)
y = np.log10(oiii_5007_flux / h_beta_flux)

inputfile = main_dir + '/' + 'pure_CAS_results.csv'

d = np.genfromtxt(inputfile, \
    names = 'dr7objid,dr8objid,C_array,A_array,S_array', \
    dtype = 'int64,int64,float64,float64,float64', \
    delimiter=',', skip_header=0)

itranslate = np.zeros(len(objid),dtype=int)

for i in np.arange(len(objid)):
    print i
    tmp = np.where(d['dr7objid'] == int(objid[i]))[0]
    if len(tmp) > 0:
        itranslate[i] = tmp[0]

translated_C = d['C_array'][itranslate]
translated_A = d['A_array'][itranslate]
translated_S = d['S_array'][itranslate]

clipped_C = np.clip((translated_C-1.0)/3.0,0.0,1.0)
clipped_A = np.clip((translated_A-0.1)/1.0,0.0,1.0)
clipped_S = translated_S

unity = elliptical.copy()
unity[:] = 1.0

color = np.transpose((clipped_C, np.zeros_like(clipped_C), np.zeros_like(clipped_C))) # Use all three colors...
color = np.transpose((np.zeros_like(clipped_C), np.zeros_like(clipped_C), clipped_A )) # Use all three colors...

#color = np.transpose((1-other, unity, 1-other)) # ...or use just green.
#color = np.transpose((1-spiral, 1-spiral, unity)) # ...or use just blue.
#color = np.transpose((unity, 1-elliptical, 1-elliptical)) # ...or use just red.

color = color.clip(0.0, 1.0)

p0 = plt.scatter(x, y, s=70, c=color)

#p0 = plt.scatter(x, y, s=50, c=color, cmap=plt.cm.coolwarm)

#plt.legend( (p0, p4), ('i=0', 'i=4'), 'upper right', shadow=True)

# Add labels to the plot.
plt.xlabel(r'$NII / H_{\alpha}$', fontsize=16)
plt.ylabel(r'$OIII / H_{\beta}$', fontsize=16)

plt.title('BPT of SDSS Galaxies', fontsize=16)

# Set the range for the plot.
plt.axis([-1.4, 0.3, -1.5,1.3])
#print(table_data.field(14))

# Now that everything is set, show the plot!
plt.show()

# PLOT.LY
#plot_url = py.plot_mpl(fig)
print
uniqObjId = f2(objid)
f = open('makelinks','w')
for i in np.arange(len(uniqObjId)):
    f.write('ln -s ' + str(uniqObjId[i]) + '.pdf ' + '{:0>4d}'.format(i) + '.pdf\n')
f.close()

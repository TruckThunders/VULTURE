import string
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
from numpy.fft import rfft, irfft
from scipy.stats import norm
from scipy import stats
np.set_printoptions(threshold=np.nan)

#===========================================================
#======================== Functions ========================
#===========================================================

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def group_consecutives(vals, step=0.0001):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = 0.
    for v in vals:
        v = round(v,5)
        expect = round(expect,5)
        if (v == expect) or (expect == 0.):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result

# ASYMMETRY CHECK
def asymmetry_check(line_redshifts, ion_spectrum):
    # Roll a ball around the line_redshifts to throw out false positives
    count = 0
    for detection in line_redshifts:
        leftup = 0
        leftdown = 0
        rightup = 0
        rightdown = 0
        # Define line center and establish max value
        center = np.where(redshift == find_nearest(redshift,detection))
        center = int(center[0])
        # Check for bigger max
        for i in range(-5,5):
            max = np.copy(scaled_MgII2796[center])
            maxcheck = np.copy(scaled_MgII2796[center+i])
            if (maxcheck > max):
                max = np.copy(maxcheck)    
        print detection
        print max 
        if ( max > 0.9 ):
            count += 1
            continue
        # range(-10,10) is how far to search left/right for asymmetry
        for i in range(-10,10):
            # Roll ball around max value and see if it checks out
            check = scaled_MgII2796[center + i]
            if i < 0:
                if check < max:
                    leftdown += 1
                if check > max: 
                    leftup += 1
            if i > 0:
                if check < max:
                    rightdown += 1
                if check > max:
                    rightup += 1
        # Delete lines based upon asymmetry around the line center
        if ((abs(leftup - rightup) > 4) & (abs(leftdown - rightdown) > 4)):
            line_redshifts = np.delete(line_redshifts,count)
            print "Deleted line (asymmetry):"
            print detection  
            count -= 1
        count += 1
    return line_redshifts

#================================================================
#======================== Initialization ========================
#================================================================

# Read in command line args
filename = str(sys.argv[1])

# Create lists
wav = []
flux = []
sig = []
err = []
cont = []

# Convert lists to numpy arrays
wav = np.array(wav)
flux = np.array(flux)
sig = np.array(sig)
err = np.array(err)
cont = np.array(cont)

# Transition Wavelengths
LYA = 1215.6
MgII2796 = 2796.353786 
MgII2803 = 2803.530982
CIV1548 = 1548.195
CIV1551 = 1550.770
LIGHTSPEED = 299792.458 # km/s

# Redshift resolution
delta_z = 0.00003 # = 9 km/s

# Read in spectra

print " "
print "========== Reading in Spectra =========="
print " "

# Ignore leading rows that are only 4 columns long
ignore = 0
with open(filename,'rb') as f:
    for line in f:
        check = len(line.split())
        if ( check == 5 ):
            ignore += 1
        else:
            print "Ignored first {} lines of code because they were dumb.".format(ignore)
            break
        
        # end else
    # end for
#end with
        
wav,dum,flux,sig,err,cont = np.loadtxt(filename,unpack=True,skiprows=ignore)

# Read in zem.dat
zem = np.loadtxt("zem.dat")

# Create flux arrays for data and error
norm_flux = flux/cont
norm_err = err/cont
norm_sig = sig/cont

#======================================================================
#======================== Set up Search Window ========================
#======================================================================

print " "
print "========== Initializing Search Arrays =========="
print " "

# Search redward of the Lya emission
if ( (LYA * (1. + zem)) > wav[0]):
    blue_limit = LYA * (1. + zem)
else:
    blue_limit = wav[0]

if ( (MgII2803 * (1. + zem)) < wav[-1]):
    red_limit = MgII2803 * (1. + zem)
else:
    red_limit = wav[-1]

search_flux = norm_flux[np.where((wav > blue_limit) & (wav < red_limit))]
search_wav = wav[np.where((wav > blue_limit) & (wav < red_limit))]
search_err = norm_err[np.where((wav > blue_limit) & (wav < red_limit))]
search_sig = norm_sig[np.where((wav > blue_limit) & (wav < red_limit))]

# Clean cosmic rays
search_flux[np.where(search_flux > 1.1)] = 1.
search_err[np.where(search_flux > 1.1)] = 0.

# Avoid Telluric features
# 9300 - 9630 feature
search_flux[np.where((search_wav > 9300.) & (search_wav < 9630.))] = 1.

# 7594 - 7700 feature
search_flux[np.where((search_wav > 7594.) & (search_wav < 7700.))] = 1.

# 6868 - 6932 feature
search_flux[np.where((search_wav > 6868.) & (search_wav < 6932.))] = 1.

# Remove negative errors
search_err[np.where(search_err < 0.)] = 0.

# Create arrays to hold the cross-correlation result (to be plotted/analyzed)
result_MgII2796 = []
np.array(result_MgII2796)
result_MgII2803 = []
np.array(result_MgII2803)
result_CIV1548 = []
np.array(result_CIV1548)
result_CIV1551 = []
np.array(result_CIV1551)

# Create arrays to hold the noise cross-correlation
noise_MgII2796 = []
np.array(noise_MgII2796)
noise_MgII2803 = []
np.array(noise_MgII2803)
noise_CIV1548 = []
np.array(noise_CIV1548)
noise_CIV1551 = []
np.array(noise_CIV1551)

# Search redshift array from z = 0 to z = zem
lower_lim = ( blue_limit / MgII2796 ) - 1.
upper_lim = ( red_limit / MgII2803 ) - 1.

# Make redshift array with delta_z resolution
redshift = arange(round(lower_lim,5),round(upper_lim,5),delta_z)

#==================================================================
#======================== Do the Filtering ========================
#==================================================================

# Define filter size (filt_size = filter width / 2)
filt_size = 5

# Find min value of correlation
min_flux = np.ones(search_flux.size+(2.*filt_size))
min_flux[0:search_flux.size] = np.copy(search_flux)
test_filt = np.ones(min_flux.size)
test_filt[search_flux.size:min_flux.size] = 0.

print "Min Value Established."
min_value = np.correlate(min_flux,test_filt) - (2.*filt_size)
print min_value

# Find max value of correlation
z = np.copy(lower_lim)
search_filt = np.ones(search_flux.size)
MgII2796central = find_nearest(search_wav,MgII2796*(1.+z))
MgII2796central = np.where(search_wav == MgII2796central)
MgII2796central = MgII2796central[0]
search_filt[MgII2796central-filt_size:MgII2796central+filt_size] = 0.
max_value = np.correlate(search_flux,search_filt)

print "Max Value Established."
max_value = np.correlate(search_flux,search_filt)
print max_value
    
# Loop redshift
for z in redshift:

    # Correlate with MgII2796 line
    search_filt = np.ones(search_flux.size)
    MgII2796central = find_nearest(search_wav,MgII2796*(1.+z))
    MgII2796central = np.where(search_wav == MgII2796central)
    MgII2796central = MgII2796central[0]
    search_filt[MgII2796central-filt_size:MgII2796central+filt_size] = 0.
    result_MgII2796.append(np.correlate(search_flux,search_filt))

    # Correlate with MgII2803 line
    search_filt = np.ones(search_flux.size)
    MgII2803central = find_nearest(search_wav,MgII2803*(1.+z))
    MgII2803central = np.where(search_wav == MgII2803central)
    MgII2803central = MgII2803central[0]
    search_filt[MgII2803central-filt_size:MgII2803central+filt_size] = 0.
    result_MgII2803.append(np.correlate(search_flux,search_filt))

# Scale the power spectrum to something we can work with (0 to 1)
scaled_MgII2796 = (result_MgII2796 - min_value) / (max_value - min_value)
scaled_MgII2803 = (result_MgII2803 - min_value) / (max_value - min_value)

scaled_MgII2796 = scaled_MgII2796.flatten()
scaled_MgII2803 = scaled_MgII2803.flatten()


#==================================================================
#================ Calculate a power spectrum error ================
#==================================================================

# Chunk the power spectrum
chunk_size = (0.03 / delta_z)
index1 = (chunk_size/2.) - 1
index2 = chunk_size/2.
index3 = chunk_size + 1
index4 = (chunk_size/2.) + 1

chunk1 = []
chunk1 = np.array(chunk1)
chunk2 = []
chunk2 = np.array(chunk2)
noise_MgII2796 = np.zeros(scaled_MgII2796.size)
noise_MgII2803 = np.copy(noise_MgII2796)

# Pad the signal array by mirroring left and right side
working_MgII2796 = np.zeros(scaled_MgII2796.size + chunk_size)
working_MgII2796[0:index1] = scaled_MgII2796[0:index1]
working_MgII2796[index2:(len(scaled_MgII2796)+index2-1)] = scaled_MgII2796
working_MgII2796[(len(scaled_MgII2796)+index4):(len(scaled_MgII2796)+index3)] = scaled_MgII2796[(len(scaled_MgII2796)-index2):-1]

working_MgII2803 = np.zeros(scaled_MgII2803.size + chunk_size)
working_MgII2803[0:index1] = scaled_MgII2803[0:index1]
working_MgII2803[index2:(len(scaled_MgII2803)+index2-1)] = scaled_MgII2803
working_MgII2803[(len(scaled_MgII2803)+index4):(len(scaled_MgII2803)+index3)] = scaled_MgII2803[(len(scaled_MgII2803)-index2):-1]

# In the error array, erase the null chunks by filling them with mirrored data
# MgII2796
zerovals_MgII2796 = np.where((working_MgII2796 < 1.e-5) & (working_MgII2796 > -1.e-5))
zerovals_MgII2796 = zerovals_MgII2796[0]
gaps_MgII2796 = group_consecutives(zerovals_MgII2796,1)

# Skip over any chunk less than 5 points in length for MgII2796
for spots in gaps_MgII2796:
    if len(spots) < 5:
        continue
    working_MgII2796[spots] = np.copy(working_MgII2796[np.subtract(spots,len(spots)).astype(int)])

# MgII2803
zerovals_MgII2803 = np.where((working_MgII2803 < 1.e-5) & (working_MgII2803 > -1.e-5))
zerovals_MgII2803 = zerovals_MgII2803[0]
gaps_MgII2803 = group_consecutives(zerovals_MgII2803,1)

# Skip over any chunk less than 5 points in length for MgII2803
for spots in gaps_MgII2803:
    if len(spots) < 5:
        continue
    working_MgII2803[spots] = np.copy(working_MgII2803[np.subtract(spots,len(spots)).astype(int)])

# Perform sigma-clipping and generate power-spectrum errors at each redshift
for counter,z in enumerate(redshift):
    flux_index = counter + chunk_size
    chunk1 = working_MgII2796[counter:flux_index]
    chunk1[np.where((chunk1 < 1.e-5) & (chunk1 > -1.e-5))] = None
    chunk2 = working_MgII2803[counter:flux_index]
    chunk2[np.where((chunk2 < 1.e-5) & (chunk2 > -1.e-5))] = None
    # Sigma-clip x times, where x --> range(0,x)
    for x in range(0,5):
        # Sigma-clip chunk1
        onesig = stats.nanstd(chunk1)
        threesig = 3. * onesig
        median = stats.nanmedian(chunk1)
        outliers = np.where(chunk1 > (threesig + median))
        for point in outliers[0]:
            chunk1[point] = None
            i = 1
            j = 1
            # Roll ball left
            while True:
                if ((point-i) < 0):
                    break
                if chunk1[point-i] < (onesig + median):
                    break
                chunk1[point-i] = None
                i += 1
            # Roll ball right
            while True:
                if ((point+j) > (len(chunk1)-1)):
                    break
                if chunk1[point+j] < (onesig + median):
                    break
                chunk1[point+j] = None
                j += 1
        # Sigma-clip chunk2
        onesig = stats.nanstd(chunk2)
        threesig = 3. * onesig
        median = stats.nanmedian(chunk2)
        outliers = np.where(chunk2 > (threesig + median))
        for point in outliers[0]:
            chunk2[point] = None
            i = 1
            j = 1
            # Roll ball left
            while True:
                if ((point-i) < 0):
                    break
                if chunk2[point-i] < (onesig + median):
                    break
                chunk2[point-i] = None
                i += 1
            # Roll ball right
            while True:
                if ((point+j) > (len(chunk2)-1)):
                    break
                if chunk2[point+j] < (onesig + median):
                    break
                chunk2[point+j] = None
                j += 1
        # End of line-sigma-clip for chunks
        
    noise_MgII2796[counter] += (stats.nanstd(chunk1) + stats.nanmedian(chunk1))
    noise_MgII2803[counter] += (stats.nanstd(chunk2) + stats.nanmedian(chunk2))

# Don't evaluate the power spectrum noise where there are gaps in the spectra
badvals_MgII2796 = np.where((scaled_MgII2796 < 1.e-5) & (scaled_MgII2796 > -1.e-5))
badvals_MgII2803 = np.where((scaled_MgII2803 < 1.e-5) & (scaled_MgII2803 > -1.e-5))

noise_MgII2796[badvals_MgII2796[0]] = 0.
noise_MgII2803[badvals_MgII2803[0]] = 0.

sigma_MgII2796 = np.divide(scaled_MgII2796,noise_MgII2796)
sigma_MgII2803 = np.divide(scaled_MgII2803,noise_MgII2803)

# Define criterion for detecting an absorption line system 
criterion_MgII = np.where((sigma_MgII2796 > 5) & (sigma_MgII2803 > 3))
final_MgII2796 = scaled_MgII2796[criterion_MgII[0]]
final_MgII2803 = scaled_MgII2803[criterion_MgII[0]]
final_redshift = redshift[criterion_MgII[0]]

# Group the adjacent redshifts to find the detection centers
grouped_redshifts = group_consecutives(final_redshift,delta_z)

# Find line centers
line_redshifts=[]
for lines in grouped_redshifts:
    center = np.mean(lines)
    line_redshifts.append(center)

print " "
print "========== Lines detected. Culling bad lines =========="
print " "

# LINE RATIO CHECK
# Throw out false positives based upon line ratio
count = 0
for detection in line_redshifts:
    center = np.where(redshift == find_nearest(redshift,detection))
    center = int(center[0])
    # Delete lines based upon line ratio
    for offset in range(3):
        if (((scaled_MgII2796[center+offset]+0.002) < scaled_MgII2803[center+offset]) | 
            ((scaled_MgII2796[center-offset]+0.002) < scaled_MgII2803[center-offset])):
            line_redshifts = np.delete(line_redshifts,count)
            print "Deleted line (MgII2803 >> MgII2796):"
            print detection
            count -= 1
            break
        elif ((scaled_MgII2803[center]/scaled_MgII2796[center]) < 0.5):
            line_redshifts = np.delete(line_redshifts,count)
            print "Deleted line (Bad doublet ratio):"
            print detection
            count -= 1
            break
    count += 1

# ASYMMETRY CHECK
# Roll a ball around the line_redshifts to throw out false positives
count = 0
for detection in line_redshifts:
    leftup = 0
    leftdown = 0
    rightup = 0
    rightdown = 0
    # Define line center and establish max value
    center = np.where(redshift == find_nearest(redshift,detection))
    center = int(center[0])
    # Check for bigger max
    for i in range(-5,5):
        max = np.copy(scaled_MgII2796[center])
        maxcheck = np.copy(scaled_MgII2796[center+i])
        if (maxcheck > max):
            max = np.copy(maxcheck)
    
    print detection
    print max 
    if ( max > 0.9 ):
        count += 1
        continue
    # range(-10,10) is how far to search left/right for asymmetry
    for i in range(-10,10):
        # Roll ball around max value and see if it checks out
        check = scaled_MgII2796[center + i]
        if i < 0:
            if check < max:
                leftdown += 1
            if check > max: 
                leftup += 1
        if i > 0:
            if check < max:
                rightdown += 1
            if check > max:
                rightup += 1
    # Delete lines based upon asymmetry around the line center
    if ((abs(leftup - rightup) > 4) & (abs(leftdown - rightdown) > 4)):
        line_redshifts = np.delete(line_redshifts,count)
        print "Deleted line (asymmetry):"
        print detection  
        count -= 1
    count += 1

print "-------------------------------------"
print "Detected lines:"
print line_redshifts
print "-------------------------------------"

# Set up array to plot detected lines
detected_MgII2796=[]
for z in line_redshifts:
    temp_z = find_nearest(final_redshift,z)
    detected_MgII2796.append(final_MgII2796[np.where(temp_z == final_redshift)])

# Plot the things

plot(redshift,scaled_MgII2796)
plot(redshift,scaled_MgII2803)
plot(final_redshift,final_MgII2796,'ko')
plot(line_redshifts,detected_MgII2796,'ro')
plot(redshift,noise_MgII2796,'r-')
plot(redshift,3.*noise_MgII2796)
plot(redshift,5.*noise_MgII2796)
plot(redshift,np.zeros(noise_MgII2796.size),'--')
show()

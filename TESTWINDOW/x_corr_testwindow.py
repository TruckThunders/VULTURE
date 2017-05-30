import string
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
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

def asymmetry_check(line_redshifts, ion_spectrum, redshift):
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
            max = np.copy(ion_spectrum[center])
            maxcheck = np.copy(ion_spectrum[center+i])
            if (maxcheck > max):
                max = np.copy(maxcheck)    
        if ( max > 0.9 ):
            count += 1
            continue
        # range(-10,10) is how far to search left/right for asymmetry
        for i in range(-10,10):
            # Roll ball around max value and see if it checks out
            check = ion_spectrum[center + i]
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

# ion_spectrum1 corresponds to the ion with the larger ionization potential
def line_ratio_check(line_redshifts, ion_spectrum1, ion_spectrum2,redshift,doublet_ratio):
    # Throw out false positives based upon line ratio
    count = 0
    for detection in line_redshifts:
        center = np.where(redshift == find_nearest(redshift,detection))
        center = int(center[0])
        # Delete lines based upon line ratio
        for offset in range(3):
            if (((ion_spectrum1[center+offset]+0.002) < ion_spectrum2[center+offset]) | 
                ((ion_spectrum1[center-offset]+0.002) < ion_spectrum2[center-offset])):
                line_redshifts = np.delete(line_redshifts,count)
                print "Deleted line (Main Transition >> Secondary):"
                print detection
                count -= 1
                break
            elif (((ion_spectrum2[center+offset]/ion_spectrum1[center+offset]) < doublet_ratio) |
                  ((ion_spectrum2[center-offset]/ion_spectrum1[center-offset]) < doublet_ratio)):
                line_redshifts = np.delete(line_redshifts,count)
                print "Deleted line (Bad doublet ratio):"
                print detection
                count -= 1
                break
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

# Transition Wavelengths and their power/noise arrays
species = {}
power = {}
noise = {}
working = {}
sigma = {}
detected = {}

#species['FeII'] = {}
#species['FeII']['1608'] = 1608.4507986
#species['FeII']['1611'] = 1611.2004037
#species['FeII']['2382'] = 2382.7640777
#species['FeII']['2600'] = 2600.1721140
#species['FeII']['2344'] = 2344.2127470
#species['FeII']['2586'] = 2586.6493120

# species['AlII'] = {}
# species['AlII']['1670'] = 1670.7886100 

# species['AlIII'] = {}
# species['AlIII']['1862'] = 1862.798581
# species['AlIII']['1854'] = 1854.708966

#species['MgI'] = {}
#species['MgI']['2852'] = 2852.9627970

species['CIV'] = {}
species['CIV']['1548'] = 1548.195
power['1548'] = []
noise['1548'] = []
working['1548'] = []
sigma['1548'] = []
detected['1548'] = []
species['CIV']['1551'] = 1550.770
power['1551'] = []
noise['1551'] = []
working['1551'] = []
sigma['1551'] = []
detected['1551'] = []

#species['MgII'] = {}
#species['MgII']['2796'] = 2796.3537860
#power['2796'] = []
#noise['2796'] = []
#working['2796'] = []
#sigma['2796'] = []
#detected['2796'] = []
#species['MgII']['2803'] = 2803.5309820
#power['2803'] = []
#noise['2803'] = []
#working['2803'] = []
#sigma['2803'] = []
#detected['2803'] = []

# Define Constants
LIGHTSPEED = 299792.458 # km/s
LYA = 1215.6 # angstrom

# Redshift resolution
#delta_z = 0.00003 # = 9 km/s
delta_z = 0.0005 # = 9 km/s

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

# Find reddest and bluest transitions
reddest_transition = 0.
bluest_transition = 100000.
for specie in species:
    for transition in species[specie]:
        rest_wavelength = species[specie][transition]
        if ( rest_wavelength > reddest_transition ):
            reddest_transition = np.copy(rest_wavelength)
        if ( rest_wavelength < bluest_transition ):
            bluest_transition = np.copy(rest_wavelength)

# Choose either the reddest_transition or the end of the spectrum
if ( (reddest_transition * (1. + zem)) < wav[-1]):
    red_limit = reddest_transition * (1. + zem)
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

# Search redshift array from z = 0 to z = zem
lower_lim = ( blue_limit / bluest_transition ) - 1.
upper_lim = ( red_limit / reddest_transition ) - 1.

if ( lower_lim > upper_lim ):
    sys.exit("No redshift range for all transitions requested.")

# Make redshift array with delta_z resolution
redshift = arange(round(lower_lim,5),round(upper_lim,5),delta_z)

#==================================================================
#======================== Do the Filtering ========================
#==================================================================

# Define filter size (filt_size = filter width / 2)
filt_size = 3

# Find min value of correlation
min_flux = np.ones(search_flux.size+(2.*filt_size))
min_flux[0:search_flux.size] = np.copy(search_flux)
test_filt = np.ones(min_flux.size)
test_filt[search_flux.size:min_flux.size] = 0.

print "Min Value Established."
min_value = np.correlate(min_flux,test_filt) - (2.*filt_size)
print min_value

# Find max value of correlation
test_filt = np.ones(search_flux.size)
central = find_nearest(search_wav,blue_limit)
central = np.where(search_wav == central)
central = central[0]
test_filt[central-filt_size:central+filt_size] = 0.
max_value = np.correlate(search_flux,test_filt)

print "Max Value Established."
max_value = np.correlate(search_flux,test_filt)
print max_value

#shifts = [-0.00009,-0.00006,-0.00003,0.00003,0.00006,0.00009,0.00042]

for shift in range(0,10):
    power['1548'] = []
    power['1551'] = []

    # Loop redshift
    for z in redshift:
        for specie in species:
            for transition in species[specie]:
                # Correlate at each transition wavelength
                rest_wavelength = species[specie][transition]
                search_filt = np.ones(search_flux.size)
                central = find_nearest(search_wav,rest_wavelength*(1.+z))
                central = np.where(search_wav == central)
                central = central[0]
                search_filt[central-filt_size:central+filt_size] = 0.
                power[transition].append(np.correlate(search_flux,search_filt))

    search_flux = np.insert(search_flux,0,1)
    #search_wav = np.insert(search_wav,0,(search_wav[0] - (search_wav[1] - search_wav[0])))

    # Scale the power spectrum to something we can work with (0 to 1)
    scaledpower = {}
    for transition in power:
        power[transition] = np.array(power[transition])
        power[transition] = power[transition].flatten()
        scaledpower[transition] = (power[transition] - min_value) / (max_value - min_value)
        
    power_z = find_nearest(redshift,1.97386)
    power_index = np.where(redshift == power_z)
    power_index = power_index[0]
    power_value = scaledpower[transition][power_index]

    print "Central = {}.".format(central)
    print "Power value = {}.".format(power_value)

    # Plot the things
    #figure()
    if ( specie == 'MgII' ):
        plot(redshift,scaledpower['2796'],label="2796")
        plot(redshift,scaledpower['2803'],label="2803")
        transition = '2796'
    if ( specie == 'CIV' ):
        plot(np.subtract(redshift,shift*0.00003),scaledpower['1548'],label="1548, shift={}".format(shift))
        #plot(redshift,scaledpower['1551'],label="1551")
        #transition = '1548'
    
legend()
show()



#==================================================================
#========================== File Output ===========================
#==================================================================

#f = open('detections.dat','w')
#for specie in species:
#    f.write("--------------------------------------------------------------------------")
#    f.write("Detected {} lines:".format(specie))
#    f.write(line_redshifts[specie])
#    f.write("--------------------------------------------------------------------------")
#f.close()

import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#from pylab import *
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

# Insert fake lines to test detection sensitivity ***Height is passed by reference
def insert_lines(spectrum,wavelengths,resolution,z,tranname,eqwidth):
    rest_wavs = []
    if (tranname == 'MgII'):
        rest_wavs.append(2796.353786)
        rest_wavs.append(2803.530982)
    if (tranname == 'CIV'):
        rest_wavs[0] = 1548.195
        rest_wavs[1] = 1550.770

    for lambda_r in rest_wavs:
        testwidth = 0.
        # Calculate gaussian parameters
        b = lambda_r * (1+z)
        c = (b / resolution) * (1. / (2.*np.sqrt(2.*np.log(2))))
        if (lambda_r == 2796.353786):
            a = eqwidth / (c * np.sqrt(2. * np.pi))
        if (lambda_r == 2803.530982):
            a = eqwidth / (1.5 * c * np.sqrt(2. * np.pi))
        npix = (resolution * 4.) / b
        half_width = int(np.ceil(npix/2.))

        middle = find_nearest(wavelengths,b)
        middle = np.where(wavelengths == middle)
        middle = middle[0]

        #print "b = {}. c = {}. a = {}.".format(b,c,a)

        for i in range(-half_width,half_width):
            spectrum[middle+i] = 1.
            spectrum[middle+i] -= a * np.exp(-((wavelengths[middle+i] - b)**2.)/(2.*c*c))
            testwidth += (1.-spectrum[middle+i])*(wavelengths[middle]-wavelengths[middle-1])

    return (a,c,spectrum)

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

#species['CIV'] = {}
#species['CIV']['1548'] = 1548.195
#power['1548'] = []
#noise['1548'] = []
#working['1548'] = []
#sigma['1548'] = []
#detected['1548'] = []
#species['CIV']['1551'] = 1550.770
#power['1551'] = []
#noise['1551'] = []
#working['1551'] = []
#sigma['1551'] = []
#detected['1551'] = []

species['MgII'] = {}
species['MgII']['2796'] = 2796.3537860
power['2796'] = []
noise['2796'] = []
working['2796'] = []
sigma['2796'] = []
detected['2796'] = []
species['MgII']['2803'] = 2803.5309820
power['2803'] = []
noise['2803'] = []
working['2803'] = []
sigma['2803'] = []
detected['2803'] = []

# Define Constants
LIGHTSPEED = 299792.458 # km/s
LYA = 1215.6 # angstrom

# Spectral resolution
spectral_res = 40000.

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

print "Searching redward of Lya."

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
#

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

# TO DO: Make different redshift arrays for each specie
# Make redshift array with delta_z resolution
redshift = np.arange(round(lower_lim,5),round(upper_lim,5),delta_z)
#redshift = np.arange(0.1,4,delta_z)


#================ INSERT FAKE LINES AND DETERMINE DETECTION LIMIT =======================

num_fake_lines = 0

#print "--------------------------------"
#print "Inserting fake absorption lines:"

# Define synthetic line redshift to insert
fake_redshift = np.arange(lower_lim+0.1,upper_lim-0.1,((upper_lim-0.1)-(lower_lim-0.1))/20.)
height=np.zeros(len(fake_redshift))
width=np.zeros(len(fake_redshift))

#print fake_redshift

# Split filename to remove .data and append .ewspec
fullname = []
fullname = filename.split('.')
ewname = fullname[0]+'.ewspec'

# Determine equivalent width detection limit
eq_wav,eq_width,eq_fivesigma = np.loadtxt(ewname,unpack=True,skiprows=ignore)

# Copy search_flux array 
fake_flux = np.copy(search_flux)

# Insert fake lines to detect
for shifts in fake_redshift:
# Define line center and establish max value
    center = np.where(eq_wav == find_nearest(eq_wav,2796.35378*(1.+shifts)))
    center = int(center[0])
    detectionlimit = np.fabs(np.mean(eq_fivesigma[center-250:center+250]))
    #print "Detection Limit = {}.".format(detectionlimit)
    if ( detectionlimit > 1 ) :
        fake_redshift = np.delete(fake_redshift,num_fake_lines)
        height = np.delete(height,num_fake_lines)
        continue
    #search_flux = insert_lines(search_flux,search_wav,40000.,shifts,'MgII',3.*detectionlimit,num_fake_lines)
    height[num_fake_lines],width[num_fake_lines],fake_flux = insert_lines(fake_flux,search_wav,spectral_res,shifts,'MgII',3.*detectionlimit)
    num_fake_lines += 1

#print "--------------------------------"

#================== END OF FAKE LINE INSERTION ROUTINES ==============================

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
min_value = np.correlate(min_flux,test_filt) - (2.*filt_size)

# Find fake min for EW conversion
fake_min = np.ones(fake_flux.size+(2.*filt_size))
fake_min[0:fake_flux.size] = np.copy(fake_flux)
test_fake = np.ones(fake_min.size)
test_fake[fake_flux.size:fake_min.size] = 0.
fake_min_value = np.correlate(fake_min,test_fake) - (2.*filt_size)

print "Min Value Established."
print min_value

# Find max value of correlation
test_filt = np.ones(search_flux.size)
central = find_nearest(search_wav,blue_limit)
central = np.where(search_wav == central)
central = central[0]
test_filt[central-filt_size:central+filt_size] = 0.
max_value = np.correlate(search_flux,test_filt)

# Find fake max for EW conversion
test_fake = np.ones(fake_flux.size)
central = find_nearest(search_wav,blue_limit)
central = np.where(search_wav == central)
central = central[0]
test_fake[central-filt_size:central+filt_size] = 0.
fake_max_value = np.correlate(fake_flux,test_fake)

print "Max Value Established."
print max_value

fake_power = []
# Loop over fake lines to get fake values for EW
for fake_z in fake_redshift:
    rest_wavelength = species['MgII']['2796']
    fake_filt = np.ones(fake_flux.size)
    central = find_nearest(search_wav,rest_wavelength*(1.+fake_z))
    central = np.where(search_wav == central)
    central = central[0]
    fake_filt[central-filt_size:central+filt_size] = 0.
    fake_power.append(np.correlate(fake_flux,fake_filt))

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

# Scale the power spectrum to something we can work with (0 to 1)
scaledpower = {}
for transition in power:
    power[transition] = np.array(power[transition])
    power[transition] = power[transition].flatten()
    scaledpower[transition] = (power[transition] - min_value) / (max_value - min_value)

# Scale fake power spectrum
fakescaled = []
fake_power = np.array(fake_power)
fake_power = fake_power.flatten()
fakescaled = (fake_power - fake_min_value) / (fake_max_value - fake_min_value)

#==================================================================
#================ Calculate a power spectrum error ================
#==================================================================

# Chunk the power spectrum
chunk_size = (0.03 / delta_z)
index1 = (chunk_size/2.) - 1
index2 = chunk_size/2.
index3 = chunk_size + 1
index4 = (chunk_size/2.) + 1

chunk = {}

# Create noise arrays
for transition in noise:
    noise[transition] = np.zeros(scaledpower[transition].size)

# Pad the signal array by mirroring left and right side
for transition in working:
    working[transition] = np.zeros(scaledpower[transition].size + chunk_size)
    working[transition][0:index1] = scaledpower[transition][0:index1]
    working[transition][index2:(len(scaledpower[transition])+index2-1)] = scaledpower[transition]
    working[transition][(len(scaledpower[transition])+index4):(len(scaledpower[transition])+index3)] = scaledpower[transition][(len(scaledpower[transition])-index2):-1]
    
    # In the error array, erase the null chunks by filling them with mirrored data
    zerovals = np.where((working[transition] < 1.e-5) & (working[transition] > -1.e-5))
    zerovals = zerovals[0]
    gaps = group_consecutives(zerovals,1)

    # Skip over any chunk less than 5 points in length for MgII2796
    for spots in gaps:
        if len(spots) < 5:
            continue
        working[transition][spots] = np.copy(working[transition][np.subtract(spots,len(spots)).astype(int)])
        
# Perform sigma-clipping and generate power-spectrum errors at each redshift
for counter,z in enumerate(redshift):
    flux_index = counter + chunk_size
    for transition in working:
        chunk[transition] = working[transition][counter:flux_index]
        chunk[transition][np.where((chunk[transition] < 1.e-5) & (chunk[transition] > -1.e-5))] = None
        # Sigma-clip x times, where x --> range(0,x)
        for x in range(0,5):
            # Sigma-clip chunk
            onesig = stats.nanstd(chunk[transition])
            threesig = 3. * onesig
            median = stats.nanmedian(chunk[transition])
            outliers = np.where(chunk[transition] > (threesig + median))
            for point in outliers[0]:
                chunk[transition][point] = None
                i = 1
                j = 1
                # Roll ball left
                while True:
                    if ((point-i) < 0):
                        break
                    if chunk[transition][point-i] < (onesig + median):
                        break
                    chunk[transition][point-i] = None
                    i += 1
                # Roll ball right
                while True:
                    if ((point+j) > (len(chunk[transition])-1)):
                        break
                    if chunk[transition][point+j] < (onesig + median):
                        break
                    chunk[transition][point+j] = None
                    j += 1
        
        noise[transition][counter] += (stats.nanstd(chunk[transition]) + stats.nanmedian(chunk[transition]))

# Don't evaluate the power spectrum noise where there are gaps in the spectra
for transition in noise:
    sigma[transition] = np.divide(scaledpower[transition],noise[transition])

# ======================== FAKE LINE/DETECTION LIMIT STUFF =================================
# Determine scaling between power and gaussian height
height = np.array(height)

A = np.vstack([height,np.ones(len(height))]).T
m,b = np.linalg.lstsq(A,fakescaled)[0]
print "m = {}   b = {}.".format(m,b)

#plot(height,fakescaled,"ro")
#plot(height,(m*height)+b)
#show()
    
#==================================================================
#============= Find the Lines and Filter the Bad ones =============
#==================================================================

criterion = {}
final = {}
final_redshift = {}
line_redshifts = {}
line_EWlim = {}
grouped_redshifts = {}
EWlimit = {}

for specie in species:
    # Define criterion for detecting an absorption line system 
    if ( specie == 'MgII' ):
        # Calculate equivalent width detection limit
        EWlimit[specie] = []
        EWlimit[specie] = np.array(EWlimit)
        # Calculate the new b in the gaussian for an unresolved line at the detection limit
        newb = species[specie]['2796'] * (1.+redshift)
        EWwidth = (newb / spectral_res) * (1. / (2.*np.sqrt(2.*np.log(2))))
        EWlimit[specie] = (np.multiply(5.*noise['2796'],EWwidth) * np.sqrt(2.*np.pi)) / m
        # Define final detection redshifts
        criterion[specie] = np.where((sigma['2796'] > 5) & (sigma['2803'] > 3))
        final['2796'] = scaledpower['2796'][criterion[specie][0]]
        final['2803'] = scaledpower['2803'][criterion[specie][0]]
        final_redshift[specie] = redshift[criterion[specie][0]]
    if ( specie == 'CIV' ):
        # Calculate equivalent width detection limit
        EWlimit[specie] = []
        EWlimit[specie] = np.array(EWlimit)
        # Calculate the new b in the gaussian for an unresolved line at the detection limit
        newb = species[specie]['1548'] * (1.+redshift)
        EWwidth = (newb / spectral_res) * (1. / (2.*np.sqrt(2.*np.log(2))))
        EWlimit[specie] = (np.multiply(5.*noise['1548'],EWwidth) * np.sqrt(2.*np.pi)) / m
        # Define final detection redshifts
        criterion[specie] = np.where((sigma['1548'] > 5) & (sigma['1551'] > 3))
        final['1548'] = scaledpower['1548'][criterion[specie][0]]
        final['1551'] = scaledpower['1551'][criterion[specie][0]]
        final_redshift[specie] = redshift[criterion[specie][0]]    
    
    # Group the adjacent redshifts to find the detection centers
    grouped_redshifts[specie] = group_consecutives(final_redshift[specie],delta_z)

    # Find line centers
    line_redshifts[specie] = []
    for lines in grouped_redshifts[specie]:
        center = np.mean(lines)
        line_redshifts[specie].append(center)

    print " "
    print "========== {} Lines detected. Culling bad lines ==========".format(specie)
    print " "

    # THROW OUT BAD DETECTIONS
    # Line Ratio Check
    if ( specie == 'MgII' ):
        line_redshifts[specie] = line_ratio_check(line_redshifts[specie],scaledpower['2796'],scaledpower['2803'],redshift,0.4)
    if ( specie == 'CIV' ):
        line_redshifts[specie] = line_ratio_check(line_redshifts[specie],scaledpower['1548'],scaledpower['1551'],redshift,0.3)

    # Figure out how many FAKE lines were detected
    #fake_detect = 0
    #for lines in line_redshifts[specie]:
    #    for fake_lines in fake_redshift:            
    #        if ( (lines < (fake_lines+0.0001)) & (lines > (fake_lines-0.0001))):
    #            print "Real line {} matched with fake line {}.".format(lines,fake_lines)
    #            fake_detect += 1
    #            break

    #for transition in scaledpower:
    # ASYMMETRY CHECK - Second arg is the ion spectrum you're looking at
    #line_redshifts[specie] = asymmetry_check(line_redshifts[specie],scaledpower[transition],redshift)

    # Calculate EW Detection limit for each real detection
    line_EWlim[specie] = []
    for lines in line_redshifts[specie]:
        fake_centers = np.where(find_nearest(redshift,lines) == redshift)
        fake_centers = fake_centers[0]
        line_EWlim[specie].append(EWlimit[specie][fake_centers])

    line_EWlim[specie] = np.array(line_EWlim[specie])
    line_EWlim[specie] = line_EWlim[specie].flatten()

    print "--------------------------------------------------------------------------"
    print "Detected {} lines:".format(specie)
    print line_redshifts[specie]
    print "Detection limit for each detection (Angstrom):"
    print line_EWlim[specie]
    print "--------------------------------------------------------------------------"

    # Set up array to plot detected lines
    for z in line_redshifts[specie]:
        temp_z = find_nearest(final_redshift[specie],z)
        for transition in final:
            detected[transition].append(final[transition][np.where(temp_z == final_redshift[specie])])

    # Plot the things
    #figure()
    #if ( specie == 'MgII' ):
    #    plot(redshift,scaledpower['2796'],label="2796")
    #    plot(redshift,scaledpower['2803'],label="2803")
    #    transition = '2796'
    #if ( specie == 'CIV' ):
    #    plot(redshift,scaledpower['1548'],label="1548")
    #    plot(redshift,scaledpower['1551'],label="1551")
    #    transition = '1548'

    #plot(final_redshift[specie],final[transition],'ko')
    #plot(line_redshifts[specie],detected[transition],'ro')
    #plot(redshift,noise[transition],'r-')
    #plot(redshift,3.*noise[transition])
    #plot(redshift,5.*noise[transition])
    #plot(redshift,np.zeros(noise[transition].size),'--')
    #legend()

#show()

#plot(redshift,EWlimit)
#xlabel('Redshift')
#ylabel('W$_r$ Detection Limit ($\AA$)')
#show()

#==================================================================
#========================== File Output ===========================
#==================================================================

# Output line detections
f = open('detections.dat','w')
for specie in species:
    f.write("--------------------------------------------------------------------------\n")
    f.write("Detected {} lines:\n".format(specie))
    for item in line_redshifts[specie]:
        f.write("%s\n" % item)
    f.write("--------------------------------------------------------------------------\n")
f.close()

# Read in atoms.dat to get atomic data
transition_library = np.genfromtxt('../atoms_chris.dat',unpack=True,dtype=None)

velwindow = 500. #km/s
velocity = {}
vel_flux = {}
vel_err = {}
ax = ['ax1','ax2','ax3','ax4','ax5','ax6','ax7','ax8','ax9','ax10']
k = 0
i = 1
n = 0

blue_limit = 3100.
red_limit = 9900.

savedtran = "First"

# Output velocity plots for detections
for specie in species:
    # Detected line redshifts
    for detection_redshift in line_redshifts[specie]:
        if not os.path.exists(str(detection_redshift)):
            os.makedirs(str(detection_redshift))
        detection_redshift = float(detection_redshift)
        i = 0
        n = 0
        axis_index = 0
        plt.figure(figsize=(8.5,11))
        for transition in transition_library:
            lambda_r = transition[1] * (1. + detection_redshift)
            lambda_min = (-(velwindow/LIGHTSPEED) + 1)*lambda_r
            lambda_max = ((velwindow/LIGHTSPEED) + 1)*lambda_r
            index_min = np.where(wav == find_nearest(wav,lambda_min))
            index_max = np.where(wav == find_nearest(wav,lambda_max))
            index_min = index_min[0]
            index_max = index_max[0]
            velocity[transition[0]] = np.zeros(index_max - index_min)
            vel_flux[transition[0]] = np.zeros(index_max - index_min)
            vel_err[transition[0]] = np.zeros(index_max - index_min)
            j = 0
            for indices in range(index_min,index_max):
                velocity[transition[0]][j] = (LIGHTSPEED * (wav[indices] * wav[indices] - lambda_r * lambda_r)) / (wav[indices] * wav[indices] + lambda_r * lambda_r)
                vel_flux[transition[0]][j] = norm_flux[indices]
                vel_err[transition[0]][j] = norm_err[indices]
                j += 1
            # End velocity for loop
            if ( i % 5 == 0 ):
                if ( i != 0 ):
                    n = 1
                if ( i == 10 ):
                    plt.savefig("{}/main_trans.eps".format(detection_redshift),format='eps')
                    plt.close('all')
                    axis_index = 0
                    n = 0
                    plt.figure(figsize=(8.5,11))
                savedtran = transition[0]
                k = 0
            ax[axis_index] = plt.subplot2grid((5,2),(k,n))
            if ( i % 10 == 0 ):
                plt.title("System at z = {}".format(detection_redshift))
            plt.xlabel(r'Velocity ($km\,s^{-1}$)',size=16)
            plt.tick_params(labelsize=15)
            plt.subplots_adjust(wspace = 0, hspace=0)
            plt.axis([-500.,500.,-0.1,1.4])
            ax[axis_index].set_yticks([0,0.5,1])
            if ( n == 1 ) :
                ax[axis_index].set_yticks([])
            plt.figtext(0.05,0.5,'Relative Flux',ha='center',va='center',size=16,rotation=90)
            if ((transition[1]*(1.+detection_redshift)) < (LYA * (1. + zem))):
                plt.plot(velocity[transition[0]],vel_flux[transition[0]],'r-',label=transition[0])
            else:
                plt.plot(velocity[transition[0]],vel_flux[transition[0]],'b-',label=transition[0])
            plt.plot(velocity[transition[0]],np.ones(vel_flux[transition[0]].size),'k--')
            plt.plot(velocity[transition[0]],np.zeros(vel_flux[transition[0]].size),'k--')
            plt.plot(velocity[transition[0]],vel_err[transition[0]],'g-')
            plt.plot([0,0],[0,2],'k--')
            plt.legend(frameon=False,prop={'size':10})
            k += 1
            i += 1
            axis_index += 1
        # End transition library loop

        plt.savefig("{}/extra_trans.eps".format(detection_redshift), format='eps')
        plt.close('all')
    # End detected line redshifts loop
    print "Saved plot files in detection redshift directories."


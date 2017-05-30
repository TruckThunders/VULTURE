import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from fzl_functions import *
from scipy import stats
import timeit
from subprocess import call
np.set_printoptions(threshold=np.nan)

#================================================================
#======================== Initialization ========================
#================================================================

if ( sys.argv[1] == 'clean' ):
    call("rm -r 0.*",shell=True)
    call("rm -r 1.*",shell=True)
    call("rm -r 2.*",shell=True)
    call("rm -r 3.*",shell=True)
    call("rm -r 4.*",shell=True)
    call("rm detections.dat",shell=True)
    sys.exit("========== Cleaned Redshift Directories And Exited ==========")

start = timeit.default_timer()

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

print "Searching redward of Lya, starting at {} A.".format(blue_limit)

emission_limit = 1. + zem - 0.01 # Emission redshift of quasar - 3000 km/s

red_limit = {}
for specie in species:
    red_limit[specie] = []
    # Choose either the reddest_transition or the end of the spectrum
    if ( specie == 'MgII' ):
        if ( (species[specie]['2803'] * emission_limit) < wav[-1]):
            red_limit[specie] = species[specie]['2803'] * emission_limit
        else:
            red_limit[specie] = wav[-1]
    if ( specie == 'CIV' ):
        if ( (species[specie]['1551'] * emission_limit) < wav[-1]):
            red_limit[specie] = species[specie]['1551'] * emission_limit
        else:
            red_limit[specie] = wav[-1]
max_red = max([red_limit['MgII'],red_limit['CIV']])

#===================== Normalize and Clean Up Spectrum ==========================

search_flux = norm_flux[np.where((wav > blue_limit) & (wav < max_red))]
search_wav = wav[np.where((wav > blue_limit) & (wav < max_red))]
search_err = norm_err[np.where((wav > blue_limit) & (wav < max_red))]
search_sig = norm_sig[np.where((wav > blue_limit) & (wav < max_red))]

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

# 6277 - 6318 feature
search_flux[np.where((search_wav > 6277.) & (search_wav < 6318.))] = 1.

# Remove negative errors
search_err[np.where(search_err < 0.)] = 0.

#=================== Make the redshift arrays (z search range) ======================

MgII_lower = (blue_limit / species['MgII']['2796']) - 1.
MgII_upper = (red_limit['MgII'] / species['MgII']['2803']) - 1
CIV_lower = (blue_limit / species['CIV']['1548']) - 1.
CIV_upper = (red_limit['CIV'] / species['CIV']['1551']) - 1.

# Make redshift array with delta_z resolution
redshift = {}
for specie in species:
    redshift[specie] = []
    if ( specie == 'MgII' ):
        lower_lim = MgII_lower
        upper_lim = MgII_upper
        #print "MgII blue_limit = {} red_limit = {}".format(blue_limit,red_limit[specie])
        #print "MgII lower_lim = {} upper_lim = {}.".format(lower_lim,upper_lim)
    if ( specie == 'CIV' ):
        lower_lim = CIV_lower
        upper_lim = CIV_upper
        #print "CIV lower_lim = {} upper_lim = {}.".format(lower_lim,upper_lim)
        #print "CIV blue_limit = {} red_limit = {}".format(blue_limit,red_limit[specie])
    redshift[specie] = np.arange(round(lower_lim,5),round(upper_lim,5),delta_z)

#================ INSERT FAKE LINES AND DETERMINE DETECTION LIMIT =======================

num_fake_lines = 0

#print "--------------------------------"
#print "Inserting fake absorption lines:"

lower_lim = min(MgII_lower,CIV_lower)
upper_lim = max(MgII_upper,CIV_upper)

# Define synthetic line redshift to insert
fake_redshift = np.arange(MgII_lower+0.11,MgII_upper-0.11,((MgII_upper-0.11)-(MgII_lower-0.11))/20.)
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
    if ( detectionlimit > 1 ) :
        fake_redshift = np.delete(fake_redshift,num_fake_lines)
        height = np.delete(height,num_fake_lines)
        continue
    height[num_fake_lines],width[num_fake_lines],fake_flux = insert_lines(fake_flux,search_wav,spectral_res,shifts,'MgII',5.*detectionlimit)
    num_fake_lines += 1

#print "--------------------------------"

#================== END OF FAKE LINE INSERTION ROUTINES ==============================

intermediate = timeit.default_timer()
print " "
print "Run time for prep work: {} mins.".format(round((intermediate-start)/60.,2))
print " "
print "========== Beginning Filtering Routines =========="
print " "
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
for specie in species:
    for z in redshift[specie]:
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

intermediate = timeit.default_timer()
print "Run time to correlation: {} mins.".format(round((intermediate-start)/60.,2))

#==================================================================
#================ Calculate a power spectrum error ================
#==================================================================

print " "
print "========== Calculating Power Spectrum Errors =========="
print " "

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

bottom = 999
top = 0
for specie in species:
    if ( redshift[specie][0] < bottom ):
        bottom = redshift[specie][0]
    if ( redshift[specie][-1] > top ):
        top = redshift[specie][-1]

total_redshift = np.arange(round(bottom,5),round(top,5),delta_z)

# Perform sigma-clipping and generate power-spectrum errors at each redshift
for specie in species:
    for counter,z in enumerate(redshift[specie]):
        flux_index = int(counter + chunk_size)
        for transition in species[specie]:
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

#==================================================================
#============= Find the Lines and Filter the Bad ones =============
#==================================================================

print " "
print "========== Finding Lines and Filtering Bad Ones =========="
print " "

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
        newb = species[specie]['2796'] * (1.+redshift[specie])
        EWwidth = (newb / spectral_res) * (1. / (2.*np.sqrt(2.*np.log(2))))
        EWlimit[specie] = (np.multiply(5.*noise['2796'],EWwidth) * np.sqrt(2.*np.pi)) / m
        # Define final detection redshifts
        criterion[specie] = np.where((sigma['2796'] > 5) & (sigma['2803'] > 3))
        final['2796'] = scaledpower['2796'][criterion[specie][0]]
        final['2803'] = scaledpower['2803'][criterion[specie][0]]
        final_redshift[specie] = redshift[specie][criterion[specie][0]]
    if ( specie == 'CIV' ):
        # Calculate equivalent width detection limit
        EWlimit[specie] = []
        EWlimit[specie] = np.array(EWlimit)
        # Calculate the new b in the gaussian for an unresolved line at the detection limit
        newb = species[specie]['1548'] * (1.+redshift[specie])
        EWwidth = (newb / spectral_res) * (1. / (2.*np.sqrt(2.*np.log(2))))
        EWlimit[specie] = (np.multiply(5.*noise['1548'],EWwidth) * np.sqrt(2.*np.pi)) / m
        # Define final detection redshifts
        criterion[specie] = np.where((sigma['1548'] > 5) & (sigma['1551'] > 3))
        final['1548'] = scaledpower['1548'][criterion[specie][0]]
        final['1551'] = scaledpower['1551'][criterion[specie][0]]
        final_redshift[specie] = redshift[specie][criterion[specie][0]]    
    
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
        line_redshifts[specie] = line_ratio_check(line_redshifts[specie],scaledpower['2796'],scaledpower['2803'],redshift[specie],0.4)
    if ( specie == 'CIV' ):
        line_redshifts[specie] = line_ratio_check(line_redshifts[specie],scaledpower['1548'],scaledpower['1551'],redshift[specie],0.3)

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
        fake_centers = np.where(find_nearest(redshift[specie],lines) == redshift[specie])
        fake_centers = fake_centers[0]
        line_EWlim[specie].append(EWlimit[specie][fake_centers])

    line_EWlim[specie] = np.array(line_EWlim[specie])
    line_EWlim[specie] = line_EWlim[specie].flatten()

    print "--------------------------------------------------------------------------"
    print "Detected {} {} lines:".format(line_redshifts[specie].size,specie)
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
    f.write("Detected {} {} lines:\n".format(line_redshifts[specie].size,specie))
    for item in line_redshifts[specie]:
        f.write("%s  %s\n" % (item,line_EWlim[specie][item]))
    f.write("--------------------------------------------------------------------------\n")
f.close()

# Read in atoms.dat to get atomic data
transition_library = np.genfromtxt('../atoms_UVES.dat',unpack=True,dtype=None)

velwindow = 500. #km/s
velocity = {}
vel_flux = {}
vel_err = {}
ax = ['ax1','ax2','ax3','ax4','ax5','ax6','ax7','ax8','ax9','ax10']
k = 0
i = 1
n = 0

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
                if ( (i % 10 == 0) & (i != 0) ):
                    plt.savefig("{}/vel_plot{}.eps".format(detection_redshift,i/10.),format='eps')
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

        plt.savefig("{}/vel_plot{}.eps".format(detection_redshift,i/10.), format='eps')
        plt.close('all')
    # End detected line redshifts loop
    print "Saved plot files in detection redshift directories."

stop = timeit.default_timer()
print "Total run time: {} mins.".format(round((stop-start)/60.,2))

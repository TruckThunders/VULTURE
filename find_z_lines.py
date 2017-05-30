import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
from fzl_functions import *
import fzl_params
from scipy import stats
import timeit
from subprocess import call
np.set_printoptions(threshold=np.nan)

# Set log options here
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,filemode='w')

# Set log file output
handler = logging.FileHandler('Runtime.log','w')
handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(handler)

try:

    #================================================================
    #======================== Initialization ========================
    #================================================================

    if ( sys.argv[1] == 'clean' ):
        call("rm -r 0.*",shell=True)
        call("rm -r 1.*",shell=True)
        call("rm -r 2.*",shell=True)
        call("rm -r 3.*",shell=True)
        call("rm -r 4.*",shell=True)
        call("rm -r 5.*",shell=True)
        call("rm -r FalsePositive",shell=True)
        call("rm -r UNCOMBINED",shell=True)
        call("rm *.log",shell=True)
        call("rm detections.dat",shell=True)
        call("rm detections_real.dat",shell=True)
        call("rm detections_combined.dat",shell=True)
        sys.exit("========== Cleaned Redshift Directories And Exited ==========\n")

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

    # Read in atoms.dat to get atomic data
    transition_library = fzl_params.transition_library

    # Transition Wavelengths and their power/noise arrays
    species = {}
    power = {}
    noise = {}
    working = {}
    sigma = {}

    # Add transitions. These are the lines you're looking for.
    search_ions = fzl_params.search_ions

    # Loop through atoms.dat and match with search_ions
    # Initialize arrays to search for transitions
    strongest_trans = {}
    second_trans = {}
    for search_trans in search_ions:
        species[search_trans] = {}
        max_oscillator = 0.
        m1 = m2 = float('-inf')
        oldwav1 = oldwav2 = " "
        for transition in transition_library:
            [atom,tranwav] = re.findall('[a-zA-Z]+|\\d+',transition[0])
            current_oscillator = transition[2]
            # Initialize arrays
            if ( atom == search_trans ):
                species[atom][tranwav] = transition[1]
                power[tranwav] = []
                noise[tranwav] = []
                working[tranwav] = []
                sigma[tranwav] = []
                # Find strongest and 2nd strongest transitions (e.g. MgII2796 is strongest MgII transition and MgII2803 is the 2nd)
                if ( transition[2] > m2 ):
                    if ( transition[2] >= m1 ):
                        m1,m2 = transition[2], m1
                        oldwav1,oldwav2 = tranwav, oldwav1
                        strongest_trans[search_trans],second_trans[search_trans] = "{}".format(tranwav),"{}".format(oldwav1)
                    else:
                        m2 = transition[2]
                        oldwav2 = tranwav
                        second_trans[search_trans] = "{}".format(tranwav)
        second_trans[search_trans] = "{}".format(oldwav2)

    # Define Constants
    LIGHTSPEED = 299792.458 # km/s
    LYA = 1215.6 # angstrom

    # Spectral resolution
    spectral_res = fzl_params.spectral_res

    # Redshift resolution
    delta_z = fzl_params.delta_z

    # Read in spectra

    logger.info("========== Reading in Spectrum ==========\n")

    # Ignore leading rows that are only 4 columns long
    ignore = 0
    with open(filename,'rb') as f:
        for line in f:
            check = len(line.split())
            if ( check != 6 ):
                ignore += 1
            else:
                #print "Ignored first {} lines of code because they were dumb.".format(ignore)
                logger.info("Ignored first {} lines of data because they were dumb.\n".format(ignore))
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

    logger.info("========== Initializing Search Arrays ==========\n")
    
    # Search redward of the Lya emission
    if ( (LYA * (1. + zem)) > wav[0]):
        blue_limit = LYA * (1. + zem)
    else:
        blue_limit = wav[0]

    logger.info("Searching redward of Lya, starting at {} A.".format(blue_limit))

    # Emission redshift of quasar - 3000 km/s
    emission_limit = 1. + zem - 0.01

    red_limit = {}
    for specie in species:
        red_limit[specie] = []
        # Choose either the reddest_transition or the end of the spectrum
        reddest = max(species[specie])
        if ( (species[specie][reddest] * emission_limit) < wav[-1]):
            red_limit[specie] = species[specie][reddest] * emission_limit
        else:
            red_limit[specie] = wav[-1]

    max_red = max(red_limit.values())

    #===================== Normalize and Clean Up Spectrum ==========================

    search_flux = norm_flux[np.where((wav > blue_limit) & (wav < max_red))]
    search_wav = wav[np.where((wav > blue_limit) & (wav < max_red))]
    search_err = norm_err[np.where((wav > blue_limit) & (wav < max_red))]
    search_sig = norm_sig[np.where((wav > blue_limit) & (wav < max_red))]

    # Clean cosmic rays/bad high flux values
    search_err[np.where(search_flux > 1.15)] = 0.
    search_sig[np.where(search_flux > 1.15)] = 0.
    search_flux[np.where(search_flux > 1.15)] = 1.

    # Avoid Telluric features, in order from strongest to weakest
    # 9300 - 9630 feature
    search_flux[np.where((search_wav > 9300.) & (search_wav < 9630.))] = 1.

    # 7594 - 7700 feature
    search_flux[np.where((search_wav > 7594.) & (search_wav < 7700.))] = 1.

    # 6868 - 6932 feature
    search_flux[np.where((search_wav > 6868.) & (search_wav < 6932.))] = 1.

    # 6277 - 6318 feature
    search_flux[np.where((search_wav > 6277.) & (search_wav < 6318.))] = 1.

    # Remove negative errors
    search_sig[np.where(search_err < 0.)] = 0.
    search_flux[np.where(search_err < 0.)] = 1.
    search_err[np.where(search_err < 0.)] = 0.

    # Remove negative flux
    search_flux[np.where(search_flux < 0.)] = 0.

    #=================== Make the redshift arrays (z search range) ======================

    upper = {}
    lower = {}

    # Assign upper and lower redshift limits for each ion
    for specie in species:
        upper[specie] = (red_limit[specie] / species[specie][max(species[specie])]) - 1.
        lower[specie] = (blue_limit / species[specie][min(species[specie])]) - 1.

    # Check to see if ions are covered in the redshift path length - if not, remove them from the search arrays.
    killflag = []
    looplist = list(species)
    for specie in looplist:
        if ( lower[specie] > upper[specie] ):
            logger.warning("***** No wavelength coverage for {} - Not searching for {} *****".format(specie,specie))
            for transition in species[specie]:
                del working[transition]
                del noise[transition]
            del species[specie]
            killflag.append(0)
        elif ( (upper[specie] - lower[specie]) < 0.02 ):
            logger.warning("***** Too little wavelength coverage for {} - Not searching for {} *****".format(specie,specie))
            for transition in species[specie]:
                del working[transition]
                del noise[transition]
            del species[specie]
            killflag.append(0)
        else:
            killflag.append(1)

    # If no wavelength coverage for any ions, kill the program and report no detections
    if ( np.count_nonzero(killflag) == 0 ):
        logger.info("========== No wavelength coverage for {} or {}, exiting with no detections ==========\n".format(search_ions[0],search_ions[1]))
        sys.exit("========== No wavelength coverage for {} or {}, exiting with no detections ==========\n".format(search_ions[0],search_ions[1]))

    # Make redshift array with delta_z resolution
    redshift = {}
    for specie in species:
        redshift[specie] = []
        logger.info("{} blue limit = {} red limit = {}".format(specie,blue_limit,red_limit[specie]))
        logger.info("{} min_z = {} max_z = {}.".format(specie,lower[specie],upper[specie]))
        redshift[specie] = np.arange(round(lower[specie],5),round(upper[specie],5),delta_z)

    # Define number of transitions you're searching for
    num_trans = len(species)

    #================ INSERT FAKE LINES AND DETERMINE DETECTION LIMIT =======================

    num_fake_lines = 0

    lower_lim = min(lower.values()) + 0.01
    upper_lim = max(upper.values()) - 0.01

    # Define synthetic line redshift to insert
    fake_redshift = np.arange(lower_lim,upper_lim,((upper_lim-lower_lim)/40.))
    height=np.zeros(len(fake_redshift))
    width=np.zeros(len(fake_redshift))

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
        truthcounter = 0
        for specie in species:
            # Make sure we're using a rest wavelength within the search window
            if ( (shifts > lower[specie]) & (shifts < upper[specie]) ):
                center = np.where(eq_wav == find_nearest(eq_wav,min(species[specie].values())*(1.+shifts)))
                center = int(center[0])
                # Make sure center is 250 units away from the beginning or end of spectrum
                if ( (center < 251) or ((len(eq_fivesigma) - center) < 251)):
                    fake_redshift = np.delete(fake_redshift,num_fake_lines)
                    height = np.delete(height,num_fake_lines)
                    width = np.delete(width,num_fake_lines)
                    break
                detectionlimit = np.fabs(np.mean(eq_fivesigma[center-250:center+250]))
                # Do not examine regions where the detection limit is outrageously bad
                if ( 5.*detectionlimit > 0.5 ):
                    fake_redshift = np.delete(fake_redshift,num_fake_lines)
                    height = np.delete(height,num_fake_lines)
                    width = np.delete(width,num_fake_lines)
                    break
                # Calculate height and width of gaussian representing 5 sigma detection of unresolved line
                height[num_fake_lines],width[num_fake_lines],fake_flux = insert_lines(fake_flux,search_wav,spectral_res,shifts,5.*detectionlimit,min(species[specie].values()))
                num_fake_lines += 1
                break
            else:
                truthcounter += 1
                if ( truthcounter == num_trans):
                    logger.info("Fake line no coverage at z = {}.".format(shifts))
                    fake_redshift = np.delete(fake_redshift,num_fake_lines)
                    height = np.delete(height,num_fake_lines)
                    width = np.delete(width,num_fake_lines)
                continue

    #================== END OF FAKE LINE INSERTION ROUTINES ==============================

    intermediate = timeit.default_timer()
    logger.info("Run time for prep work: {} mins.\n".format(round((intermediate-start)/60.,2)))
    logger.info("========== Beginning Filtering Routines ==========\n")

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

    logger.info("Min Value Established: {}".format(min_value))

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

    logger.info("Max Value Established: {}".format(max_value))

    fake_power = []
    # Loop over fake lines to get fake values for EW
    for fake_z in fake_redshift:
        # Make sure we're using a rest wavelength within the search range
        for specie in species:
            if ( (fake_z > lower[specie]) & (fake_z < upper[specie]) ):
                rest_wavelength = min(species[specie].values())
                break
            else:
                continue
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
    logger.info("Run time to correlation: {} mins.\n".format(round((intermediate-start)/60.,2)))

    #==================================================================
    #================ Calculate a power spectrum error ================
    #==================================================================

    logger.info("========== Calculating Power Spectrum Errors ==========\n")

    # Chunk the power spectrum
    chunk_size = (0.03 / delta_z)
    #chunk_size = 1000.
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
        #working[transition][index2:(len(scaledpower[transition])+index2)] = scaledpower[transition]
        working[transition][index2:(len(scaledpower[transition])+index2-1)] = scaledpower[transition]
        working[transition][(len(scaledpower[transition])+index4):(len(scaledpower[transition])+index3)] = scaledpower[transition][(len(scaledpower[transition])-index2):-1]

        # In the error array, erase the null chunks by filling them with mirrored data
        zerovals = np.where((working[transition] < 1.e-5) & (working[transition] > -1.e-5))
        zerovals = zerovals[0]
        gaps = group_consecutives(zerovals,1)

        # Skip over any chunk less than 5 points in length for transition
        for spots in gaps:
            if len(spots) < 5:
                continue
            working[transition][spots] = np.copy(working[transition][np.subtract(spots,len(spots)).astype(int)])

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
    logger.info("Params for detection limit fit: m = {}   b = {}.\n".format(m,b))

    #==================================================================
    #============= Find the Lines and Filter the Bad ones =============
    #==================================================================

    logger.info("========== Finding Lines and Filtering Bad Ones ==========\n")

    criterion = {}
    final = {}
    final_redshift = {}
    line_redshifts = {}
    line_EWlim = {}
    line_sigma = {}
    grouped_redshifts = {}
    EWlimit = {}

    for specie in species:
        # Define criterion for detecting an absorption line system 
        EWlimit[specie] = []
        EWlimit[specie] = np.array(EWlimit)
        # Calculate the new b in the gaussian for an unresolved line at the detection limit
        newb = species[specie][strongest_trans[specie]] * (1.+redshift[specie])
        EWwidth = (newb / spectral_res) * (1. / (2.*np.sqrt(2.*np.log(2))))
        EWlimit[specie] = (np.multiply(5.*noise[strongest_trans[specie]],EWwidth) * np.sqrt(2.*np.pi)) / m
        # Define final detection redshifts
        criterion[specie] = np.where((sigma[strongest_trans[specie]] > 5) & (sigma[second_trans[specie]] > 3))
        final[strongest_trans[specie]] = scaledpower[strongest_trans[specie]][criterion[specie][0]]
        final[second_trans[specie]] = scaledpower[second_trans[specie]][criterion[specie][0]]
        final_redshift[specie] = redshift[specie][criterion[specie][0]]

        # Group the adjacent redshifts to find the detection centers
        grouped_redshifts[specie] = group_consecutives(final_redshift[specie],delta_z)

        # Find line centers
        line_redshifts[specie] = []
        for lines in grouped_redshifts[specie]:
            center = np.mean(lines)
            line_redshifts[specie].append(center)

        logger.info("========== {} Lines detected. Culling bad lines ==========\n".format(specie))

        # THROW OUT BAD DETECTIONS AND FLAG WEAK DETECTIONS IF QUESTIONABLE
        # Bad chunk of spectrum check
        # Thow out detections where the S/N < 2
        counter = 0
        for detection in line_redshifts[specie]:
            rest_wavelength = min(species[specie].values())
            central = find_nearest(search_wav,rest_wavelength*(1.+detection))
            central = np.where(search_wav == central)
            central = central[0]
            # Delete lines if sigma > continuum
            if (search_sig[central] > 1.):
                line_redshifts[specie] = np.delete(line_redshifts[specie],counter)
                logger.info("Deleted {} line (Noise > Signal): {}".format(specie,detection))
                counter -= 1
            # Delete lines if the error is more than 1/2 the continuum
            elif ((1./search_sig[central]) < 2.):
                line_redshifts[specie] = np.delete(line_redshifts[specie],counter)
                logger.info("Deleted {} line (S/N < 2): {}".format(specie,detection))
                counter -= 1
            counter += 1                
        
        # Line Ratio Check
        dub_ratio = fzl_params.dub_ratio
        line_redshifts[specie] = line_ratio_check(specie,line_redshifts[specie],noise[strongest_trans[specie]],scaledpower[strongest_trans[specie]],scaledpower[second_trans[specie]],redshift[specie],dub_ratio)

        # Calculate EW Detection limit for each real detection
        line_EWlim[specie] = []
        line_sigma[specie] = []
        for lines in line_redshifts[specie]:
            fake_centers = np.where(find_nearest(redshift[specie],lines) == redshift[specie])
            fake_centers = fake_centers[0]
            line_EWlim[specie].append(EWlimit[specie][fake_centers])
            line_sigma[specie].append(sigma[strongest_trans[specie]][fake_centers])

        # Flatten arrays for easy output
        line_EWlim[specie] = np.array(line_EWlim[specie])
        line_EWlim[specie] = line_EWlim[specie].flatten()
        line_sigma[specie] = np.array(line_sigma[specie])
        line_sigma[specie] = line_sigma[specie].flatten()
        line_redshifts[specie] = np.array(line_redshifts[specie])
        line_redshifts[specie] = line_redshifts[specie].flatten()

        # Report the detected lines to log file and terminal
        logger.info("--------------------------------------------------------------------------")
        logger.info("Detected {} {} lines:".format(line_redshifts[specie].size,specie))
        logger.info(line_redshifts[specie])
        logger.info("Detection limit for each detection (Angstrom):")
        logger.info(line_EWlim[specie])
        logger.info("Line sigma:")
        logger.info(line_sigma[specie])
        logger.info("--------------------------------------------------------------------------\n")

    #==================================================================
    #========================== File Output ===========================
    #==================================================================

    # Output ewlimit spectrum
    for specie in species:
        ewfile = open("{}_ewlim.data".format(specie),"w")
        for i,z in enumerate(redshift[specie]):
            ewfile.write("{}   {}\n".format(z, EWlimit[specie][i]))

    # Output line detections
    f = open('detections.dat','w')
    for specie in species:
        f.write("Detected {} {} lines:\n".format(line_redshifts[specie].size,specie))
        f.write("--------------------------------------------------------------------------\n")
        f.write("Redshift         EW_limit         Power_Rating\n")
        for num,item in enumerate(line_redshifts[specie]):
            f.write("%10s   %6s   %4s\n" % (item,line_EWlim[specie][num],line_sigma[specie][num]))
        f.write("--------------------------------------------------------------------------\n")
    f.close()

    # Plotting routines. Can edit velocity window in params file
    # Outputs plots based upon your atoms.dat input file so you can see all the transitions.
    # Red = in the Lya forest
    # Blue = redward of Lya forest
    velwindow = fzl_params.velwindow
    velocity = {}
    vel_flux = {}
    vel_err = {}
    ax = ['ax1','ax2','ax3','ax4','ax5','ax6','ax7','ax8','ax9','ax10']
    k = 0
    i = 1
    n = 0

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
                        plt.savefig("{}/vel_plot{}.eps".format(detection_redshift,int(i/10.)),format='eps')
                        plt.close('all')
                        axis_index = 0
                        n = 0
                        plt.figure(figsize=(8.5,11))
                    k = 0
                ax[axis_index] = plt.subplot2grid((5,2),(k,n))
                if ( i % 10 == 0 ):
                    plt.title("System at z = {}".format(detection_redshift))
                plt.xlabel(r'Velocity ($km\,s^{-1}$)',size=16)
                plt.tick_params(labelsize=15)
                plt.subplots_adjust(wspace = 0, hspace=0)
                plt.axis([-velwindow,velwindow,-0.1,1.4])
                ax[axis_index].set_yticks([0,0.5,1])
                if ( n == 1 ) :
                    ax[axis_index].set_yticks([])
                plt.figtext(0.05,0.5,'Relative Flux',ha='center',va='center',size=16,rotation=90)
                if ((transition[1]*(1.+detection_redshift)) < (LYA * (1. + zem))):
                    plt.step(velocity[transition[0]],vel_flux[transition[0]],'r-',label=transition[0])
                else:
                    plt.step(velocity[transition[0]],vel_flux[transition[0]],'b-',label=transition[0])
                plt.plot(velocity[transition[0]],np.ones(vel_flux[transition[0]].size),'k--')
                plt.plot(velocity[transition[0]],np.zeros(vel_flux[transition[0]].size),'k--')
                plt.step(velocity[transition[0]],vel_err[transition[0]],'g-')
                plt.plot([0,0],[0,2],'k--')
                plt.legend(frameon=False,prop={'size':10})
                k += 1
                i += 1
                axis_index += 1
            # End transition library loop

            plt.savefig("{}/vel_plot{}.eps".format(detection_redshift,int(i/10.)), format='eps')
            plt.close('all')
        # End detected line redshifts loop
        logger.info("Saved plot files in detection redshift directories.")

    stop = timeit.default_timer()
    logger.info("Total run time: {} mins.".format(round((stop-start)/60.,2)))

except Exception, e:
    logger.error('Critical Error in {}, Program Exited.'.format(str(sys.argv[1])),exc_info=True)


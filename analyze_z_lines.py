# PYTHON CODE TO REBIN VELOCITY PLOT SPECTRA,
# OUTPUT VELOCITY DATA FILES, AND
# OUTPUT ***.EWREG FILES FOR SYSANAL

import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from fzl_functions import *
from scipy import stats
from subprocess import call
np.set_printoptions(threshold=np.nan)

# Set log options here
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,filemode='w')

# Set log file output
handler = logging.FileHandler('Analysis.log','w')
handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(handler)

try:

    # Make sure visual inspection has happened
    if ( os.path.isfile("detections_real.dat") == False ):
        sys.exit("No detections_real.dat. Please visually inspect directory using check_z_lines.py.")

    # Read in zem.dat
    zem = np.loadtxt("zem.dat")

    readfile = "detections_real.dat"

    # Read in command line args
    filename = str(sys.argv[1])

    # Create lists
    wav = []
    flux = []
    sig = []
    err = []
    cont = []

    # Define Constants
    LIGHTSPEED = 299792.458 # km/s
    LYA = 1215.6 # angstrom

    # Convert lists to numpy arrays
    wav = np.array(wav)
    flux = np.array(flux)
    sig = np.array(sig)
    err = np.array(err)
    cont = np.array(cont)

    # Array definitions
    redshift = {}
    ewlim = {}
    power_rating = {}

    redshift['CIV'] = []
    ewlim['CIV'] = []
    power_rating['CIV'] = []

    redshift['MgII'] = []
    ewlim['MgII'] = []
    power_rating['MgII'] = []

    delimiter = "CIV"

    logger.info("========== Reading in Spectrum ==========\n")

    # Ignore leading rows that are only 4 columns long
    ignore = 0
    with open(filename,'rb') as f:
        for line in f:
            check = len(line.split())
            if ( check == 5 ):
                ignore += 1
            else:
                #print "Ignored first {} lines of code because they were dumb.".format(ignore)
                logger.info("Ignored first {} lines of data because they were dumb.\n".format(ignore))
                break

            # end else
        # end for
    #end with

    wav,dum,flux,sig,err,cont = np.loadtxt(filename,unpack=True,skiprows=ignore)

    norm_flux = flux / cont
    norm_sig = sig / cont
    norm_err = err / cont

    # Read in detections.dat and store the information
    logger.info("========== Reading in detections file ==========\n")
    header = 0
    with open(readfile,'rb') as f:
        alltext = f.readlines()
        for row in alltext:
            words = row.split()
            check = len(words)
            if ( check == 1 ):
                logger.info("Next because: {}".format(row))
            elif ( check == 4 ):
                if ( words[2] == "MgII" ):
                    delimiter = "MgII"
            elif ( words[0] == "Redshift" ):
                logger.info("Next because: {}".format(row))
            else:
                logger.info("Found good row.")
                redshift[delimiter].append(words[0])
                ewlim[delimiter].append(float(words[1]))
                power_rating[delimiter].append(float(words[2]))
                
    logger.info("========== Detections read in. Analyzing. ==========\n")

    # Read in atoms.dat to get atomic data
    transition_library = np.genfromtxt('/home/achernar-data/mathes/UVES_SQUAD/atoms_UVES.dat',unpack=True,dtype=None)

    velwindow = 500. #km/s
    velocity = {}
    vel_flux = {}
    vel_err = {}
    ax = ['ax1','ax2','ax3','ax4','ax5','ax6','ax7','ax8','ax9','ax10']
    k = 0
    i = 1
    n = 0

    # Output velocity files for detections
    for specie in redshift:
        sysvel = [0] * len(redshift[specie])
        # Detected line redshifts
        for systems in redshift[specie]:
            systems = float(systems)
            i = 0
            n = 0
            axis_index = 0
            plt.figure(figsize=(8.5,11))
            call(["cp", "/achernar-data/mathes/UVES_SQUAD/sysanal.inp", "{}/.".format(systems)])
            # Write ions.table and zabs.fzl for each system redshift directory
            zabsfile = open("{}/zabs.fzl".format(systems),"w")
            if ( specie == "MgII" ):
                zabsfile.write("  %.7f   %.4f  MgII2796\n" % (systems,(1.+systems)*2796.352))
            if ( specie == "CIV" ):
                zabsfile.write("  %.7f   %.4f  CIV1548\n" % (systems,(1.+systems)*1548.195))
            # Close zabs.fzl
            zabsfile.close()
            # Loop over transitions in atoms_UVES.dat and output velocity .data files, .ewreg files, and velocity plots
            for transition in transition_library:
                # Calculate min/max wavelength regions and create the arrays for the velocity data
                lambda_r = transition[1] * (1. + systems)
                lambda_min = (-(velwindow/LIGHTSPEED) + 1)*lambda_r
                lambda_max = ((velwindow/LIGHTSPEED) + 1)*lambda_r
                index_min = np.where(wav == find_nearest(wav,lambda_min))
                index_max = np.where(wav == find_nearest(wav,lambda_max))
                index_mid = np.where(wav == find_nearest(wav,lambda_r))
                index_min = index_min[0]
                index_max = index_max[0] + 1
                index_mid = index_mid[0]
                velocity[transition[0]] = np.zeros(index_max - index_min)
                vel_flux[transition[0]] = np.zeros(index_max - index_min)
                vel_err[transition[0]] = np.zeros(index_max - index_min)
                # Write velocity .data file for systems/transition[0]
                velfile = open("{}/{}.data".format(systems,transition[0]),"w")
                j = 0
                for indices in range(index_min,index_max):
                    if ( (index_max - index_min) < 3 ):
                        break
                    velocity[transition[0]][j] = (LIGHTSPEED * (wav[indices] * wav[indices] - lambda_r * lambda_r)) / (wav[indices] * wav[indices] + lambda_r * lambda_r)
                    vel_flux[transition[0]][j] = norm_flux[indices]
                    vel_err[transition[0]][j] = norm_sig[indices]
                    velfile.write("%10s   %10s   %10s  %10s  %10s  1.000\n" % (wav[indices],velocity[transition[0]][j],vel_flux[transition[0]][j],vel_err[transition[0]][j],vel_err[transition[0]][j]))
                    j += 1
                # End Velocity for loop
                velfile.close()
                logger.info("Wrote velocity file for {}/{}.".format(systems,transition[0]))
                call(["cp", "{}/{}.data".format(systems,transition[0]), "{}/{}".format(systems,transition[0])])
                
                # If secondary transition, don't bother with the ball rolling
                # Instead, use MgII2796 or CIV1548 as the velocity bounds for the 2803/1551 detections
                if ( (transition[0] == "MgII2803") or (transition[0] == "CIV1551")):
                    position = index_mid - roll_left
                    detection_min = wav[position]
                    vel_min = (LIGHTSPEED * (wav[position] * wav[position] - lambda_r * lambda_r)) / (wav[position] * wav[position] + lambda_r * lambda_r)
                    position = index_mid + roll_right
                    vel_max = (LIGHTSPEED * (wav[position] * wav[position] - lambda_r * lambda_r)) / (wav[position] * wav[position] + lambda_r * lambda_r)
                    detection_max = wav[position]
                else:
                    # Roll ball around mid point to trace down to 1-sigma detection limit for sysanal.ewreg output
                    # Roll left
                    roll_left = 1
                    while True:
                        position = index_mid - roll_left
                        if (position < 0):
                            vel_min = [0.]
                            detection_min = wav[index_mid]
                            roll_left = 0
                            break
                        elif norm_flux[position] > (1. - norm_sig[position]):
                            detection_min = wav[position]
                            vel_min = (LIGHTSPEED * (wav[position] * wav[position] - lambda_r * lambda_r)) / (wav[position] * wav[position] + lambda_r * lambda_r)
                            break
                        roll_left += 1
                    # Roll right
                    roll_right = 1
                    while True:
                        position = index_mid + roll_right
                        if (position >= norm_flux.size):
                            vel_max = [0.]
                            detection_max = wav[index_mid]
                            roll_right = 0
                            break
                        elif norm_flux[position] > (1. - norm_sig[position]):
                            detection_max = wav[position]
                            vel_max = (LIGHTSPEED * (wav[position] * wav[position] - lambda_r * lambda_r)) / (wav[position] * wav[position] + lambda_r * lambda_r)
                            break
                        roll_right += 1
                # Plotting time! Lots of magic here. Makes 2 column 5 row velocity plots for spectra.
                if ( i % 5 == 0 ):
                    if ( i != 0 ):
                        n = 1
                    if ( (i % 10 == 0) & (i != 0) ):
                        plt.savefig("{}/vel_plot{}.eps".format(systems,int(i/10.)),format='eps')
                        logger.info("Wrote {}/vel_plot{}.eps".format(systems,int(i/10.)))
                        plt.close('all')
                        axis_index = 0
                        n = 0
                        plt.figure(figsize=(8.5,11))
                    k = 0
                ax[axis_index] = plt.subplot2grid((5,2),(k,n))
                if ( i % 10 == 0 ):
                    plt.title("System at z = {}".format(systems))
                plt.xlabel(r'Velocity ($km\,s^{-1}$)',size=16)
                plt.tick_params(labelsize=15)
                plt.subplots_adjust(wspace = 0, hspace=0)
                plt.axis([-500.,500.,-0.1,1.4])
                ax[axis_index].set_yticks([0,0.5,1])
                if ( n == 1 ) :
                    ax[axis_index].set_yticks([])
                plt.figtext(0.05,0.5,'Relative Flux',ha='center',va='center',size=16,rotation=90)
                if ((transition[1]*(1.+systems)) < (LYA * (1. + zem))):
                    plt.step(velocity[transition[0]],vel_flux[transition[0]],'r-',label=transition[0],where="mid")
                else:
                    plt.step(velocity[transition[0]],vel_flux[transition[0]],'b-',label=transition[0],where="mid")
                plt.plot(velocity[transition[0]],np.ones(vel_flux[transition[0]].size),'k--')
                plt.plot(velocity[transition[0]],np.zeros(vel_flux[transition[0]].size),'k--')
                plt.step(velocity[transition[0]],vel_err[transition[0]],'g-',where="mid")
                plt.plot([0,0],[0,2],'k--')
                plt.legend(frameon=False,prop={'size':10})
                k += 1
                i += 1
                axis_index += 1

                # Write .ewreg file for systems/transition[0]
                sysfile = open("{}/{}.ewreg".format(systems,transition[0]),"w")
                sysfile.write("  Reg  Sub  ypos       Lambda-     Lambda+     Vel-      Vel+     Lambda0   Detect\n")
                sysfile.write("    1    1  1.25  %12.4f%12.4f%12.4f%12.4f%12.4f   T" % (detection_min[0],detection_max[0],vel_min[0],vel_max[0],transition[1]))
                sysfile.close()

                # Loop from CIV1548 middle pixel in absorption to find optical depth center
                if ( (specie == 'CIV') and (transition[0] == "CIV1548") ):
                    #logger.info("Entered CIV1548 optical depth loop.")
                    move = 0
                    left = False
                    right = False
                    while True:
                        tauleft = 0.
                        tauright = 0.
                        for index in range(index_min,index_mid):
                            if ( norm_flux[index] > 0. ):
                                tauleft += -np.log(norm_flux[index])
                            else:
                                tauleft += -np.log(norm_sig[index])
                        for index in range(index_mid,index_max):
                            if ( norm_flux[index] > 0. ):
                                tauright += -np.log(norm_flux[index])
                            else:
                                tauright += -np.log(norm_sig[index])
                        #logger.info("tauright = {}  tauleft = {}".format(tauright,tauleft))
                        if ( tauleft > tauright ):
                            #logger.info("Entered tauleft > tauright")
                            index_mid -= 1
                            left = True
                            move = -1
                        else: 
                            #logger.info("Entered tauleft < tauright")
                            index_mid += 1
                            right = True
                            move = 1
                        if ( (right == True) and (left == True) ):
                            if ( move == -1 ):
                                index_mid += 1
                            elif ( move == 1 ):
                                index_mid -= 1
                            else: 
                                logger.error('Left and Right = True, but move = 0. Error.',exec_info=True)
                            #logger.info("Found Tau-weighted mean redshift at {}, z_tau = {}.".format(wav[index_mid],(wav[index_mid]/1548.195) - 1.))
                            break

                    # Find the optical depth weighted mean and output zabs.dat file
                    z_tau = (wav[index_mid]/1548.195) - 1.

                    syszabsfile = open("{}/zabs.dat".format(systems),"w")
                    syszabsfile.write("  %.7f   %.4f  CIV1548\n" % (z_tau,(1. + z_tau)*1548.195))
                    syszabsfile.close()

                # Loop from MgII2796 middle pixel in absorption to find optical depth center
                if ( (specie == 'MgII') and (transition[0] == "MgII2796") ):
                    #logger.info("Entered MgII2796 optical depth loop.")
                    move = 0
                    left = False
                    right = False
                    while True:
                        tauleft = 0.
                        tauright = 0.
                        for index in range(index_min,index_mid):
                            if ( norm_flux[index] > 0. ):
                                tauleft += -np.log(norm_flux[index])
                            else:
                                tauleft += -np.log(norm_sig[index])
                        for index in range(index_mid,index_max):
                            if ( norm_flux[index] > 0. ):
                                tauright += -np.log(norm_flux[index])
                            else:
                                tauright += -np.log(norm_sig[index])
                        if ( tauleft > tauright ):
                            #logger.info("Entered tauleft > tauright")
                            index_mid -= 1
                            left = True
                            move = -1
                        else: 
                            #logger.info("Entered tauright > tauleft")
                            index_mid += 1
                            right = True
                            move = 1
                        if ( (right == True) and (left == True) ):
                            if ( move == -1 ):
                                index_mid += 1
                            elif ( move == 1 ):
                                index_mid -= 1
                            else: 
                                logger.error('Left and Right = True, but move = 0. Error.',exec_info=True)
                            #logger.info("Found Tau-weighted mean redshift at {}, z_tau = {}.".format(wav[index_mid],(wav[index_mid]/2796.352) - 1.))
                            break

                    # Find the optical depth weighted mean and output zabs.dat file
                    z_tau = (wav[index_mid]/2796.352) - 1.

                    syszabsfile = open("{}/zabs.dat".format(systems),"w")
                    syszabsfile.write("  %.7f   %.4f  MgII2796\n" % (z_tau,(1.+z_tau)*2796.352))
                    syszabsfile.close()

            # End transition library loop
            
            ionsfile = open("{}/ions.table".format(systems),"w")

            if (specie == "CIV"):
                ionsfile.write("CIV1548  1  1.0  0.0  5.0\n")
                ionsfile.write("CIV1551  1  1.0  0.0  3.0\n")
                
            if (specie == "MgII"):
                ionsfile.write("MgII2796  1  1.0  0.0  5.0\n")
                ionsfile.write("MgII2803  1  1.0  0.0  3.0\n")

            ionsfile.close()

            plt.savefig("{}/vel_plot{}.eps".format(systems,int(i/10.)), format='eps')
            logger.info("Wrote {}/vel_plot{}.eps".format(systems,int(i/10.)))
            plt.close('all')

        # End detected line redshifts loop
        logger.info("Saved {} velocity data files in detection redshift directories.".format(specie))

except Exception, e:
    logger.error('Critical Error. Program Exited.',exc_info=True)


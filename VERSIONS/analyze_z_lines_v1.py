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
handler = logging.FileHandler('Runtime.log','w')
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

            # 1.) Trace line down to 1 sigma to get width
            # 2.) Output velocity plots of +/- 500 km/s

            # Lya.data example
            # wav           vel            flux             err             err          norm_cont (1)
            #1392.8040  -1999.0869      6.354745E-01    3.998853E-02    3.998853E-02    1.000000E+00
            #1392.8140  -1996.9488      6.585073E-01    4.064885E-02    4.064885E-02    1.000000E+00
            #1392.8240  -1994.8108      6.443548E-01    4.055932E-02    4.055932E-02    1.000000E+00
            #1392.8339  -1992.6940      7.399658E-01    4.219596E-02    4.219596E-02    1.000000E+00

            # 3.) Output sysanal files - ION.ewreg

            # Lya.ewreg example
            #  Reg  Sub  ypos       Lambda-     Lambda+       Vel-        Vel+      Lambda0   Detect
            #    1    1  1.25     1399.4641   1400.1421   -575.1008   -430.1386   1215.6701   T
            #    2    1  1.25     1400.3515   1400.9796   -385.3671   -251.0739   1215.6701   T
            #    3    1  1.25     1401.1690   1402.3056   -210.5786     32.4362   1215.6701   T
            #    4    1  1.25     1402.5050   1403.1431     75.0696    211.5008   1215.6701   T

    # Output velocity files for detections
    for specie in redshift:
        # Detected line redshifts
        for systems in redshift[specie]:
            systems = float(systems)
            i = 0
            n = 0
            axis_index = 0
            plt.figure(figsize=(8.5,11))
            for transition in transition_library:
                velfile = open("{}.data".format(transition),"w")
                sysfile = open("{}.ewreg".format(transition),"w")
                lambda_r = transition[1] * (1. + systems)
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
                    velfile.write("%10s   %10s   %10s  %10s  %10s  1.000\n" % (wav[index_min + indices],velocity[transition[0]][j],vel_flux[transition[0]][j],vel_err[transition[0]][j],vel_err[transition[0]][j]))
                # End Velocity for loop
                if ( i % 5 == 0 ):
                    if ( i != 0 ):
                        n = 1
                    if ( (i % 10 == 0) & (i != 0) ):
                        plt.savefig("{}/vel_plot{}.eps".format(systems,int(i/10.)),format='eps')
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
                plt.step(velocity[transition[0]],np.ones(vel_flux[transition[0]].size),'k--',where="mid")
                plt.step(velocity[transition[0]],np.zeros(vel_flux[transition[0]].size),'k--',where="mid")
                plt.step(velocity[transition[0]],vel_err[transition[0]],'g-',where="mid")
                plt.step([0,0],[0,2],'k--',where="mid")
                plt.legend(frameon=False,prop={'size':10})
                k += 1
                i += 1
                axis_index += 1
                sysfile.write("  1   1   1.25   %10s  %10s  %10s  %10s  %10s  T\n" % (lambda_min,lambda_max,vel[index_min],vel[index_max],lambda_r))
            velfile.close()
            sysfile.close()
            # End transition library loop

        # End detected line redshifts loop
        logger.info("Saved velocity data files in detection redshift directories.")



except Exception, e:
    logger.error('Critical Error in {}, Program Exited.'.format(str(sys.argv[1])),exc_info=True)


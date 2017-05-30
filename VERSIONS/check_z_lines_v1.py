import string
import numpy as np
import sys
import os
import timeit
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

# Function recursive function to ask user for input and return whether the detection is real (True) or fake (False)
def questionlogic(answer):

    finalanswer = True

    # Yes, good detection
    if ( answer == "y" ):
        logger.info("Good detection. Keeping things.")

    # More, look at Al/Mn/Dust transitions
    elif ( answer == "m" ):
        call("gv {}/vel_plot1.eps & gv {}/vel_plot2.eps".format(redshift[specie][count],redshift[specie][count]),shell=True)
        answer = raw_input("Is this a good detection? (y/n/m/a/q/r)")
        answer.lower()
        finalanswer = questionlogic(answer)

    # All, look at all 3 vel plots together
    elif ( answer == "a" ):
        call("gv {}/vel_plot1.eps & gv {}/vel_plot2.eps & gv {}/vel_plot3.0.eps".format(redshift[specie][count],redshift[specie][count],redshift[specie][count]),shell=True)
        answer = raw_input("Is this a good detection? (y/n/m/a/q/r)")
        answer.lower()
        finalanswer = questionlogic(answer)

    # No, false positive
    elif ( answer == "n" ):
        logger.info("Bad detection. {} moved to FalsePositive directory.".format(redshift[specie][count]))
        finalanswer = False

    # Quit program
    elif ( answer == "q"):
        # Open file to write
        outfile = open(writefile,'w')

        logger.info("Writing detections_real.dat")
        for line in redshift:
            outfile.write("Detected {} {} lines:\n".format(len(redshift[line]),line))
            outfile.write("--------------------------------------------------------------------------\n")
            outfile.write("Redshift         EW_limit         Power_Rating\n")
            for num,item in enumerate(redshift[line]):
                outfile.write("%10s   %6s   %4s\n" % (item,ewlim[line][num],power_rating[line][num]))
            outfile.write("--------------------------------------------------------------------------\n")
        outfile.close()

        sys.exit("Quit command detected. Exiting.")

    # Reset directory
    elif ( answer == "r"):
        call(["mv","detections_real.dat","detections_real.old"])
        call("mv FalsePositive/* .",shell=True)
        sys.exit("Quitting and resetting.")

    # Display help menu
    elif ( answer == "h"):
        print " "
        print "y = yes, good detection. Moving on"
        print "n = no, not a good detection. Remove from detections_real.dat and move to FalsePositive"
        print "m = more, show Al/Mn/Dust transitions"
        print "a = all, show all available transitions"
        print "q = quit where you are and write detections_real.dat"
        print "r = reset the directory"
        print " "
        answer = raw_input("Is this a good detection? (y/n/m/a/q/r)")
        answer.lower()
        finalanswer = questionlogic(answer)

    # Bad input, try again
    else:
        logger.info("Bad input detected. Try again, this time with feeling.")
        answer = raw_input("Is this a good detection? (y/n/m/a/q/r)")
        answer.lower()
        finalanswer = questionlogic(answer)
    
    return finalanswer

# Do stuff
try: 

    # If no filtering has been done before, make the FalsePositive directory
    if ( os.path.exists("FalsePositive") == False ):
        call(["mkdir","FalsePositive"])

    readfile = "detections.dat"

    # If filtering has been done before, save the old file as *.old
    if ( os.path.isfile("detections_real.dat") == True ):
        answer = raw_input("Check has been run before. Start over? (y/n)")
        answer.lower()
        if ( answer == "y" ):
            logger.info("Starting over from scratch. Reading in detections.dat")
            call(["mv","detections_real.dat","detections_real.old"])
            call("mv FalsePositive/* .",shell=True)
            readfile = 'detections.dat'
        elif ( answer == "n" ):
            logger.info("Checking twice, using pre-checked sample. Reading in detections_old.dat")
            call(["mv","detections_real.dat","detections_real.old"])
            readfile = 'detections_real.old'
        else:
            logger.info("No valid input detected. Defaulting to starting over from scratch.")
            call(["mv","detections_real.dat","detections_real.old"])
            call("mv FalsePositive/* .",shell=True)
            readfile = 'detections.dat'

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

    # Output file name
    writefile = "detections_real.dat"

    delimiter = "CIV"

    # Read in detections.dat and store the information
    logger.info("========== Reading in file ==========\n")
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

    for specie in redshift:
        count = 0
        for i,systems in enumerate(redshift[specie]):
            finalanswer = True
            if ( power_rating[specie][count] < 20. ):
                # Visually inspect the detection
                logger.info("Checking directory {} because power = {}".format(redshift[specie][count],power_rating[specie][count]))
                call("gv {}/vel_plot1.eps & gv {}/vel_plot3.0.eps".format(redshift[specie][count],redshift[specie][count]),shell=True)

                # Ask user if good detection
                print "y = yes, n = no, m = more, a = all, q = quit, r = reset dir, h = help menu"
                answer = raw_input("Is this a good detection? (y/n/m/a/q/r)")
                answer.lower()

                finalanswer = questionlogic(answer)

                if ( finalanswer == False ):
                    call("mv {} FalsePositive/{}".format(redshift[specie][count],redshift[specie][count]),shell=True)
                    redshift[specie] = np.delete(redshift[specie],count)
                    ewlim[specie] = np.delete(ewlim[specie],count)
                    power_rating[specie] = np.delete(power_rating[specie],count)
                    count -= 1

            count += 1

    logger.info("Final redshift list: {}".format(redshift))
    logger.info("Final ewlimit list: {}".format(ewlim))
    logger.info("Final power rating list: {}".format(power_rating))

    # Open file to write
    outfile = open(writefile,'w')

    logger.info("Writing detections_real.dat")
    for specie in redshift:
        outfile.write("Detected {} {} lines:\n".format(len(redshift[specie]),specie))
        outfile.write("--------------------------------------------------------------------------\n")
        outfile.write("Redshift         EW_limit         Power_Rating\n")
        for num,item in enumerate(redshift[specie]):
            outfile.write("%10s   %6s   %4s\n" % (item,ewlim[specie][num],power_rating[specie][num]))
        outfile.write("--------------------------------------------------------------------------\n")
    outfile.close()
    
# Stuff broke
except Exception, e:
    logger.error('Critical Error. Program Exited.',exc_info=True)


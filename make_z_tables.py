# PYTHON CODE TO READ IN ANALYSIS DATA
# AND OUTPUT THE COMBINED SAMPLE DATA
# INTO TABLES:
#  1.) ews.table
#  2.) vel.table
#  3.) aod.table
# CALL AS: python make_z_tables.py list_of_JNAMES.txt
# CALL FROM: spectra parent directory (list_of_JNAMES.txt directory)

import string
import numpy as np
import sys
import os
from subprocess import call
np.set_printoptions(threshold=np.nan)

# Set log options here
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,filemode='w')

# Set log file output
handler = logging.FileHandler('Tables.log','w')
handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(handler)


try:

    # Read in list and store JNAME directories
    listfile = str(sys.argv[1])
    readfile = "detections_combined.dat"
    infile_ews = "sysanal.ews"
    infile_vel = "sysanal.vmom"
    infile_aod = "sysanal.aod"
    infile_zabs = "zabs.dat"
    outfile_MgII = "total_MgII.table"
    outfile_CIV = "total_CIV.table"
    header = "%16s%16s%16s%20s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n" % ("JNAME","ZABS","ZABSWAV","EWLIM","EW_R","SIG","VBAR","SIG","WIDTH","SIG","ASYM","SIG","TAU","SIG-","SIG+","LOGN","SIG-","SIG+")

    MgII_table = []
    CIV_table = []
    
    jnames = np.loadtxt(listfile,unpack=True,dtype="string")

    counter = 0
    for directory in jnames:
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
        os.chdir(directory)
        logger.info("=========== In directory {}. ===========".format(directory))
        if ( os.path.isfile(readfile) == False):
            logger.info("======= {} doesn't exist, skipping directory {} =======".format(readfile,directory))
            os.chdir("../../")
            counter += 1
            continue
        # Read in detections_combined.dat
        with open(readfile,'rb') as f:
            alltext = f.readlines()
            for row in alltext:
                words = row.split()
                check = len(words)
                if ( check == 1 ):
                    logger.info("Next because:{}".format(row))
                elif ( check == 4 ):
                    if ( words[2] == "MgII" ):
                        delimiter = "MgII"
                elif ( words[0] == "Redshift" ):
                    logger.info("Next because: {}".format(row))
                else:
                    logger.info("Found good row.")
                    redshift[delimiter].append(words[0])
                    ewlim[delimiter].append(words[1])

        logger.info("============= Detections read in. Mining Analysis Data. ============\n")

        # Go into directories and read in sysanal files
        for specie in redshift:
            count = 0
            for systems in redshift[specie]:
                if ( os.path.isdir(systems) == False):
                    logger.info("{} does not exist. Copying from UNCOMBINED.".format(systems))
                    call(["cp","-r","UNCOMBINED/{}".format(systems),"."])
                os.chdir(systems)
                # Check if Sysanal files exist
                if ( os.path.isfile(infile_ews) == True):
                    # Read in EWS file
                    #logger.info("Reading in EW file.")
                    num,ew_o,ew_o_sig,ew_r,ew_r_sig,dum,dum,dum,tranname = np.genfromtxt(infile_ews,unpack=True,skip_header=1,skip_footer=1,dtype='string')
                else:
                    print "No ews file, skipping {}".format(systems)
                    os.chdir("../")
                    continue
                if ( os.path.isfile(infile_vel) == True):
                    # Read in VEL file
                    #logger.info("Reading in VEL file.")
                    num,vbar,vbar_sig,width,width_sig,asym,asym_sig,tranname = np.genfromtxt(infile_vel,unpack=True,skip_header=1,skip_footer=1,dtype='string')
                else:
                    print "No vel file, skipping {}".format(systems)
                    os.chdir("../")
                    continue
                if ( os.path.isfile(infile_aod) == True):
                    # Read in AOD file
                    #logger.info("Reading in AOD file.")
                    num,tau,tau_sig_down,tau_sig_up,logn,logn_sig_down,logn_sig_up,tranname = np.genfromtxt(infile_aod,unpack=True,skip_header=1,skip_footer=1,dtype='string')
                else:
                    print "No aod file, skipping {}".format(systems)
                    os.chdir("../")
                    continue
                if ( os.path.isfile(infile_zabs) == True):
                    # Read in zabs.dat
                    zabs,zabswav,dum = np.genfromtxt(infile_zabs,unpack=True,dtype='string')
                else:
                    print "No zabs.dat, skipping {}".format(systems)
                    os.chdir("../")
                    continue

                #logger.info("Read in all sysanal files.")

                # Combine relevant data into data table arrays
                # JNAME ZABS EW_R EW_R_SIG VBAR VBAR_SIG WIDTH WIDTH_SIG ASYM ASYM_SIG TAU TAU_SIG_- TAU_SIG_+ LOGN LOGN- LOGN+
                (instrument,jname) = jnames[counter].split('/')
                if ( tranname == "MgII2796" ): #18 columns
                    MgII_table.append("%16s%16s%16s%20s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s" % (jname,zabs,zabswav,ewlim[specie][count],ew_r,ew_r_sig,vbar,vbar_sig,width,width_sig,asym,asym_sig,tau,tau_sig_down,tau_sig_up,logn,logn_sig_down,logn_sig_up))
                if ( tranname == "CIV1548" ):
                    CIV_table.append("%16s%16s%16s%20s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s" % (jname,zabs,zabswav,ewlim[specie][count],ew_r,ew_r_sig,vbar,vbar_sig,width,width_sig,asym,asym_sig,tau,tau_sig_down,tau_sig_up,logn,logn_sig_down,logn_sig_up))

                os.chdir("../")
                count += 1

            # End for systems loop
        
        counter += 1

        os.chdir("../../")

        # End for specie in redshift loop

    # End for directory in jnames loop

    MgII_out = open(outfile_MgII,'w')
    MgII_out.write(header)
    for row in MgII_table:
        MgII_out.write("{}\n".format(row))

    CIV_out = open(outfile_CIV,'w')
    CIV_out.write(header)
    for row in CIV_table:
        CIV_out.write("{}\n".format(row))

except Exception, e:
    logger.error('Critical Error. Program Exited.',exc_info=True)

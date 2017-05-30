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
handler = logging.FileHandler('gwz.log','w')
handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(handler)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

try:

    # Read in list and store JNAME directories
    readfile = {}
    listfile = str(sys.argv[1])
    readfile["MgII"] = "MgII_ewlim.data"
    readfile["CIV"] = "CIV_ewlim.data"
    outfile_MgII = "gwz_MgII.table"
    outfile_CIV = "gwz_CIV.table"

    jnames = np.loadtxt(listfile,unpack=True,dtype="string")

    redshift = {}
    EW_limit = {}

    counter = 0
    for directory in jnames:
        # Split instrument/JNAME
        instrument,jname = directory.split("/")

        os.chdir(directory)
        logger.info("=========== In directory {}. ===========".format(directory))

        # Array definitions
        redshift[directory] = {}
        EW_limit[directory] = {}

        for specie in readfile:
            # Check if MgII_ewlim.data and CIV_ewlim.data exist
            if ( os.path.isfile(readfile[specie]) == False):
                logger.info("======= {} doesn't exist. =======".format(readfile[specie]))
                continue
            if ( os.path.isfile(readfile[specie]) == False):
                logger.info("======= {} doesn't exist. =======".format(readfile[specie]))
                continue

            redshift[directory][specie] = [] 
            EW_limit[directory][specie] = []

            # Read in specie_ewlim.data files and store in big array
            with open(readfile[specie],'rb') as f:
                alltext = f.readlines()
                for row in alltext:
                    z,ewlim = row.split()
                    redshift[directory][specie].append(float(z))
                    EW_limit[directory][specie].append(float(ewlim) / (1.+float(z)))

        # End for specie in readfile

        logger.info("============= EWLIM files read in and stored for {}. ============\n".format(directory))

        os.chdir("../../")

    # End for directory in jnames

    logger.info("============= Done reading in all files. Calculating g(w,z) now. =============\n")

    w_grid = np.arange(0.005,0.2,0.005)
    z_grid = np.arange(0.13,4.89,0.01)

    gwz_MgII = {}
    gwz_CIV = {}
    for w in w_grid:
        w = round(w,4)
        gwz_MgII[w] = {}
        gwz_CIV[w] = {}
        for z in z_grid:
            z = round(z,4)
            gwz_MgII[w][z] = 0
            gwz_CIV[w][z] = 0
            for directory in redshift:
                for specie in redshift[directory]:
                    if ( (z > min(redshift[directory][specie])) & (z < max(redshift[directory][specie]))):
                        
                        redshift[directory][specie] = np.array(redshift[directory][specie])
                        z_val = find_nearest(redshift[directory][specie],z)
                        z_index = np.where(redshift[directory][specie] == z_val)[0]

                        if ( (EW_limit[directory][specie][z_index] < w) and (specie == "MgII")):
                            gwz_MgII[w][z] += 1
                        if ( (EW_limit[directory][specie][z_index] < w) and (specie == "CIV")):
                            gwz_CIV[w][z] += 1
                        

    logger.info("=========== Calculated g(w,z). ===============\n")

    sys.exit()

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




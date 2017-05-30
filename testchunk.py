import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
from fzl_functions import *
from scipy import stats
import timeit
from subprocess import call
np.set_printoptions(threshold=np.nan)

# Read in atoms.dat to get atomic data
# EDIT: Path length for the location of your atoms.dat file
transition_library = np.genfromtxt('/home/achernar-data/mathes/UVES_SQUAD/atoms_UVES.dat',unpack=True,dtype=None)

# Transition Wavelengths and their power/noise arrays
species = {}
power = {}
noise = {}
working = {}
sigma = {}

# EDIT: Add transitions. These are the lines you're looking for. 
search_transitions = ['CIV', 'MgII']

# Loop through atoms.dat and match with search_transitions
# Initialize arrays to search for transitions
for search_trans in search_transitions:
    species[search_trans] = {}
    for transition in transition_library:
        [atom,tranwav] = re.findall('[a-zA-Z]+|\\d+',transition[0])
        if ( atom == search_trans ):
            species[atom][tranwav] = transition[1]
            power[tranwav] = []
            noise[tranwav] = []
            working[tranwav] = []
            sigma[tranwav] = []
                

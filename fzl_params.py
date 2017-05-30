# Path to atomic data file
transition_library = np.genfromtxt('/home/achernar-data/mathes/VULTURE/UVES_SQUAD/atoms_UVES.dat',unpack=True,dtype=None)

# Ionic transitions to find
search_ions = ['CIV', 'MgII']

# Spectral resolution
spectral_res = 45000.

# Redshift resolution
delta_z = 0.00001

# Can tune how strict you are about the doublet ratio: if weak_line/strong_line < doublet_ratio, throw out.
# Higher values of dub_ratio = throw out more lines. DO NOT EXCEED 0.5 (unphysical)
# Can comment lines out if you don't want this constraint applied.
dub_ratio = 0.3

# Can edit velocity window for output plots
velwindow = 500. #km/s


import string
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from numpy.fft import rfft, irfft
from scipy.stats import norm
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
        v = round(v,4)
        expect = round(expect,4)
        print v,expect
        if (v == expect) or (expect == 0.):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result

#================================================================
#======================== Initialization ========================
#================================================================

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

# Transition Wavelengths
LYA = 1215.6
MgII2796 = 2796.353786 
MgII2803 = 2803.530982
CIV1548 = 1548.195
CIV1551 = 1550.770
LIGHTSPEED = 299792.458 # km/s
# +/- 0.5\AA = 8 pix

# Read in spectra
wav,dum,flux,sig,err,cont = np.loadtxt("J220852-194359.data",unpack=True)
#eqwav,eqwidth,eqsig = np.loadtxt("J005758-264314.ewspec",unpack=True)
#wav,dum,flux,sig,err,cont = np.loadtxt("J014333-391700.data",unpack=True)

# Read in zem.dat
zem = np.loadtxt("zem.dat")

# Create flux arrays for data and error
norm_flux = flux/cont
norm_err = err/cont
norm_sig = sig/cont

#======================================================================
#======================== Set up Search Window ========================
#======================================================================

# Search redward of the Lya emission
blue_limit = LYA * (1. + zem)
red_limit = MgII2803 * (1. + zem)
search_flux = norm_flux[np.where((wav > blue_limit) & (wav < red_limit))]
search_wav = wav[np.where((wav > blue_limit) & (wav < red_limit))]
search_err = norm_err[np.where((wav > blue_limit) & (wav < red_limit))]
search_sig = norm_sig[np.where((wav > blue_limit) & (wav < red_limit))]

# Clean cosmic rays
search_flux[np.where(search_flux > 1.1)] = 1.

# Avoid Telluric features
# 7594 - 7700 feature
search_flux[np.where((search_wav > 7594.) & (search_wav < 7700.))] = 1.
# 6868 - 6929 feature
search_flux[np.where((search_wav > 6868.) & (search_wav < 6929.))] = 1.

# Remove negative errors
search_err[np.where(search_err < 0.)] = 0.

# Create arrays to hold the cross-correlation result (to be plotted/analyzed)
result_MgII2796 = []
np.array(result_MgII2796)
result_MgII2803 = []
np.array(result_MgII2803)
result_CIV1548 = []
np.array(result_CIV1548)
result_CIV1551 = []
np.array(result_CIV1551)

# Create arrays to hold the noise cross-correlation
noise_MgII2796 = []
np.array(noise_MgII2796)
noise_MgII2803 = []
np.array(noise_MgII2803)
noise_CIV1548 = []
np.array(noise_CIV1548)
noise_CIV1551 = []
np.array(noise_CIV1551)

# Search redshift array from z = 0 to z = zem
redshift = arange(0,zem,0.0001)

# Normalize search_flux
#search_flux = (search_flux - np.mean(search_flux))/(np.std(search_flux) * len(search_flux))
#search_err = (search_err - np.mean(search_err))/(np.std(search_err) * len(search_err))
search_flux = (search_flux - 1.)/(np.std(search_flux) * len(search_flux))
search_err = (search_err - np.mean(search_err))/(np.std(search_err) * len(search_err))

#==================================================================
#======================== Do the Filtering ========================
#==================================================================

# Loop redshift, filter is 40 pixels in width (+/- 20)
for z in redshift:

    search_filt = np.ones(search_flux.size)

    MgII2796central = find_nearest(search_wav,MgII2796*(1.+z))
    MgII2796central = np.where(search_wav == MgII2796central)
    MgII2796central = MgII2796central[0]
    search_filt[MgII2796central-20:MgII2796central+20] = 0.
    search_filt = (search_filt - np.mean(search_filt))/(np.std(search_filt))
    result_MgII2796.append(np.correlate(search_flux,search_filt))
    noise_MgII2796.append(np.correlate(search_err,search_filt))

    # Do it all again 
    search_filt = np.ones(search_flux.size)
    MgII2803central = find_nearest(search_wav,MgII2803*(1.+z))
    MgII2803central = np.where(search_wav == MgII2803central)
    MgII2803central = MgII2803central[0]
    search_filt[MgII2803central-20:MgII2803central+20] = 0.
    search_filt = (search_filt - np.mean(search_filt))/(np.std(search_filt))
    result_MgII2803.append(np.correlate(search_flux,search_filt))
    noise_MgII2803.append(np.correlate(search_err,search_filt))

    #CIV1548central = find_nearest(search_wav,CIV1548*(1.+z))
    #CIV1548central = np.where(search_wav == CIV1548central)
    #CIV1548central = CIV1548central[0]

    #CIV1551central = find_nearest(search_wav,CIV1551*(1.+z))
    #CIV1551central = np.where(search_wav == CIV1551central)
    #CIV1551central = CIV1551central[0]

    #search_filt[MgII2803central-10:MgII2803central+10] = 0.
    #search_filt[CIV1548central-10:CIV1548central+10] = 0.
    #search_filt[CIV1551central-10:CIV1551central+10] = 0.

result_MgII2796 = np.divide(result_MgII2796,noise_MgII2796)
result_MgII2803 = np.divide(result_MgII2803,noise_MgII2803)
criterion_MgII = np.where((result_MgII2796 > 4) & (result_MgII2803 > 2))
final_MgII = result_MgII2796[criterion_MgII[0]]
final_redshift = redshift[criterion_MgII[0]]
print final_redshift

# Group the adjacent redshifts to find the detection centers
grouped_redshifts = group_consecutives(final_redshift)

# Find line centers
line_redshifts=[]
rounded_redshifts=[]
for lines in grouped_redshifts:
    # Skip over single element 'detections
    if len(lines) == 1:
        continue
    center = np.mean(lines)
    line_redshifts.append(round(center,4))

print line_redshifts

# Align with common zero point and do chi-sqr? Or something.
for z in line_redshifts:
    velocity = ((LIGHTSPEED * z * z) + (2. * LIGHTSPEED * z)) / ((z * z) + (2. * z) + 2.)
    

# Fit gaussians to the perceived detections and see if they're real
#param = norm.fit(cut)
#fit_window=[]
#fit_params=[]
#for element in line_redshifts:
#    MgII2796central = find_nearest(wav,MgII2796*(1. + element))
#    MgII2796central = np.where(wav == MgII2796central)
#    MgII2796central = MgII2796central[0]
#    fit_wav = wav[MgII2796central-25:MgII2796central+25]
#    fit_flux = norm_flux[MgII2796central-40:MgII2796central+40]
#    plot(fit_wav,fit_flux)
#    show()

# Plot the things
plot(redshift,result_MgII2796)
plot(redshift,result_MgII2803)
plot(final_redshift,final_MgII,'ro')
xlim(0.5,2.8)
ylim(-5,50)
xlabel('Redshift')
ylabel('Result')
grid(True)
#savefig("Result.png")
show()

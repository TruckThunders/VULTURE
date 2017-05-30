import numpy as np
# Set log options here
import logging
func_logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,filemode='a')

# Set log file output
func_handler = logging.FileHandler('False_positive.log','w')
func_handler.setLevel(logging.INFO)

# Create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
func_handler.setFormatter(formatter)

# Add the handlers to the logger
func_logger.addHandler(func_handler)

#===========================================================
#======================== Functions ========================
#===========================================================

# Check if all elements in a list are equal
def checkEqual(iterator):
    try:
        iterator = iter(iterator)
        first = next(iterator)
        return all(first == rest for rest in iterator)
    except StopIteration:
        return True

def exitlog(qsoname,stop):
    if ( stop == 0 ):
        open("../start_log.txt","a").write("%s \n" % qsoname)
    else:
        open("../end_log.txt","a").write("%s \n" % qsoname)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def group_consecutives(vals, step=0.0001):
    """Return list of consecutive lists of numbers from vals (number list)."""
    run = []
    result = [run]
    expect = 0.
    for v in vals:
        v = round(v,5)
        expect = round(expect,5)
        if (v == expect) or (expect == 0.):
            run.append(v)
        else:
            run = [v]
            result.append(run)
        expect = v + step
    return result

def asymmetry_check(line_redshifts, ion_spectrum, redshift):
    # Roll a ball around the line_redshifts to throw out false positives
    count = 0
    for detection in line_redshifts:
        leftup = 0
        leftdown = 0
        rightup = 0
        rightdown = 0
        # Define line center and establish max value
        center = np.where(redshift == find_nearest(redshift,detection))
        center = int(center[0])
        # Check for bigger max
        for i in range(-5,5):
            max = np.copy(ion_spectrum[center])
            maxcheck = np.copy(ion_spectrum[center+i])
            if (maxcheck > max):
                max = np.copy(maxcheck)    
        if ( max > 0.9 ):
            count += 1
            continue
        # range(-10,10) is how far to search left/right for asymmetry
        for i in range(-10,10):
            # Roll ball around max value and see if it checks out
            check = ion_spectrum[center + i]
            if i < 0:
                if check < max:
                    leftdown += 1
                if check > max: 
                    leftup += 1
            if i > 0:
                if check < max:
                    rightdown += 1
                if check > max:
                    rightup += 1
        # Delete lines based upon asymmetry around the line center
        if ((abs(leftup - rightup) > 4) & (abs(leftdown - rightdown) > 4)):
            line_redshifts = np.delete(line_redshifts,count)
            func_logger.info("Deleted line (asymmetry): {}".format(detection))
            count -= 1
        count += 1
    return line_redshifts

# ion_spectrum1 corresponds to the ion with the larger ionization potential
# *** if weak_line/strong_line < doublet_ratio, throw out. 
def line_ratio_check(specie, line_redshifts, ion_spectrum1, ion_spectrum2,redshift,doublet_ratio):
    # Throw out false positives based upon line ratio
    count = 0
    for detection in line_redshifts:
        center = np.where(redshift == find_nearest(redshift,detection))
        center = int(center[0])
        # Delete lines based upon line ratio
        for offset in range(3):
            if ((ion_spectrum1[center+offset] > 0.95) & (ion_spectrum2[center+offset] > 0.95) |
                (ion_spectrum1[center-offset] > 0.95) & (ion_spectrum2[center-offset] > 0.95)):
                func_logger.info("Both members of {} doublet saturated: {}".format(specie,detection))
            elif (((ion_spectrum1[center+offset]+0.002) < ion_spectrum2[center+offset]) | 
                  ((ion_spectrum1[center-offset]+0.002) < ion_spectrum2[center-offset])):
                line_redshifts = np.delete(line_redshifts,count)
                #print "Deleted line (Main Transition >> Secondary): {}".format(detection)
                func_logger.info("Deleted {} line (Main Transition < Secondary): {}".format(specie,detection))
                count -= 1
                break
            elif (((ion_spectrum2[center+offset]/ion_spectrum1[center+offset]) < doublet_ratio) |
                  ((ion_spectrum2[center-offset]/ion_spectrum1[center-offset]) < doublet_ratio)):
                line_redshifts = np.delete(line_redshifts,count)
                #print "Deleted line (Bad doublet ratio): {}".format(detection)
                func_logger.info("Deleted {} line (Bad doublet ratio): {}".format(specie,detection))
                count -= 1
                break
        count += 1
    return line_redshifts

# Insert fake lines to test detection sensitivity
def insert_lines(spectrum,wavelengths,resolution,z,eqwidth,lambda_r):

    testwidth = 0.
    # Calculate gaussian parameters
    b = lambda_r * (1+z)
    c = (b / resolution) * (1. / (2.*np.sqrt(2.*np.log(2))))
    a = eqwidth / (c * np.sqrt(2. * np.pi))
    npix = (resolution * 4.) / b
    half_width = int(np.floor(npix/2.))
    
    middle = find_nearest(wavelengths,b)
    middle = np.where(wavelengths == middle)
    middle = middle[0]
    
    for i in range(-half_width,half_width):
        spectrum[middle+i] = 1.
        spectrum[middle+i] -= a * np.exp(-((wavelengths[middle+i] - b)**2.)/(2.*c*c))
        testwidth += (1.-spectrum[middle+i])*(wavelengths[middle]-wavelengths[middle-1])

    return (a,c,spectrum)

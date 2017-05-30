# CALL BY USING:                                                    #
#  perl parfile.pl [WAV_MIN]  [WAV_MAX]                             #
#                   $ARGV[0]   $ARGV[1]                             #

if ( $ARGV[0] eq '' ) { 

    die "FRAUD DETECTED! No Minimum Wavelength"; 
    
}

elsif ( $ARGV[1] eq '') {

    die "FRAUD DETECTED! No Maximum Wavelength"; 

}

$wmin = $ARGV[0];
$wmax = $ARGV[1];

# Write specsynth parameter file
open(PAR,">specsynthUVES.par");
	
print PAR "C***************** input parameters for specsynth *****************************\n";
print PAR "$wmin   wave_min - minimum wavelength\n";
print PAR "$wmax   wave_max - maximum wavelength\n";
print PAR "0.005 dwave - linear dispersion of sampling in angstroms\n";
print PAR "40000.  R_fac - the spectrograph resolution R=lambda/(Delta lambda)\n";
print PAR "1.14    slit - the spectrograph slit width in arcseconds\n";
print PAR "1       conflag - toggle convolution 1=ON 0=OFF\n";
print PAR "3.0     conwindo -  \# of inst. profile sigma to include in convolution\n";
print PAR "3       resfac - (odd integer) sample rate factor for convolution arrays\n";
print PAR "0.     snr - signal to noise ratio\n";

close(PAR);

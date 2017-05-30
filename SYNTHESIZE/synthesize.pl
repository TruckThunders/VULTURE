# CALL BY USING:                                                    #
#  perl synthesize.pl 

# Get first line of line.inp for system redshift $zsys
open my $file, '<', "lines.inp"; 
my $zsys = <$file>; 
close $file;

print "zsys = $zsys\n";

# Read the rest of lines.inp to get line parameters
open(INPUT,"<lines.inp");
@LINES= <INPUT>;
close(INPUT);
chomp(@LINES);
shift(@LINES);

# Get atomic data from atoms.dat
open(ATOMIC,"<atoms.dat");
@ATOMS= <ATOMIC>;
close(ATOMIC);
chomp(@ATOMS);

foreach (@LINES) {

    ($tranname, $redshift, $column, $bparam) = split();

    foreach (@ATOMS) {

	#print "Searching atoms.dat\n";

	($maintran,$centralwave,$oscillator,$dum,$dum,$dum) = split;

	if ( $tranname eq $maintran ) { 
	    
	    print "Found $tranname at $centralwave\n";
	    last;
	    
	}
	
    }    
    
    $wmin = ($centralwave - 4.) * (1. + $zsys);
    $wmax = ($centralwave + 4.) * (1. + $zsys);
    
    print "Writing specsynthUVES\.par using $wmin and $wmax\n";

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
    print PAR "20.     snr - signal to noise ratio\n";
    
    close(PAR);

# Run Specsynth for a given line

    system("/home/achernar/mathes/PROGRAMS/specsynth/specsynth2 lines.inp specsynthUVES.par lines");
    system("cp lines\.out $maintran\.data");
    
}

system("cat *\.data > smashed\.synth");

#!/usr/local/bin/perl

#*******************************************************************************
# Check for Fits and Dump directories and make them if they don't exist
#*******************************************************************************

unless (-d "Fits") { 
    system("mkdir Fits");
    system("cp *.fits Fits/.");
}
unless (-d "Dump") { system("mkdir Dump");}
unless (-e "quickdump.pro") { system("cp /achernar-data/mathes/UVES_SQUAD/quickdump.pro .");}

#*******************************************************************************
# FILENAME_up.fits does not exist - No recent UVES_popler run
#*******************************************************************************

unless (grep -f, glob '*up.fits') {

    print "_up.fits file does not exist. Dumping .fits file.\n";

# Create fits.list with new fit file
    system("ls *.fits > fits.list");

# Use IDL script to Dump .spec, .sig, .err, and .cont files
    
    system("idl -e quickdump");
    
# Read in quasar from fits.list
    open(A,"<fits.list");
    @FILES = <A>;
    close(A);
    chomp(@FILES);
    $fitsfile = @FILES[0];
    
# Split Jname and .fits extension
    ($qsoname,$extension) = split(/\./,$fitsfile);
    
    chdir("Dump");
    
    open(A,"<$qsoname\.spec");
    @SPEC = <A>;
    close(A);
    chomp(@SPEC);
    
    open(B,"<$qsoname\.cont");
    @CONT = <B>;
    close(B);
    chomp(@CONT);
    
    open(C,"<$qsoname\.sig");
    @SIG = <C>;
    close(C);
    chomp(@SIG);
    
    open(D,"<$qsoname\.err");
    @ERR = <D>;
    close(D);
    chomp(@ERR);
    
# Continuum read
    $i = 0;
    foreach $line (@CONT){
	($dum,$cont[$i]) = split(" ",$line);
	$i++;
    }
    
# Spectrum read and un-normalize flux
    $i = 0;
    foreach $line (@SPEC){
	($wav[$i], $flux[$i]) = split(" ",$line);
	$flux[$i] *= $cont[$i];
	$wav[$i] = 10.**$wav[$i];
	$i++;
    }
    
# Sigma read and normalize sigma
    $i = 0;
    foreach $line (@SIG){
	($dum,$sig[$i]) = split(" ",$line);
	$sig[$i] = -1. if ($sigma[$i] < -1.);
	$sig[$i] *= $cont[$i] unless ($sig[$i] == -1.);
	$i++;
    }
    
# Error read and normalize error
    $i = 0;
    foreach $line (@ERR){
	($dum,$err[$i]) = split(" ",$line);
	$err[$i] = -1. if ($err[$i] < -1.);
	$err[$i] *= $cont[$i] unless ($err[$i] == -1.);
	$i++;
    }
    
    chdir("../");
    
    $i=0;
    open(Z,">$qsoname\.data");
    foreach $line (@SPEC){
	print Z "$wav[$i] 0.00 $flux[$i] $sig[$i] $err[$i] $cont[$i]\n";
	$i++;
    }
    close(Z);

    print "$qsoname successfully combined.\n";
    
}

#*******************************************************************************
# FILENAME_up.fits exists - Recent UVES_popler run
#*******************************************************************************

else {

    print "_up.fits file exists. Dumping *up.fits file.\n";

# Create fits.list with new fit file
    system("ls *up.fits > fits.list");
    
# Use IDL script to Dump .spec, .sig, .err, and .cont files
    
    system("idl -e quickdump");
    
# Read in quasar from fits.list
    open(A,"<fits.list");
    @FILES = <A>;
    close(A);
    chomp(@FILES);
    $fitsfile = @FILES[0];
    
# Split Jname and .fits extension
    ($upname,$extension) = split(/\./,$fitsfile);
    ($qsoname,$up) = split(/\_/,$upname);
    
    print "$qsoname $upname \n";
    
    chdir("Dump");
    
    open(A,"<$upname\.spec");
    @SPEC = <A>;
    close(A);
    chomp(@SPEC);
    
    open(B,"<$upname\.cont");
    @CONT = <B>;
    close(B);
    chomp(@CONT);
    
    open(C,"<$upname\.sig");
    @SIG = <C>;
    close(C);
    chomp(@SIG);
    
    open(D,"<$upname\.err");
    @ERR = <D>;
    close(D);
    chomp(@ERR);
    
# Continuum read
    $i = 0;
    foreach $line (@CONT){
	($dum,$cont[$i]) = split(" ",$line);
	$i++;
    }
    
# Spectrum read and un-normalize flux
    $i = 0;
    foreach $line (@SPEC){
	($wav[$i], $flux[$i]) = split(" ",$line);
	$flux[$i] *= $cont[$i];
	$wav[$i] = 10.**$wav[$i];
	$i++;
    }
    
# Sigma read and normalize sigma
    $i = 0;
    foreach $line (@SIG){
	($dum,$sig[$i]) = split(" ",$line);
	$sig[$i] = -1. if ($sigma[$i] < -1.);
	$sig[$i] *= $cont[$i] unless ($sig[$i] == -1.);
	$i++;
    }
    
# Error read and normalize error
    $i = 0;
    foreach $line (@ERR){
	($dum,$err[$i]) = split(" ",$line);
	$err[$i] = -1. if ($err[$i] < -1.);
	$err[$i] *= $cont[$i] unless ($err[$i] == -1.);
	$i++;
    }
    
    chdir("../");
    
    $i=0;
    open(Z,">$qsoname\.data");
    foreach $line (@SPEC){
	print Z "$wav[$i] 0.00 $flux[$i] $sig[$i] $err[$i] $cont[$i]\n";
	$i++;
    }
    close(Z);
    
    print "$qsoname successfully combined.\n";
    
# Nuke _up files
    
#print "$upname\.fits $qsoname\.fits\n";
    
    system("mv $upname\.fits $qsoname\.fits");
    
}

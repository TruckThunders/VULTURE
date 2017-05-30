#!/usr/local/bin/perl                                                          #
#                                                                              #
#      USAGE: batchcaller.pl list_of_dirs.txt [clean/check/count]              #
#                                                                              #
#             No [] arg = find_z_lines (detect MgII/CIV lines in spectra)      #
#             [clean] = remove detections from directory to start over         #
#             [check] = verify detections using check_z_lines                  #
#             [count] = count number of detections for MgII/CIV                #
#                                                                              #
################################################################################

#$pythonpath = "/home/local/Canopy_64bit/User/bin/python";   # Dept. CentOS
$pythonpath = "/opt/local/bin/python";   # Mac OS
#$findpath = "/home/achernar/mathes/PROGRAMS/x_corr/find_z_lines.py";
$findpath = "/Users/Nigel/SOFTWARE/x_corr/find_z_lines.py";
#$checkpath = "/home/achernar/mathes/PROGRAMS/x_corr/check_z_lines.py";
$checkpath = "/Users/Nigel/SOFTWARE/x_corr/check_z_lines.py";

# (1) quit unless we have the correct number of command-line args
$num_args = $#ARGV + 1;
if ($num_args != 1 && $num_args != 2) {
    print "\nUsage: batchcaller.pl list_of_dirs.txt [clean/check/count]\n";
    exit;
}

if ($num_args == 2) {
    $dirlist = $ARGV[0];
    $commandflag = $ARGV[1];

    # Clean out directories
    if ( $commandflag eq "clean" ){
	print "Cleaning directories\n";

	my $filename = "$dirlist";
	open ( my $fh, '<:encoding(UTF-8)', $filename)
	    or die "Could not open file $filename $!";

	while (my $row = <$fh>) {
	    chomp $row;
	    chdir("$row");
	    print "\n*************************************************\n";
	    print "         In directory $row.                \n";
	    print "*************************************************\n";
	    system("$pythonpath $findpath clean");
	    chdir("../");
	}
    }
    # Check if the detections automatically found by find_z_lines are real
    if ( $commandflag eq "check" ) {
	print "Checking detections\n";

	my $filename = "$dirlist";
	open ( my $fh, '<:encoding(UTF-8)', $filename)
	    or die "Could not open file '$filename' $!";

	while (my $row = <$fh>) {
	    chomp $row;
	    chdir("$row");
	    print "\n*************************************************\n";
	    print "         In directory $row.                \n";
	    print "*************************************************\n";
	    system("$pythonpath $checkpath");
	    chdir("../");
	}
    }
    if ( $commandflag eq "count" ) {
	print "Counting detections\n";

	$totalcountCIV = 0;
	$totalcountMgII = 0;

	my $filename = "$dirlist";
	open ( my $fh, '<:encoding(UTF-8)', $filename)
	    or die "Could not open file '$filename' $!";

	while (my $row = <$fh>) {
	    chomp $row;
	    chdir("$row");
	    open(A,"<detections.dat");
	    while (<A>){
		($dum,$count,$delimiter,$dum) = split(" ",$_);
		if ( $delimiter eq "CIV" ){
		    $totalcountCIV += $count;
		}
		if ( $delimiter eq "MgII" ){
		    $totalcountMgII += $count;
		}
	    }
	    close(A);
	    chdir("../");
	}
	print "Detected $totalcountCIV CIV lines.\n";
	print "Detected $totalcountMgII MgII lines.\n\n";
	$scaledcountCIV = 0.75 * $totalcountCIV;
	$scaledcountMgII = 0.75 * $totalcountMgII;
	print "Expected $scaledcountCIV CIV lines.\n";
	print "Expected $scaledcountMgII MgII lines.\n\n";
    }
    else { die "\nUsage: batchcaller.pl list_of_dirs.txt [clean,check,count]\n"; }
}

# Run find_z_lines (automatic detection of absorption lines)
if ($num_args == 1) {
    $dirlist = $ARGV[0];
    print "Detecting lines in directories listed in $dirlist\n";

    my $filename = "$dirlist";
    open ( my $fh, '<:encoding(UTF-8)', $filename)
	or die "Could not open file $filename $!";

    while (my $row = <$fh>) {
	chomp $row;
	chdir("$row");
	print "\n*************************************************\n";
	print "         In directory $row.                \n";
	print "*************************************************\n";
	system("$pythonpath $findpath $row\.data");
	chdir("../");
    }
}



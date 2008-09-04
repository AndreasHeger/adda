# calculate selectivity based on the
# output of OutputStatisticsClustering
# note: several repeat units of a PFAM domain might be 
#	within a single domain (for example 00818, ice nucleation repeat)
# The ratio of anno/total would be zero.

use strict;

my %anno;
my %total;
my $family;

my $last_did = 0;

my $t = 0;
my $n = 0;

print <<END_OF_TEXT;
# Selecitivity of a clustering
# FAMILY:	family that is compared
# ANNO:	        number of domains in cluster annotated by assigned reference family
# TOTAL:	number of total annotated domains in cluster
# SENSI:	sensitivity
# REFF:		assigned reference family
# FAMILY\tANNO\tTOTAL\tSENSI\tREFF
END_OF_TEXT


while (<STDIN>) {
    if (/# master/) {
	print $_; 
	next;
    }

    if (/# reference/) {
	print $_; 
	next;
    }

    next if (/^\#/);
    chop();

    my ($did, $nunits, $aunits, $nseqs, $aseqs, $length, 
	$runits, $tunits, $rseqs, $tseqs, $alength, $aovl, $anno1, $anno2) = split(/\t/);
    
    if ($did) {
	
	my $family = $anno1."\t".$anno2;
	next unless ($runits);
	my $r = $runits / $aunits;
	if ($r > 1) {$r = 1};
	printf ("$did\t$runits\t$aunits\t%5.2f\t$family\n", $r);	
	$t += $r;
	$n++;
    }
}


print "# average selectivity: " . $t/$n . "\n" ;



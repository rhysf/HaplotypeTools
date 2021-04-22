#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF from step 1 of HaplotypeTools> -a <alignment file from step 2 of HaplotypeTools>\n
Optional: -s\tSample number [0]\n
Outputs:  -o\tOutput [opt_v-phased-opt_s]\n";
our($opt_a, $opt_o, $opt_s, $opt_v);
getopt('aosv');
die $usage unless ($opt_a && $opt_v);
if(!defined $opt_s) { $opt_s = 0; }
if(!defined $opt_o) { $opt_o = "$opt_v-phased-$opt_s"; }

# Phase
my $phased_VCF_hash = vcfphase::phase_reads($opt_a, $opt_v, $opt_s);

# Print
open my $ofh, '>', $opt_o or die "Cannot open $opt_o : $!\n";
foreach my $positions(sort { $a <=> $b } keys %{$phased_VCF_hash}) {
	my $line = $$phased_VCF_hash{$positions};
	print $ofh "$line\n";
}
close $ofh;

#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -a <alignment file from VCF_and_BAM_to_phased_VCF.pl> -v <VCF subset from VCF_and_BAM_to_phased_VCF.pl> -c <contig>
Outputs: -o\tOutput [opt_v-phased]\n";
our($opt_a, $opt_c, $opt_g, $opt_o, $opt_p, $opt_q, $opt_v);
getopt('acov');
die $usage unless ($opt_a && $opt_v && $opt_c);
if(!defined $opt_o) { $opt_o = "$opt_v-phased"; }

# Phase
my $phased_VCF_hash = vcfphase::phase_reads($opt_a, $opt_v, $opt_c);

# Print
open my $ofh, '>', $opt_o or die "Cannot open $opt_o : $!\n";
foreach my $positions(sort { $a <=> $b } keys %{$phased_VCF_hash}) {
	my $line = $$phased_VCF_hash{$positions};
	print $ofh "$line\n";
}
close $ofh;

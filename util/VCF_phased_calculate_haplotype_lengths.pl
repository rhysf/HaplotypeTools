#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;
use lib "/home/unix/rfarrer/perl5/lib/perl5/";
use Hash::Merge qw(merge);

### rfarrer@broadinstitue.org

# Opening commands 
my $usage = "Usage: perl $0 <Phased VCF file (with optional additional VCFs separated by comma)>\n";
die $usage if (@ARGV != 1);
my ($vcf_file) = @ARGV;
my @vcf_files = split /,/, $vcf_file;

# Go through each VCF get phase info
my (%counts, %types_found);
foreach my $file(@vcf_files) {
	my ($counts_temp, $types_found_temp) = vcflines::parse_VCF($file, 'none', 'all', 'n');	
	%counts = %{merge(\%counts, $counts_temp)};
	%types_found = %{merge(\%types_found, $types_found_temp)};
}
#print "Dumper:\n";
#print Dumper(%types_found);

# Merge phase group lengths
my %phase_group_lengths;
foreach my $file(@vcf_files) {
	warn "Calculating summary from $file...\n";
	foreach my $contig(sort keys %{$types_found{'Phase_group_lengths'}{$file}}) {
		foreach my $length(sort keys %{$types_found{'Phase_group_lengths'}{$file}{$contig}}) {
			my $tally = $types_found{'Phase_group_lengths'}{$file}{$contig}{$length};
			#warn "making hash of hap length $length, tally = $tally\n";
			$phase_group_lengths{$length} += $tally;
		}
	}
}

# Summary 1
my ($min, $max, $mean, $total, $num) = vcflines::summarise_phases(\%phase_group_lengths);
warn "Min length\t$min\n";
warn "Max length\t$max\n";
warn "Mean length\t$mean\n";
warn "Mean total\t$total\n";
warn "Number of haplotypes\t$num\n\n";

# Summary 2
print "HapLength\tTally\n";
foreach my $length(sort keys %phase_group_lengths) {
	my $tally = $phase_group_lengths{$length};	
	print "$length\t$tally\n";
}

#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;
use Hash::Merge qw(merge);
use Getopt::Std;
use Data::Dumper;

### r.farrer@exeter.ac.uk

# Opening commands
my $usage = "Usage: perl $0 -v <Phased VCF file (with optional additional VCFs separated by comma)>\n
Optional: -s <sample names (with optional additional sample names separated by comma)> []
	  -e Output results for every sample (y/n) [n]\n";
our($opt_e, $opt_v, $opt_s);
getopt('evs');
die $usage unless ($opt_v);
if(!defined $opt_e) { $opt_e = 'n'; }
my @vcf_files = split /,/, $opt_v;

# Go through each VCF get phase info
my $merger = Hash::Merge->new('LEFT_PRECEDENT');
my (%counts, %types_found);
foreach my $file(@vcf_files) {
	die "Cannot open $file : $!" unless (-e $file);

	# If samples are specified, only process those wherever they are found
	if((defined $opt_s) && ($opt_e eq 'n')) {
		my @sample_names = split /,/, $opt_s;
		foreach my $sample_name(@sample_names) {
			my ($sample_number, $found) = vcflines::VCF_and_sample_name_to_sample_number($file, $sample_name);
			next if($found eq 0);
			warn "Found $sample_name in $file...\n";
			my ($counts_temp, $types_found_temp) = vcflines::parse_VCF($file, 'none', 'all', 'n', $sample_name);	
			%counts = $merger->merge(\%counts, $counts_temp);
			%types_found = $merger->merge(\%types_found, $types_found_temp);
		}
	} 
	# Everything
	elsif($opt_e ne 'n') {
		# Get sample names
		my $all_sample_names = vcflines::VCF_to_sample_names($file);
		foreach my $sample_numbers(sort keys %{$all_sample_names}) {
			my $sample_name = $$all_sample_names{$sample_numbers};
			warn "Found $sample_name in $file...\n";
			my ($counts_temp, $types_found_temp) = vcflines::parse_VCF($file, 'none', 'all', 'n', $sample_name);
			%counts = $merger->merge(\%counts, $counts_temp);
			%types_found = $merger->merge(\%types_found, $types_found_temp);
		}
	}
	# Otherwise the first sample is fine
	else {
		my ($counts_temp, $types_found_temp) = vcflines::parse_VCF($file, 'none', 'all', 'n', 'WGS');
		%counts = $merger->merge(\%counts, $counts_temp);
		%types_found = $merger->merge(\%types_found, $types_found_temp);
	}
}
#print "Dumper:\n";
#print Dumper(%types_found);

# Merge phase group lengths
my %phase_group_lengths;
foreach my $file(@vcf_files) {
	warn "Calculating summary from $file...\n";
	foreach my $sample(sort keys %{$types_found{'Phase_group_lengths'}{$file}}) {
		foreach my $contig(sort keys %{$types_found{'Phase_group_lengths'}{$file}{$sample}}) {
			foreach my $length(sort keys %{$types_found{'Phase_group_lengths'}{$file}{$sample}{$contig}}) {
				my $tally = $types_found{'Phase_group_lengths'}{$file}{$sample}{$contig}{$length};
				#warn "making hash of hap length $length, tally = $tally\n";
				$phase_group_lengths{$length} += $tally;
			}
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

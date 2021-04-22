#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;

# rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 <VCF file 1> > outfile.vcf\n";
die $usage if (@ARGV != 1);
my ($vcf_file1) = @ARGV;

# Go through VCF 
open my $fh, '<', $vcf_file1 or die "Cannot open $vcf_file1 : $!\n";
warn "Reading $vcf_file1...\n";
VCF: while (my $line = <$fh>) {
   	chomp $line;
	my ($VCF_line) = vcflines::read_VCF_lines($line);

	# Print header
	if($$VCF_line{'header'} eq'Y') {
		print "$line\n";
		next VCF;
	}

	# Make everything a phased group called supercontig name
	if($$VCF_line{'0phased'}) { die "line is already phased: $line\n"; }
	if($$VCF_line{'base_type0'} eq 'heterozygous') { die "Found heterozygous position. Script intended only for haploid VCFs: $line\n"; }

	$$VCF_line{'info'} = "PHASE=$$VCF_line{'supercontig'}";

	# print
	my $new_line = join("\t", $$VCF_line{'supercontig'},$$VCF_line{'position'},$$VCF_line{'id'},$$VCF_line{'reference_VCF_format'},$$VCF_line{'consensus_VCF_format'},$$VCF_line{'cons_qual'},$$VCF_line{'filter'},$$VCF_line{'info'},$$VCF_line{'format'},$$VCF_line{'sample_info'});
	print "$new_line\n";
}
close $fh;

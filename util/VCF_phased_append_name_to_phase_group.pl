#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;

### r.farrer@exeter.ac.uk

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

	# Prepend supercontig name to phase group number
	if($$VCF_line{'0phased'}) {
		die "Phase group not found (PID0): $line\n" if(!defined $$VCF_line{'PID0'});
		my $phase_group = $$VCF_line{'PID0'};
		$$VCF_line{'PID0'} = "$$VCF_line{'supercontig'}\_$phase_group";
		my @format_parts = split /:/, $$VCF_line{'format'};
		my @sample_info_parts = split /:/, $$VCF_line{'sample_info'};
		my $new_sample_info;
		for(my $i=0; $i<scalar(@format_parts); $i++) {
			if($format_parts[$i] eq 'PID') {
				$new_sample_info .= "$$VCF_line{'PID0'}:";
			} else {
				$new_sample_info .= "$sample_info_parts[$i]:";
			}
		}
		$new_sample_info =~ s/:$//;
		$$VCF_line{'sample_info'} = $new_sample_info;
	}

	# print
	my $new_line = join("\t", $$VCF_line{'supercontig'},$$VCF_line{'position'},$$VCF_line{'id'},$$VCF_line{'reference_VCF_format'},$$VCF_line{'consensus_VCF_format'},$$VCF_line{'cons_qual'},$$VCF_line{'filter'},$$VCF_line{'info'},$$VCF_line{'format'},$$VCF_line{'sample_info'});
	print "$new_line\n";
}
close $fh;

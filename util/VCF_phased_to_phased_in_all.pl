#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;
use read_FASTA;
use genome_array;

### rfarrer@broadinsitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF's (separated by comma)> -f <reference FASTA> > outfile_PIA\n
Optional: -c\tLength cut-off for minimum haplotype to be considered [10]
          -p\tPhased in any (1), Phased in all (2) [1]\n";
our($opt_c, $opt_f, $opt_p, $opt_v);
getopt('cfpv');
die $usage unless ($opt_v && $opt_f);
if(!defined $opt_c) { $opt_c = 10; }
if(!defined $opt_p) { $opt_p = 1; }
my @files = split /,/, $opt_v;
foreach(@files) { die "file $_ not readable: $!\n" if(! -e $_); }
die "file $opt_f not readable : $!\n" if(! -e $opt_f);

# Phase group size and metrics etc.
my (%phase_group_lengths, %covered_regions, %haplotype_regions);
my $cumulative_size;
my @phase_group_lengths_array;

# Save sequences and lengths
my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($opt_f);
my ($seq_lengths) = fastafile::fasta_id_to_seq_length_hash($opt_f);

# Create array of zeros representing genome
my $genome_coverage = genomearray::make_genome_hash_array_from_seq_hash($sequences);

# Go through each VCF get info 
foreach my $file(@files) {

	# Haplotype local variables
	my $phase_info;

	# Save genomic coordinates of haplotypes
	my $supercontig;
	my $phased_positions = 0;

	open my $fh, '<', $file or die "Cannot open $file\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		next VCF1 if($$VCF_line{'next'} eq 1);

		# New supercontig
		if((!defined $supercontig) || ($supercontig ne $$VCF_line{'supercontig'})) {
			if(defined $supercontig) { 
				warn "\tfound $phased_positions phased positions.\n";
				$phased_positions = 0;
			}
			$supercontig = $$VCF_line{'supercontig'};
			warn "Processing $supercontig...\n";
		}

		# Phased only
		next VCF1 if(!defined $$VCF_line{'phased'}); 
		$phased_positions++;

		# Phase metrics
		($phase_info) = vcflines::get_phase_metrics($phase_info, $$VCF_line{'PID0'}, $$VCF_line{'position'}, $$VCF_line{'supercontig'}); 
		#warn "$line ... Previous info = $$phase_info{'previous_group'}, $$phase_info{'previous_first_position'}, $$phase_info{'previous_last_position'}\n";
		next VCF1 if(!defined $$phase_info{'save'});

		# New phase group = push haplotype into array
		my $phase_group_to_save = $$phase_info{'save'}{'group'};
		my $supercontig_of_phase_group_to_save = $$phase_info{'save'}{'supercontig'};
		my $first_phase_group_position = $$phase_info{'save'}{'first_position'};
		my $last_phase_group_position = $$phase_info{'save'}{'last_position'};
		my $length_of_last_phase_group = ($last_phase_group_position - $first_phase_group_position);	

		# New phase group = push haplotype into array
		push (@phase_group_lengths_array, $length_of_last_phase_group);
		$cumulative_size += $length_of_last_phase_group;
		$phase_group_lengths{$length_of_last_phase_group}++;
		#warn "New haplotype = $supercontig_of_phase_group_to_save $first_phase_group_position - $last_phase_group_position\n";
		$haplotype_regions{$supercontig_of_phase_group_to_save}{$first_phase_group_position}{$last_phase_group_position} = 1;

		# fill array					
		for(my $i=$first_phase_group_position; $i<=$last_phase_group_position;$i++) { 
			if($opt_p eq 1) { @{$$genome_coverage{$$VCF_line{'supercontig'}}}[$i] = '1'; }
			else { @{$$genome_coverage{$$VCF_line{'supercontig'}}}[$i]++; }
		}
		$length_of_last_phase_group = 0;
	}
	close $fh;
}
warn "Identified haplotypes on " . scalar(keys(%haplotype_regions)) . " supercontigs\n";

# Summarise and print
warn "Converting arrays into haplotype locations...\n";
my ($PIA, $sum_tally_haplotypes) = (0, 0);
my ($PIA_excluded, $sum_tally_haplotypes_excluded) = (0, 0);
my $number_of_files_needed = 1;
if($opt_p eq 2) { $number_of_files_needed = scalar(@files); }
CONTIGS: foreach my $contig(sort keys %{$genome_coverage}) {
	my @test = @{$$genome_coverage{$contig}};
	my $length = $$seq_lengths{$contig};

	my $am_i_searchin =0;
	my $start_of_no_alignment;

	for(my $i=0; $i<$length;$i++) {
		# A phased region found or continues
		if($test[$i] eq $number_of_files_needed) {
			## Found start
			if($am_i_searchin eq 0) { $start_of_no_alignment = ($i+1); }
			$am_i_searchin=1;
		}

		# End of phased region or still nothing found
		if((($opt_p eq 1) && (($test[$i] eq 0) || ($i eq ($length-1)))) ||
		   (($opt_p eq 2) && (($test[$i] < $number_of_files_needed) || ($i eq ($length-1))))) {
			$am_i_searchin=0;
			if(defined $start_of_no_alignment) {
				my $current = ($i+1);
				my $hap_length = ($current - $start_of_no_alignment);
				if($hap_length > $opt_c) {
					$PIA+=$hap_length;
					$sum_tally_haplotypes++;
					print "$contig\t$start_of_no_alignment\t$current\n";
				} else {
					$PIA_excluded+=$hap_length;
					$sum_tally_haplotypes_excluded++;
				}
			}
			undef $start_of_no_alignment;
		}
	}
}

# Summary
warn "Setting: $opt_p (1=phased in any, 2=phased in all)\n";
warn "PIA excluded (<$opt_c nt) = $PIA_excluded\n";
warn "Number of haplotypes excluded (<$opt_c nt) = $sum_tally_haplotypes_excluded\n\n";

warn "Phased in All (PIA) (nt) = $PIA\n";
warn "Number of haplotypes = $sum_tally_haplotypes\n";

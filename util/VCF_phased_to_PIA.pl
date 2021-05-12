#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;
use read_FASTA;
use genome_array;
use Data::Dumper;

### rfarrer@broadinsitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF's (separated by comma)> -f <reference FASTA>\n
Optional: -c\tLength cut-off for minimum haplotype to be considered [10]
          -p\tPhased in any (1), Phased in all (2) [1]
          -s\tSample names to restrict analysis too (separated by comma) [all]
	  -t\tPhase Tag (HaplotypeTools uses PID, Whatshap uses PS or HP etc) [PID]\n
Outputs:  -u\tSummary [opt-v-PIA-p-Opt_p-c-Opt_c-s-Opt_s.summary]
          -o\tOutput [opt-v-PIA-p-Opt_p-c-Opt_c-s-Opt_s.tab]\n";
our($opt_c, $opt_f, $opt_o, $opt_p, $opt_s, $opt_t, $opt_u, $opt_v);
getopt('cfopstuv');
die $usage unless ($opt_v && $opt_f);
if(!defined $opt_c) { $opt_c = 10; }
if(!defined $opt_p) { $opt_p = 1; }
if(!defined $opt_s) { $opt_s = 'all'; }
if(!defined $opt_t) { $opt_t = 'PID'; }
my @files = split /,/, $opt_v;
foreach(@files) { die "file $_ not readable: $!\n" if(! -e $_); }
if(scalar(@files) eq 1) {
	if(!defined $opt_u) { $opt_u = "$files[0]-PIA-p-$opt_p-c-$opt_c-s-$opt_s.summary"; }
	if(!defined $opt_o) { $opt_o = "$files[0]-PIA-p-$opt_p-c-$opt_c-s-$opt_s.tab"; }
}
else {
	if(!defined $opt_u) { $opt_u = "$files[0]-plus_other_VCFs-PIA-p-$opt_p-c-$opt_c-s-$opt_s.summary"; }
	if(!defined $opt_o) { $opt_o = "$files[0]-plus_other_VCFs-PIA-p-$opt_p-c-$opt_c-s-$opt_s.tab"; }
}
die "file $opt_f not readable : $!\n" if(! -e $opt_f);
my %samples;
my $number_of_samples_wanted = 0;
if($opt_s ne 'all') { 
	my @samples_array = split /,/, $opt_s; 
	$number_of_samples_wanted = scalar(@samples_array);
	foreach my $sample(@samples_array) {
		$samples{$sample} = 1;
	}
} 
warn "$0: Settings -v $opt_v -f $opt_f -c $opt_c -p $opt_p -s $opt_s -t $opt_t -u $opt_u -o $opt_o\n";

# Phase group size and metrics etc.
my (%phase_group_lengths, %covered_regions, %haplotype_regions);
my $sum_number_of_samples = 0;
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
	my %isolate_names;

	# Save genomic coordinates of haplotypes
	my $supercontig;
	my $phased_positions = 0;
	my $count_once_per_VCF = 0;

	open my $fh, '<', $file or die "Cannot open $file\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		if($$VCF_line{'isolate_names'}) { %isolate_names = %{$$VCF_line{'isolate_names'}}; }
		next VCF1 if($$VCF_line{'next'} eq 1);
		my $contig = $$VCF_line{'supercontig'};
		my $position = $$VCF_line{'position'};
		my $number_of_samples = $$VCF_line{'number_of_samples'};

		# How many samples in this file (count once)
		if($count_once_per_VCF eq 0) {
			$count_once_per_VCF = 1;
			$sum_number_of_samples += $number_of_samples;
		}

		# New supercontig
		if((!defined $supercontig) || ($supercontig ne $$VCF_line{'supercontig'})) {
			if(defined $supercontig) { 
				warn "\tfound $phased_positions phased positions.\n";
				$phased_positions = 0;
			}
			$supercontig = $$VCF_line{'supercontig'};
			warn "Processing $supercontig...\n";
		}

		# For each samples (if multi VCF)
		SAMPLE: for(my $i=0; $i < $number_of_samples; $i++) {

			# restrict to not all samples
			my $isolate_name = $isolate_names{$i};
			if(%samples) {
				next SAMPLE if(!defined $samples{$isolate_name});
			}
			my $sample_info_name = "sample_info$i";
			#my $phased_id = ('PID' . $i);
			my $phased_id = ($opt_t . $i);

			# Phased only
			next SAMPLE if(!defined $$VCF_line{$phased_id}); 
			next SAMPLE if($$VCF_line{$phased_id} eq '.');
			$phased_positions++;

			# Phase metrics
			($phase_info) = vcflines::get_phase_metrics($phase_info, $$VCF_line{$phased_id}, $position, $contig); 
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
	}
	close $fh;
}
warn "Identified haplotypes on " . scalar(keys(%haplotype_regions)) . " supercontigs\n";

# Output files
open my $ofh1, '>', $opt_o or die "Cannot open $opt_o : $!\n";
open my $ofh2, '>', $opt_u or die "Cannot open $opt_u : $!\n";

# Summarise and print
warn "Converting arrays into haplotype locations...\n";
my ($PIA, $sum_tally_haplotypes) = (0, 0);
my ($PIA_excluded, $sum_tally_haplotypes_excluded) = (0, 0);
my $number_of_files_needed = 1;
if($opt_p eq 2) { 
	# if sample name defined as well as asking for phased in all
	if($opt_s ne 'all') { $number_of_files_needed = $number_of_samples_wanted; }
	# otherwise everything
	else { $number_of_files_needed = $sum_number_of_samples; }
}
CONTIGS: foreach my $contig(sort keys %{$genome_coverage}) {
	my @test = @{$$genome_coverage{$contig}};
	my $length = $$seq_lengths{$contig};

	my $am_i_searchin =0;
	my $start_of_new_haplotype;

	for(my $i=0; $i<$length;$i++) {
		# A phased region found or continues
		if($test[$i] eq $number_of_files_needed) {
			## Found start
			if($am_i_searchin eq 0) { $start_of_new_haplotype = $i; }
			$am_i_searchin=1;
		}

		# End of phased region or still nothing found
		if((($opt_p eq 1) && (($test[$i] eq 0) || ($i eq ($length-1)))) ||
		   (($opt_p eq 2) && (($test[$i] < $number_of_files_needed) || ($i eq ($length-1))))) {
			$am_i_searchin=0;
			if(defined $start_of_new_haplotype) {
				my $current = ($i+1);
				my $hap_length = ($current - $start_of_new_haplotype);
				if($hap_length > $opt_c) {
					$PIA+=$hap_length;
					$sum_tally_haplotypes++;
					print $ofh1 "$contig\t$start_of_new_haplotype\t$current\n";
				} else {
					$PIA_excluded+=$hap_length;
					$sum_tally_haplotypes_excluded++;
				}
			}
			undef $start_of_new_haplotype;
		}
	}
}

# Summary
print $ofh2 "Setting: $opt_p (1=phased in any, 2=phased in all)\n";
print $ofh2 "PIA (nt) = $PIA\n";
print $ofh2 "Number of haplotypes = $sum_tally_haplotypes\n\n";
print $ofh2 "PIA excluded (<$opt_c nt) = $PIA_excluded\n";
print $ofh2 "Number of haplotypes excluded (<$opt_c nt) = $sum_tally_haplotypes_excluded\n";

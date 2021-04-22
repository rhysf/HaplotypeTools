#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;
use VCF_phase;
use File::Basename;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF-phased-Sample_Number (from VCF_phase.pl)>
Optional: -c\tCut-off percent reads supporting phase group [90]
          -m\tMinimum read depth overlapping two heterozygous SNPs [4]
          -z\tVerbose for error checking (y/n) [n]
Output:   -o\tOutput VCF [opt_v-and-assigned.tab]
          -l\tTallies for percent of reads agreeing with phases [opt_v-Tally-agree-disagree-opt_m-MD-opt_c-pc_cutoff.tab]\n";
our($opt_c, $opt_l, $opt_m, $opt_o, $opt_s, $opt_v, $opt_z);
getopt('clmosvz');
die $usage unless ($opt_v);
if(!defined $opt_c) { $opt_c = 90; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_o) { $opt_o = $opt_v . '-and-assigned.tab'; }
if(!defined $opt_l) { $opt_l = $opt_v .'-Tally-agree-disagree-' . $opt_m . '-MD-' . $opt_c . '-pc_cutoff.tab'; }
if(!defined $opt_z) { $opt_z = 'n'; }
$opt_c /= 100;
die "VCF needs to be made by VCF_phase.pl. Name will be similar to sample.vcf-contig1-0-10000-phased-0" if($opt_v !~ /-(.+)-(\d+)-(\d+)-phased-(\d+)$/);

# File name parts
my $VCF = fileparse($opt_v);
my ($contig, $start_window, $stop_window, $sample_number) = $VCF =~ /-(.+)-(\d+)-(\d+)-phased-(\d+)$/;

# Save VCF to hash
warn "Save $opt_v to memory...\n";
my $VCF_position_to_line = vcflines::VCF_to_position_to_line_hash($opt_v);

# Save heterozygous positions to array
my ($base_type_id, $sample_id, $genotype_id) = ("base_type$sample_number", "sample_info$sample_number", "GT$sample_number");
my @VCF_positions;
foreach my $positions(sort { $a <=> $b } keys %{$VCF_position_to_line}) {
	my $line = $$VCF_position_to_line{$positions};
	my $VCF_line = vcflines::read_VCF_lines($line);
	my $base_type = $$VCF_line{$base_type_id};

	# Only save heterozygous positions
	next if($base_type ne 'heterozygous');
	push(@VCF_positions, $positions);
}

# Go through VCF positions
my %saved_phase_agree_and_disagree;
my $phase_switch_genotype = 0;
my $phase_block = 0;
my $previous_VCF_entry;
warn "Validate and assign phase groups...\n";
LINES: for(my $i=0; $i<scalar(@VCF_positions); $i++) {
	my $position = $VCF_positions[$i];

	# Save VCF and phase info
	my $line = $$VCF_position_to_line{$position};
	if($opt_z eq 'y') { warn "\n\n$position) $line\n"; }
	my $VCF_line = vcflines::read_VCF_lines($line);
	my $number_of_samples = $$VCF_line{'number_of_samples'};
	my $phase_summary = vcfphase::VCF_phased_id_read_numbers_to_summary($$VCF_line{'id'});

	# Positions without or not enough read info saved (same as depth cutoff essentially - except no data for foreach loop)
	if(!defined $$phase_summary{'top_two_haplotype_depth'}) {
		if($opt_z eq 'y') { warn "top_two_haplotype_depth not defined from $line\n"; }

		# If previous VCF line is defined, remove the read numbers (unless it's been phased)
		if(defined $previous_VCF_entry) {
			my $previous_VCF_line = vcflines::read_VCF_lines($previous_VCF_entry);
			if($$previous_VCF_line{'format'} !~ m/PID/) {
				$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'remove_read_numbers', 'n/a', $sample_number, 'n/a', $opt_z);
				my $previous_position = &VCF_line_to_position($previous_VCF_entry);
				$$VCF_position_to_line{$previous_position} = $previous_VCF_entry;
			}
		}

		# Update current line
		$line = &VCF_update_phase($line, 'insufficient_read_depth', 'n/a', $sample_number, 'n/a', $opt_z);
		$$VCF_position_to_line{$position} = $line;	

		# reset and new phase group
		undef $previous_VCF_entry;
		$phase_switch_genotype = 0;
		$phase_block++; 
		next LINES;
	}

	# Depth cutoff
	if($opt_z eq 'y') { warn "Check depth...\n"; }
	my $count = 0;
	foreach my $read_count(sort { $b <=> $a } keys %{$$phase_summary{'read_count_to_phase_group'}}) {
		# only check for 2 haplotypes
		last if($count eq 2);
		$count++;

		# if either are below opt_m, make new line
		if($read_count < $opt_m) {
			if($opt_z eq 'y') { warn "read count $read_count < $opt_m\n"; }

			# If previous VCF line is defined, remove the read numbers (unless it's been phased)
			if(defined $previous_VCF_entry) {
				my $previous_VCF_line = vcflines::read_VCF_lines($previous_VCF_entry);
				if($$previous_VCF_line{'format'} !~ m/PID/) {
				
					$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'remove_read_numbers', 'n/a', $sample_number, 'n/a', $opt_z);
					my $previous_position = &VCF_line_to_position($previous_VCF_entry);
					$$VCF_position_to_line{$previous_position} = $previous_VCF_entry;
				}
			}

			# Update current line
			$line = &VCF_update_phase($line, 'insufficient_read_depth', 'n/a', $sample_number, 'n/a', $opt_z);
			$$VCF_position_to_line{$position} = $line;	

			# reset and new phase group
			undef $previous_VCF_entry;
			$phase_switch_genotype = 0;
			$phase_block++; 
			next LINES;
		}
	}

	# Save phase cutoff percent
	my $percent_phase_agree = ($$phase_summary{'top_two_haplotype_depth'} / $$phase_summary{'read_count'}) * 100;
	my $rounded_percent = sprintf("%.2f", $percent_phase_agree);
	$saved_phase_agree_and_disagree{$rounded_percent}++;
	if($opt_z eq 'y') { warn "Phase percent = $rounded_percent ($$phase_summary{'top_two_haplotype_depth'} / $$phase_summary{'read_count'}) * 100) - check cutoff...\n"; }

	# Filter on phase cutoff
	if($$phase_summary{'top_two_haplotype_depth'} < ($$phase_summary{'read_count'} * $opt_c)) {
		if($opt_z eq 'y') { warn "haplotype depth $$phase_summary{'top_two_haplotype_depth'} < ($$phase_summary{'read_count'} * $opt_c)\n"; }

		# If previous VCF line is defined, remove the read numbers (unless it's been phased)
		if(defined $previous_VCF_entry) {
			my $previous_VCF_line = vcflines::read_VCF_lines($previous_VCF_entry);
			if($$previous_VCF_line{'format'} !~ m/PID/) {
				$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'remove_read_numbers', 'n/a', $sample_number, 'n/a', $opt_z);
				my $previous_position = &VCF_line_to_position($previous_VCF_entry);
				$$VCF_position_to_line{$previous_position} = $previous_VCF_entry;
			}
		}

		# Update current line
		$line = &VCF_update_phase($line, 'base_below_cutoff', 'n/a', $sample_number, 'n/a', $opt_z);
		$$VCF_position_to_line{$position} = $line;	

		# reset and new phase group
		undef $previous_VCF_entry;
		$phase_switch_genotype = 0;
		$phase_block++; 
		next LINES;
	}

	# Save the first VCF position of the contig and move on to next VCF position
	if(!defined $previous_VCF_entry) {
		$previous_VCF_entry = $line; 
		$phase_switch_genotype = 0;
		$phase_block++;
		if($opt_z eq 'y') { warn "Saving as the first phase_block: $phase_block\n"; }
		next LINES;
	}

	# I have 2 positions. Save VCF and phase info for line 2
	if($opt_z eq 'y') { warn "Two lines found. Save previous VCF to try and phase...\n"; }
	my $previous_VCF_line = vcflines::read_VCF_lines($previous_VCF_entry);
	my $previous_position = $$previous_VCF_line{'position'};
	my $previous_id = $$previous_VCF_line{'id'};
	my $previous_phase_summary = vcfphase::VCF_phased_id_read_numbers_to_summary($previous_id);

	# Read overlap between 2 positions (PG1->PG2)
	my ($phase_match, $phase_matches) = vcfphase::VCF_phase_summaries_to_match($previous_phase_summary, $phase_summary);

	# Some consecutive hets don't have any reads overlapping - they might not be close enough
	# Also are both alleles predominantly hitting two separate alleles (e.g. some positions have a mix of reads linking ref to ref and ref to consensus - which in turn lead to a het call... but bad evidence of phase)
	if($phase_matches < 2) {
		if($opt_z eq 'y') { 
			warn "phase match < 2 (no read overlaps or conflicting info). Remove reads from previous VCF\n";
			warn Dumper($phase_match);
		}

		# Remove read numbers from previous line (unless it's been phased)
		if($$previous_VCF_line{'format'} !~ m/PID/) {
			$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'remove_read_numbers', 'n/a', $sample_number, 'n/a', $opt_z);
			$$VCF_position_to_line{$previous_position} = $previous_VCF_entry;
		}

		# reset and new phase group
		$previous_VCF_entry = $line;
		$phase_switch_genotype = 0;
		$phase_block++; 
		next LINES;
	}

	# Create new genotypes (still might need reordering based on previous switching
	my $updated_previous_GT;
	my $updated_GT;
	my $switch_gt_flag = 0;
	foreach my $PG1(keys %{$phase_match}) {
		my $PG2 = $$phase_match{$PG1};
		$updated_previous_GT .= "$PG1|";
		$updated_GT .= "$PG2|";

		# Switch genotype?
		if($PG1 ne $PG2) { $switch_gt_flag = 1; }
	}
	$updated_previous_GT =~ s/\|$//;
	$updated_GT =~ s/\|$//;
	#warn "new) $updated_previous_GT and $updated_GT\n";

	# Phaseable
	if($opt_z eq 'y') { warn "PHASEABLE!\n"; }

	# Update id, format and sample info with phase group for previous line
	$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'remove_read_numbers', 'n/a', $sample_number, 'n/a', $opt_z);
	if($opt_z eq 'y') { warn "previous VCF format = $$previous_VCF_line{'format'}\n"; }
	if($$previous_VCF_line{'format'} !~ m/PID/) {
		if($opt_z eq 'y') { warn "Phase previous line\n"; }

		# Switch GT?
		if(($phase_switch_genotype % 2) ne 0) {
			if($updated_previous_GT =~ m/(\d)[\/|\|](\d)/) {
				$updated_previous_GT = ($2 . '|' . $1);
			}
			if($opt_z eq 'y') { warn "Switch GT = $updated_previous_GT\n"; }
		}

		# Phase
		$previous_VCF_entry = &VCF_update_phase($previous_VCF_entry, 'phase', ":$contig-$sample_number-$start_window-$phase_block", $sample_number, $updated_previous_GT, $opt_z);
		$$VCF_position_to_line{$previous_position} = $previous_VCF_entry;	
	}

	# Phase switch genotypes for current line
	if($switch_gt_flag eq 1) { $phase_switch_genotype++; }

	# Update current line to include phase group (but leave read numbers in for phasing to new lines after)
	# if even don't switch, if odd do switch
	if(($phase_switch_genotype % 2) ne 0) {
		if($updated_GT =~ m/(\d)[\/|\|](\d)/) {
			$updated_GT = ($2 . '|' . $1);
		}
		if($opt_z eq 'y') { warn "Switch GT = $updated_GT\n"; }
	}

	# Phase
	$line = &VCF_update_phase($line, 'phase', ":$contig-$sample_number-$start_window-$phase_block", $sample_number, $updated_GT, $opt_z);
	if($opt_z eq 'y') { warn "save $position -> $line\n"; }
	$$VCF_position_to_line{$position} = $line;

	# Save read info for comparing with next read
	$previous_VCF_entry = $line;
	$previous_position = $position;
}

# Print phased or phase_info lines (those updated)
open my $ofh, '>', $opt_o or die "Cannot open $opt_o : $!";
VCFPOSITIONS: foreach my $positions(sort { $a <=> $b } keys %{$VCF_position_to_line}) {
	my $line = $$VCF_position_to_line{$positions};

	# save space. Only print supercontig, position, and sample
	my @bits = split /\t/, $line;
	my $supercontig = $bits[0];
	my $position = $bits[1];
	my $id = $bits[2];
	next VCFPOSITIONS if($id !~ m/phase/);
	my $sample_info = $bits[($sample_number + 9)];
	#print $ofh "$line\n";
	print $ofh "$supercontig\t$position\t$id\t$sample_info\n";
}
close $ofh;

# Print Percent Phased agree tally
#warn "Saving percent phased agree tally file...\n";
open my $ofh2, '>', $opt_l or die "Cannot open $opt_l: $!";
print $ofh2 "Percent Agree\tTally\n";
foreach my $percent(sort { $a <=> $b } keys %saved_phase_agree_and_disagree) {
	my $tally = $saved_phase_agree_and_disagree{$percent};
	print $ofh2 "$percent\t$tally\n";
}
close $ofh2;

sub VCF_line_to_position {
	my $line = $_[0];
	my $VCF_line = vcflines::read_VCF_lines($line);
	my $position = $$VCF_line{'position'};
	return $position;
}

sub VCF_update_phase {
	my ($line, $task, $phase_group, $sample_number, $new_genotype, $verbose) = @_;
	my $VCF_line = vcflines::read_VCF_lines($line);

	# Remove read numbers
	if($task eq 'remove_read_numbers') { $$VCF_line{'id'} = ''; }

	# Info = insufficient_read_depth
	if($task eq 'insufficient_read_depth') { $$VCF_line{'id'} = 'phase_info=sample' . $sample_number . '-haplotype_DP_below_cutoff;'; }

	# Info = base_below_cutoff
	if($task eq 'base_below_cutoff') { $$VCF_line{'id'} = 'phase_info=sample' . $sample_number . '-base_below_cutoff;'; }

	# Info = phase_below_cutoff
	if($task eq 'phase_below_cutoff') { $$VCF_line{'id'} = 'phase_info=sample' . $sample_number . '-phase_below_cutoff;'; }

	# Phase
	if($task eq 'phase') {
		$$VCF_line{'id'} = 'sample-' . $sample_number . '-phased;';
		$$VCF_line{'format'} .= ":PID"; # Needed so lines are valid for read_VCF

		# Update every sample (necessary so lines are valid for read_VCF, even though not all samples are printed)
		my $sample_id = "sample_info$sample_number";
		my $genotype_id = "GT$sample_number";
		my $number_of_samples = $$VCF_line{'number_of_samples'};
		for(my $i=0; $i<$number_of_samples; $i++) {
			if($sample_number eq $i) { 
				my @sample_id_parts = split /:/, $$VCF_line{$sample_id};
				my $genotype = $$VCF_line{$genotype_id};
				my $new_sample_id = '';
				foreach my $sample_id_part(@sample_id_parts) {
					if($sample_id_part eq $genotype) { 
						$sample_id_part = $new_genotype;
						#if($sample_id_part =~ m/(\d)[\/|\|](\d)/) {
						#	$sample_id_part = ($2 . '|' . $1);
						#}
					}
					$new_sample_id .= "$sample_id_part:";
				}
				$new_sample_id =~ s/:$//;
				$$VCF_line{$sample_id} = ($new_sample_id . $phase_group); 
			}
			else { $$VCF_line{"sample_info$i"} .= ":."; }
		}
	}

	$line = vcflines::VCF_line_make($VCF_line);
	if($verbose eq 'y') { warn "New line: $line\n"; }
	return ($line);
}


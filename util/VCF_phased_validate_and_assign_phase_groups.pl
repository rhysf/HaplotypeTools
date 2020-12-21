#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_VCF;
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF-phased (from VCF_and_BAM_to_phased_VCF)>
Optional: -c\tCut-off percent reads supporting phase group [90]
          -m\tMinimum read depth overlapping two heterozygous SNPs [4]
Output:   -o\tOutput VCF [opt_v-and-assigned.VCF]
	  -l\tTallies for percent of reads agreeing with phases [opt_v-Tally-agree-disagree-opt_m-MD-opt_c-pc_cutoff.tab]\n";
our($opt_c, $opt_l, $opt_m, $opt_o, $opt_v);
getopt('clmov');
die $usage unless ($opt_v);
if(!defined $opt_c) { $opt_c = 90; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_o) { $opt_o = $opt_v . '-and-assigned.VCF'; }
if(!defined $opt_l) { $opt_l = $opt_v .'-Tally-agree-disagree-' . $opt_m . '-MD-' . $opt_c . '-pc_cutoff.tab'; }
$opt_c /= 100;
die "VCF needs to be made by VCF_and_BAM_to_phased_VCF. Name will be similar to unphased.vcf-contig1-0-10000-phased" if($opt_v !~ /-(.+)-(\d+)-(\d+)-phased$/);

# Save VCF to hash
my $VCF_position_to_line = vcflines::VCF_to_position_to_line_hash($opt_v);

# File name parts
my ($contig, $start_window, $stop_window) = $opt_v =~ /-(.+)-(\d+)-(\d+)-phased$/;

# Go through each VCF position in this contig
my %saved_phase_agree_and_disagree;
my $previous_VCF_entry;
my @previous_read_numbers;
my $phase_block = 0;
VCFPOSITIONS: foreach my $positions(sort { $a <=> $b } keys %{$VCF_position_to_line}) {
	my $line = $$VCF_position_to_line{$positions};
	my $VCF_line = vcflines::read_VCF_lines($line);

	# Only sort phased data here
	if ($line !~ m/\|/) {

		# If something is saved in memory, remove read numbers
		if(defined $previous_VCF_entry) {
			my ($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($previous_VCF_entry, $phase_block, $start_window, 0, 5);
			$$VCF_position_to_line{$position} = $newline;
		}
		next VCFPOSITIONS;	
	}

	# Count reads (dis)agreeing with phase (and remove from id)
	my @read_numbers = split /\;/, $$VCF_line{'id'};
	my ($new_id, $read_count, $phase_agree, $phase_change, $other_bases) = vcfphase::VCF_phased_id_read_numbers_to_phase_pc($$VCF_line{'id'});
	$$VCF_line{'id'} = $new_id;	
	$line = &VCF_line_make($VCF_line);

	# Min depth (replace read numbers and save)
	if ($read_count < $opt_m) {
		my ($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($line, $phase_block, $start_window, 0, 4);
		$$VCF_position_to_line{$position} = $newline;
		#warn "MD < $opt_m: $newline\n";
		next VCFPOSITIONS;
	}

	# Different base voids phasing (replace read numbers and save)
	if($other_bases > ($read_count * $opt_c)) {
		my ($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($line, $phase_block, $start_window, 0, 3);
		$$VCF_position_to_line{$position} = $newline;
		#warn "Other bases $other_bases > $read_count * $opt_c: $newline\n";
		next VCFPOSITIONS;
	}

	# Check if position appears below cut-off
	if($phase_agree ne $read_count) { warn "$phase_agree/$read_count\n"; }
	my $next_entry = 0;
	if($phase_agree < ($read_count * $opt_c)) { 
		warn "phase_agree ($phase_agree) < $read_count * $opt_c\n";
		# Check for alternative phase
		if($phase_change > ($read_count * $opt_c)) { 
			warn "Phase group change for: $line\n"; 
			if($$VCF_line{'GT0'} =~ m/^0\|1/) {
				my $old_ref = $$VCF_line{'reference_VCF_format'};
				$$VCF_line{'reference_VCF_format'} = $$VCF_line{'consensus_VCF_format'};
				$$VCF_line{'consensus_VCF_format'} = $old_ref;
				$line = &VCF_line_make($VCF_line);
				warn "New line = $line\n";
			}
			elsif($$VCF_line{'GT0'} =~ m/^1\|2/) {
				my @alternative_consensus = split /,/, $$VCF_line{'consensus_VCF_format'};
				$$VCF_line{'consensus_VCF_format'} = ($alternative_consensus[1] . ',' . $alternative_consensus[0]);
				$line = &VCF_line_make($VCF_line);
				warn "New line = $line\n";
			}
			my $temp_agree = $phase_agree;
			$phase_agree = $phase_change;
			$phase_change = $temp_agree;
		}
		# Otherwise remove phase
		else {
			my ($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($line, $phase_block, $start_window, 0, 2);
			$$VCF_position_to_line{$position} = $newline;
			warn "Remove phase: $newline\n";
			$next_entry = 1;
		}
	}
	my $percent_phase_agree = 0;
	if(($phase_agree + $phase_change) eq 0) { warn "Phase agree = $phase_agree, phase change = $phase_change, over $line\n"; }
	else { $percent_phase_agree = (($phase_agree / ($phase_agree + $phase_change)) * 100); }
	my $rounded_percent = sprintf("%.2f", $percent_phase_agree);
	$saved_phase_agree_and_disagree{$rounded_percent}++;
	next VCFPOSITIONS if($next_entry eq 1);

	# Save the first VCF position of the contig and move on to next VCF position
	if(!defined $previous_VCF_entry) {
		$previous_VCF_entry = $line; 
		@previous_read_numbers = split /\;/, $$VCF_line{'id'};
		$phase_block++;
		#warn "Saving as the first phase_block: $phase_block\n";
		next VCFPOSITIONS;
	}

	# Clean up read numbers with phase block
	my $found_overlap = 0;
	LOOKFOROVERLAP: foreach my $current_read_numbers(@read_numbers) {
		foreach my $previous_read_numbers(@previous_read_numbers) {
			if($current_read_numbers eq $previous_read_numbers) {
				$found_overlap = 1;
				last LOOKFOROVERLAP;
			}
		}
	}

	# Edit previous line with phase block number and save
	my ($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($previous_VCF_entry, $phase_block, $start_window, $found_overlap, 1);
	$phase_block = $new_phase_block;
	$$VCF_position_to_line{$position} = $newline;
	#warn "1) Saving $position -> $newline\n";

	# Edit current line with phase block number in case it is the end of the phase block if overlapping with last read
	($newline, $new_phase_block, $position) = vcfphase::replace_read_number_with_phase($line, $phase_block, $start_window, $found_overlap, 1);
	$$VCF_position_to_line{$position} = $newline;
	$line = $newline;
	#warn "2) Saving $position -> $newline\n";

	# Save read info for comparing with next read
	$previous_VCF_entry = $line;
	@previous_read_numbers = split /\;/, $$VCF_line{'id'};
}

# Print new VCFs VCF_position_to_line
open my $ofh, '>', $opt_o or die "Cannot open $opt_o : $!";
VCFPOSITIONS: foreach my $positions(sort { $a <=> $b } keys %{$VCF_position_to_line}) {
	my $line = $$VCF_position_to_line{$positions};
	print $ofh "$line\n";
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

sub VCF_line_make {
	my $VCF_line = $_[0];
	my $line = join "\t", 
	$$VCF_line{'supercontig'}, 
	$$VCF_line{'position'}, 
	$$VCF_line{'id'}, 
	$$VCF_line{'reference_VCF_format'}, 
	$$VCF_line{'consensus_VCF_format'}, 
	$$VCF_line{'cons_qual'}, 
	$$VCF_line{'filter'}, 
	$$VCF_line{'info'}, 
	$$VCF_line{'format'}, 
	$$VCF_line{'sample_info0'};
	return $line;
}

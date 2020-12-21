#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -a <VCF file 1> -b <VCF file 2>\n
Optional: -o Outfile extension [HaplotypeTools_VCF1_vs_VCF2]>\n";
our($opt_a, $opt_b, $opt_o);
getopt('abo');
die $usage unless($opt_a && $opt_b);
if(!defined $opt_o) { $opt_o = 'HaplotypeTools_VCF1_vs_VCF2'; }

# Outfiles
my $of_summary = ($opt_a . '-' . $opt_o . '-Phase-Summary.tab');
my $of_overlaps = ($opt_a . '-' . $opt_o . '-Phase-Overlaps.tab');
my $of_positions = ($opt_a . '-' . $opt_o . '-Phase-Positions.tab');

open my $ofh1, '>', $of_summary or die "Cannot open $of_summary: $!\n";
open my $ofh2, '>', $of_overlaps or die "Cannot open $of_overlaps: $!\n";
open my $ofh3, '>', $of_positions or die "Cannot open $of_positions: $!\n";
print $ofh2 "Overlapping Phase Group\tContig\tPosition\tVCF1\tVCF2\tIdentity\n";
print $ofh3 "Overlapping Phase Group\tContig\tPosition\tVCF1\tVCF2\tIdentity\n";

# Save phase groups from VCFs
my ($phased_hets1) = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_a, $ofh1);
my ($phased_hets2) = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_b, $ofh1);

# Compare phase groups
my $found_in_both = 0;
my ($total_overlap, $total_same, $total_cross) = (0,0,0);
PHASEGROUP1: foreach my $phase_group1(sort {$a<=>$b} keys %{$phased_hets1}) {
	#warn "$phase_group1 / " . scalar(keys(%phased_hets1)) . "\n";
	PHASEGROUP2: foreach my $phase_group2(sort {$a<=>$b} keys %{$phased_hets2}) {
		foreach my $contig(sort keys %{$$phased_hets1{$phase_group1}}) {
			next PHASEGROUP2 if(!exists $$phased_hets2{$phase_group2}{$contig});

			# Same contig in the 2 different phase groups
			my $found = 0;
			POSITIONS1: foreach my $positions(sort keys %{$$phased_hets1{$phase_group1}{$contig}}) {
				next POSITIONS1 if(!exists $$phased_hets2{$phase_group2}{$contig}{$positions});
				$found++;

				# Look for overlaps of at least 2
				if($found eq 2) {
					my ($lines, $cross_overs, $tally_overlap, $tally_same, $tally_cross) = &make_printable_lines($found_in_both, $contig, $$phased_hets1{$phase_group1}{$contig}, $$phased_hets2{$phase_group2}{$contig});
					$total_overlap += $tally_overlap;
					$total_same += $tally_same;
					$total_cross += $tally_cross;
					print $ofh2 "$lines";
					print $ofh3 "$cross_overs";
					$found_in_both++;
					next PHASEGROUP1;
				}
			}
		}
	}
}
close $ofh2;
close $ofh3;

# Calculate percents
my ($same_percent, $cross_percent) = (0, 0);
if(($total_overlap ne 0) && ($total_same ne 0)) { $same_percent = (($total_same / ($total_same + $total_cross)) * 100); }
if(($total_overlap ne 0) && ($total_cross ne 0)) { $cross_percent = (($total_cross / ($total_same + $total_cross)) * 100); }
my $total = ($total_same + $total_cross);

# Summary
warn "Overlapping phase groups = $found_in_both\n";
warn "Overlapping phased positions = $total\n";
warn "Overlapping phased positions (SAME) = $total_same / $total = $same_percent\%\n";
warn "Overlapping phased positions (CROSS) = $total_cross / $total = $cross_percent\%\n";
print $ofh1 "Overlapping phase groups = $found_in_both\n";
print $ofh1 "Overlapping phased positions = $total_overlap\n";
print $ofh1 "Overlapping phased positions (SAME) = $total_same / $total = $same_percent\%\n";
print $ofh1 "Overlapping phased positions (CROSS) = $total_cross / $total = $cross_percent\%\n";
close $ofh1;

sub make_printable_lines {
	my ($overlapping_group_number, $contig, $PHG1HASH, $PHG2HASH) = @_;
	my ($group_of_lines, $cross_over_lines) = ('', '');
	my $overlap_found;
	my ($tally_overlap, $tally_same, $tally_cross) = (0,0,0);
	
	POS: foreach my $positions(sort {$a<=>$b} keys %{$PHG1HASH}) {
		my $PHG1DATA = $$PHG1HASH{$positions};
		if(!exists $$PHG2HASH{$positions}) { 
			$group_of_lines .= "$overlapping_group_number\t$contig\t$positions\t$PHG1DATA\t-\t-\n"; 
			next POS;
		}

		my $PHG2DATA = $$PHG2HASH{$positions};
		$tally_overlap++;
		my $same = '?';

		# After defining.
		if(defined $overlap_found) {
			# Started in same order
			if($overlap_found eq 1) {
				my @parts1 = split /\|/, $PHG1DATA;
				my @parts2 = split /\|/, $PHG2DATA;
				if($PHG1DATA eq $PHG2DATA) { 
					#warn "same\n";
					$same = 'Y'; 
					$tally_same++;
				}
				elsif(($parts1[0] eq $parts2[1]) && ($parts1[1] eq $parts2[0])) { 
					#warn "cross\n";
					$same = 'Cross-over'; 
					$tally_cross++;

					# Change order based on cross-over
					$overlap_found = 2;
				}
				else { $same = 'Different-bases'; }
			}

			# Started in different order
			else {
				#warn "not same starting order.. so check this the other way round\n";
				my @parts1 = split /\|/, $PHG1DATA;
				my @parts2 = split /\|/, $PHG2DATA;
				if(($parts1[0] eq $parts2[1]) && ($parts1[1] eq $parts2[0])) { 
					#warn "same\n";
					$same = 'Y'; 
				}
				elsif($PHG1DATA eq $PHG2DATA) { 
					#warn "cross\n";
					$same = 'Cross-over'; 
					$tally_cross++;

					# Change order based on cross-over
					$overlap_found = 1;
				}
				else { $same = 'Different-bases'; }
			}
		}
			
		# Define order
		if(!defined $overlap_found) {
			#warn "init\n";
			# Same order to start
			if($PHG1DATA eq $PHG2DATA) { 
				#warn "same\n";
				$same = 'Y'; 
				$tally_same++;
				$overlap_found = 1;
			}
			# Not same order to start
			else {
				#warn "not same order\n";
				my @parts1 = split /\|/, $PHG1DATA;
				my @parts2 = split /\|/, $PHG2DATA;
				if(($parts1[0] eq $parts2[1]) && ($parts1[1] eq $parts2[0])) {
					$same = 'Y';
					$tally_same++;
					$overlap_found = 2;
				}
				else {
					#warn "Error Init: Different bases?: file 1) $PHG1DATA and file 2) $PHG2DATA. Overlap not defined\n";	
					$same = 'Different-bases(un-init)';
				}
			}
		}
		if($same =~ m/Cross/) { $cross_over_lines .= "$overlapping_group_number\t$contig\t$positions\t$PHG1DATA\t$PHG2DATA\t$same\n"; }
		$group_of_lines .= "$overlapping_group_number\t$contig\t$positions\t$PHG1DATA\t$PHG2DATA\t$same\n";
	}
	return ($group_of_lines, $cross_over_lines, $tally_overlap, $tally_same, $tally_cross);
}

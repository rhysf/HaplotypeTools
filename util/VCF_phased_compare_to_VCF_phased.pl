#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;
use VCF_phase;

### r.farrer@exeter.ac.uk

# Opening commands 
my $usage = "Usage: perl $0 -a <VCF file 1> -b <VCF file 2> || -a <VCF file 1> -c <Sample ID 1> -d <Sample ID 2>\n
Optional: -o Outfile extension [HaplotypeTools_VCF1_vs_VCF2]>\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_o);
getopt('abcdo');
die $usage unless(($opt_a && $opt_b) || ($opt_a && $opt_c && $opt_d));
if(!defined $opt_o) { $opt_o = 'HaplotypeTools_VCF1_vs_VCF2'; }

# Outfiles
my $of_summary = "$opt_a-$opt_o-Phase-Summary.tab";
my $of_overlaps = "$opt_a-$opt_o-Phase-Overlaps.tab";
my $of_positions = "$opt_a-$opt_o-Phase-Positions.tab";
my $of_refined = "$opt_a-$opt_o-Phase-Positions-Refined.tab";

open my $ofh1, '>', $of_summary or die "Cannot open $of_summary: $!\n";
open my $ofh2, '>', $of_overlaps or die "Cannot open $of_overlaps: $!\n";
open my $ofh3, '>', $of_positions or die "Cannot open $of_positions: $!\n";
open my $ofh4, '>', $of_refined or die "Cannot open $of_refined: $!\n";

print $ofh2 "Overlapping_Phase_Group\tContig\tPosition\tVCF1\tVCF2\tIdentity\n";
print $ofh3 "Overlapping_Phase_Group\tContig\tPosition\tVCF1\tVCF2\tIdentity\n";
print $ofh4 "Overlapping_Phase_Group\tContig\tPosition\tCategory\n";

# Load phased hets
my ($phased_hets1, $phased_hets2);
# Two separate VCFs
if($opt_a && $opt_b) {
	$phased_hets1 = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_a, $ofh1, 0);
	$phased_hets2 = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_b, $ofh1, 0); 
}
# One VCF
else {
	# Find sample numbers
	my ($sample_number1, $sample_found1) = vcflines::VCF_and_sample_name_to_sample_number($opt_a, $opt_c);
	my ($sample_number2, $sample_found2) = vcflines::VCF_and_sample_name_to_sample_number($opt_a, $opt_d);
	$phased_hets1 = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_a, $ofh1, $sample_number1);
	$phased_hets2 = vcfphase::VCF_phased_to_phase_group_contig_pos_to_bases($opt_a, $ofh1, $sample_number2); 
}

# Compare phase groups between samples
warn "Comparing phase groups...\n";
my ($found_in_both, $total_overlap, $total_same, $total_cross) = (0,0,0,0);
my %crossovers_by_group;

PHASEGROUP1: foreach my $pg1(sort keys %{$phased_hets1}) {
	#warn "$pg1 / " . scalar(keys(%phased_hets1)) . "\n";
	PHASEGROUP2: foreach my $pg2(sort keys %{$phased_hets2}) {
		foreach my $contig(sort keys %{$$phased_hets1{$pg1}}) {
			next PHASEGROUP2 if(!exists $$phased_hets2{$pg2}{$contig});

			# Collect shared positions
			my @shared_positions = grep {
				exists $$phased_hets2{$pg2}{$contig}{$_}
			} keys %{$$phased_hets1{$pg1}{$contig}};

			# Require at least two heterozygous overlaps to compare phase
			next unless scalar(@shared_positions) >= 2;

			# Call the crossover analysis for this shared block
			my ($lines, $cross_overs, $overlap, $same, $cross) =
				make_printable_lines(
				    $found_in_both,
				    $contig,
				    $$phased_hets1{$pg1}{$contig},
				    $$phased_hets2{$pg2}{$contig}
			);

			# Update totals
			$total_overlap += $overlap;
			$total_same    += $same;
			$total_cross   += $cross;

			# Output
			print $ofh2 $lines;
			print $ofh3 $cross_overs;

			# Record crossovers for classification
            foreach my $line (split /\n/, $cross_overs) {
                my @cols = split /\t/, $line;
                push @{ $crossovers_by_group{$pg1}{$contig} }, $cols[2] if defined $cols[2];
            }

			$found_in_both++;
		}
	}
}
close $ofh2;
close $ofh3;

# Classify crossovers
my $distance_threshold = 1000;  # <1kb = likely gene conversion
my $support_threshold  = 3;     # ≥3 phased hets per side
my ($count_meiosis, $count_conversion) = (0, 0);

foreach my $pg (sort keys %crossovers_by_group) {
    foreach my $contig (sort keys %{ $crossovers_by_group{$pg} }) {
        my @positions = sort { $a <=> $b } @{ $crossovers_by_group{$pg}{$contig} };

        for (my $i = 0; $i < @positions; $i++) {
            my $pos = $positions[$i];

            # nearest distances
            my $dist_prev = ($i > 0) ? $pos - $positions[$i - 1] : undef;
            my $dist_next = ($i < $#positions) ? $positions[$i + 1] - $pos : undef;
            my $min_dist;
            $min_dist = $dist_prev if defined $dist_prev;
            $min_dist = $dist_next if defined $dist_next && (!defined $min_dist || $dist_next < $min_dist);

            # default classification
            my $class = "Cross-over_unresolved";

            if (defined $min_dist && $min_dist < $distance_threshold) {
                $class = "Crossover_likely_gene_conversion";
                $count_conversion++;
            } else {
                my $left_support  = count_phased_sites($pg, $contig, $pos, -1, $phased_hets1, $phased_hets2);
                my $right_support = count_phased_sites($pg, $contig, $pos,  1, $phased_hets1, $phased_hets2);
                if ($left_support >= $support_threshold && $right_support >= $support_threshold) {
                    $class = "Crossover_meiosis_candidate";
                    $count_meiosis++;
                }
            }
            print $ofh4 join("\t", $pg, $contig, $pos, $class), "\n";
        }
    }
}
close $ofh4;

# Percentages
my $total_crossovers = $total_cross || 1; # avoid division by zero
my $meiosis_percent    = sprintf("%.2f", ($count_meiosis / $total_crossovers) * 100);
my $conversion_percent = sprintf("%.2f", ($count_conversion / $total_crossovers) * 100);

# Original stats
my ($same_percent, $cross_percent) = (0, 0);
if ($total_overlap && $total_same) { $same_percent = sprintf("%.2f", ($total_same / ($total_same + $total_cross)) * 100); }
if ($total_overlap && $total_cross) { $cross_percent = sprintf("%.2f", ($total_cross / ($total_same + $total_cross)) * 100); }
my $total = ($total_same + $total_cross);

# Write summary to screen and file
warn  "Overlapping phase groups = $found_in_both\n";
warn  "Overlapping phased positions = $total\n";
warn  "Overlapping phased positions (SAME) = $total_same / $total = $same_percent%\n";
warn  "Overlapping phased positions (CROSS) = $total_cross / $total = $cross_percent%\n";
warn  "Crossovers (Likely Gene Conversion) = $count_conversion ($conversion_percent%)\n";
warn  "Crossovers (Meiosis Candidates)     = $count_meiosis ($meiosis_percent%)\n";

print $ofh1 "Overlapping phase groups = $found_in_both\n";
print $ofh1 "Overlapping phased positions = $total_overlap\n";
print $ofh1 "Overlapping phased positions (SAME) = $total_same / $total = $same_percent%\n";
print $ofh1 "Overlapping phased positions (CROSS) = $total_cross / $total = $cross_percent%\n";
print $ofh1 "Crossovers (Likely Gene Conversion) = $count_conversion ($conversion_percent%)\n";
print $ofh1 "Crossovers (Meiosis Candidates)     = $count_meiosis ($meiosis_percent%)\n";
close $ofh1;

# Helper: count phased sites nearby
sub count_phased_sites {
    my ($pg, $contig, $pos, $dir, $hets1, $hets2) = @_;
    my $count = 0;
    my $limit = 10000; # search window

    foreach my $ref ($hets1, $hets2) {
        next unless exists $$ref{$pg}{$contig};
        my $sites = $$ref{$pg}{$contig};

        foreach my $p (sort { $a <=> $b } keys %$sites) {
            next if ($dir == -1 && $p >= $pos);
            next if ($dir ==  1 && $p <= $pos);
            next if abs($p - $pos) > $limit;
            $count++;
            last if $count >= 3;
        }
    }
    return $count;
}

sub make_printable_lines {
	my ($group_num, $contig, $PHG1, $PHG2) = @_;
	my ($lines, $cross_lines) = ('', '');
	my $phase_state;  # 1 = same, 2 = inverted
	my ($overlap, $same, $cross) = (0,0,0);
	my $overlap_found;

	POS: foreach my $pos(sort {$a<=>$b} keys %{$PHG1}) {
		my $v1 = $$PHG1{$pos};
		unless (exists $$PHG2{$pos}) {
			$lines .= "$group_num\t$contig\t$pos\t$v1\t-\t-\n"; 
			next POS;
		}

		my $v2 = $$PHG2{$pos};
		$overlap++;

		my @a1 = split /\|/, $v1;
        my @a2 = split /\|/, $v2;
        my $status = 'Different-bases';

        # Determine phase relationship
        if (!defined $phase_state) {
        	if ($v1 eq $v2) { $phase_state = 1; $status = 'Y'; $same++; }
            elsif ($a1[0] eq $a2[1] && $a1[1] eq $a2[0]) { $phase_state = 2; $status = 'Y'; $same++; }
            else { $status = 'Different-bases(un-init)'; }
        } else {
        	if ($phase_state == 1) {
        		if ($v1 eq $v2) { $status = 'Y'; $same++; }
                elsif ($a1[0] eq $a2[1] && $a1[1] eq $a2[0]) { $status = 'Cross-over'; $cross++; $phase_state = 2; }
        	} else { # inverted phase
        		if ($a1[0] eq $a2[1] && $a1[1] eq $a2[0]) { $status = 'Y'; $same++; }
                elsif ($v1 eq $v2) { $status = 'Cross-over'; $cross++; $phase_state = 1; }
        	}
        }

		$lines .= "$group_num\t$contig\t$pos\t$v1\t$v2\t$status\n";
		$cross_lines .= "$group_num\t$contig\t$pos\t$v1\t$v2\t$same\n" if $status =~ /Cross/;
	}
	return ($lines, $cross_lines, $overlap, $same, $cross);
}
warn "Finished.\n";
#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_FASTA;
use read_Tab;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: $0 -a <PIA_1 FASTA> -b <PIA_2 FASTA> -c <IRMS HET.fasta> -d <IRMS HET.details>\n
Optional: -p\tprint to stdout (s) or (f) [s]
          -o\tFile to print to if opt_p=f [opt_a-and-2.accuracy]
	  -s\tSummary file to print to if opt_p=f [opt_a-and-2.accuracy.summary]
	  -v\tVerbose for testing (y/n) [n]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_o, $opt_p, $opt_s, $opt_v);
getopt('abcdopsv');
die $usage unless ($opt_a && $opt_b && $opt_c && $opt_d);
if(!defined $opt_p) { $opt_p = 's'; }
if(!defined $opt_o) { $opt_o = ($opt_a . '-and-2.accuracy'); }
if(!defined $opt_s) { $opt_s = ($opt_a . '-and-2.accuracy.summary'); }
if(!defined $opt_v) { $opt_v = 'n'; }
warn "$0: Settings -a $opt_a -b $opt_b -c $opt_c -d $opt_d -p $opt_p -o $opt_o -s $opt_s\n";

# Save sequences
my $fasta1 = fastafile::fasta_to_struct($opt_a);
my $fasta2 = fastafile::fasta_to_struct($opt_b);
my $fasta3 = fastafile::fasta_to_struct($opt_c);

# Save tab file
my $irms_details = tabfile::save_columns_to_column_hash($opt_d, 0, 1, 2, 3);

# Count sum IRMS included
my $total_number_variants = `wc -l < $opt_d`;
chomp($total_number_variants);
$total_number_variants -= 1;
die "Not sufficient IRMS lines found: $total_number_variants\n" if($total_number_variants < 1);

my $header = "Haplotype_identity\tNum_difference_in_hap\ttrue_positive_phase\tfalse_positive_phase\tfalse_positive_variant_call\tfalse_negative\tswitch_errors\tswitch_error_rate\n";
my $ofh;
my $ofh2;
if($opt_p eq 's') { print $header; }
else {
	open $ofh, '>', $opt_o or die "Cannot open $opt_o : $!\n";
	open $ofh2, '>', $opt_s or die "Cannot open $opt_s : $!\n";
	print $ofh $header;
}

# For each pair or sequences, how closely do they match the contig and homologous contig (e.g. supercont1.1 or supercont1.1_homologous)
my %adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives; 
my $positions_seen_total = 0;
my $switch_errors_total = 0;
my $count = 0;
ID: foreach my $id(sort keys %{$$fasta1{'seq'}}) {
	if($opt_v eq 'y') { warn "$id\n"; }
	die "Nothing found for $id in $opt_b\n" if(!defined $$fasta2{'seq'}{$id});
	my $seq1 = $$fasta1{'seq'}{$id};
	my $seq2 = $$fasta2{'seq'}{$id};

	# Save contig and haplotype coords from ID
	my @id_parts = split /\_/, $id;
	my $contig = '';
	for(my $i=0; $i<(scalar(@id_parts) - 2); $i++) {
		$contig .= ($id_parts[$i] . "_");
	}
	$contig =~ s/\_$//;
	my $start = $id_parts[scalar(@id_parts) - 2];
	my $end = $id_parts[scalar(@id_parts) - 1];

	# Save truth set
	my $contig_hom = ($contig . "_homologous");
	die "Cannot find $contig in $opt_c\n" if(!defined $$fasta3{'seq'}{$contig});
	die "Cannot find $contig_hom in $opt_c\n" if(!defined $$fasta3{'seq'}{$contig_hom});
	my $seq3 = $$fasta3{'seq'}{$contig};
	my $seq4 = $$fasta3{'seq'}{$contig_hom};
	my $truth_set_start = ($start - 1);
	my $seq3_haplotype = substr $seq3, $truth_set_start, ($end - $truth_set_start);
	my $seq4_haplotype = substr $seq4, $truth_set_start, ($end - $truth_set_start);

	# Summary
	if($opt_v eq 'y') { 
		warn "contig = $contig\n";
		warn "start = $start\n";
		warn "end = $end\n";
		#warn "seq 1 =\n$seq1\n";
		#warn "seq 2 =\n$seq2\n";
		#warn "truth set seq3 =\n$seq3_haplotype\n";
		#warn "truth set seq4 =\n$seq4_haplotype\n";
	}

	# Pull IRMS details from those regions (pos -> nt -> introduced het)
	my ($IRMS_positions, $num_mutations) = &save_number_mutations_in_region($irms_details, $contig, $truth_set_start, $end);
	my ($IRMS_positions_homologous, $num_mutations_homologous) = &save_number_mutations_in_region($irms_details, $contig_hom, $truth_set_start, $end);
	my %IRMS_positions = (%{$IRMS_positions}, %{$IRMS_positions_homologous});
	if($opt_v eq 'y') {
		#print "irms positions\n";
		#print Dumper($IRMS_positions);
		#warn "num mutation = $num_mutations\n";
		
		#print "irms homologous positions\n";
		#print Dumper($IRMS_positions_homologous);
		#warn "num mutation homolgous = $num_mutations_homologous\n";

		#print "\ncombined\n";
		#print Dumper(%IRMS_positions);
	}
	next ID if(($num_mutations + $num_mutations_homologous) eq 0);
	my $num_differences_in_hap = ($num_mutations + $num_mutations_homologous);

	# Standard accuracy
	my $TP = 0;
	my $FN = 0;
	my $false_positive = 0;
	my $FP_phase = 0;

	# Switch errors
	my $switch_errors = 0;
	my $switch_for_errors;

	# QAN50
	my $positions_seen_for_sub_haplotype_with_no_switch_errors_total_for_haplotype = 0;
	my $positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
	my $span_start;
	my $span_last_seen;
	my $TP_in_sub_haplotype = 0;

	# Go through each position that has been changed
	my $hap_match;
	POSITIONS: foreach my $positions(sort { $a <=> $b } keys %IRMS_positions) {
		$positions_seen_for_sub_haplotype_with_no_switch_errors_total_for_haplotype++;
		$positions_seen_for_sub_haplotype_with_no_switch_errors++;
		$positions_seen_total++;

		# Pull out bases of all four haplotypes (2 real, 2 truth)
		my $hap_position = ($positions - $truth_set_start);
		my $hap_base1 = substr $seq1, $hap_position, 1;
		my $hap_base2 = substr $seq2, $hap_position, 1;
		my $hap_base3 = substr $seq3_haplotype, $hap_position, 1;
		my $hap_base4 = substr $seq4_haplotype, $hap_position, 1;
		if(!defined $span_start) { $span_start = $hap_position; }
		$span_last_seen = $hap_position;

		# Summary
		if($opt_v eq 'y') {
			warn "\nIRMS position = $positions\n";
			warn "Hap position = $hap_position\n";
			warn "Hap bases = \n";
			warn "\t1 = $hap_base1\n"; 
			warn "\t2 = $hap_base2\n";
			warn "\t3 = $hap_base3\n"; 
			warn "\t4 = $hap_base4\n";
		}

		# Define which haplotypes (1 and 2) match either normal or homologous simulated haplotypes/truth set
		if(!defined $hap_match) {
			# 1 = 3, 2 = 4
			if(($hap_base1 eq $hap_base3) && ($hap_base2 eq $hap_base4)) {
				if($opt_v eq 'y') { warn "Define hap: 1 = 3, 2 = 4\n"; }
				#warn "Define hap: 1 = 3, 2 = 4\n";
				$hap_match = 1;
				$switch_for_errors = 1;
				$TP++;
				$TP_in_sub_haplotype++;
			}
			# 1 = 4, 2 = 3
			elsif(($hap_base1 eq $hap_base4) && ($hap_base2 eq $hap_base3)) {
				if($opt_v eq 'y') { warn "Define hap: 1 = 4, 2 = 3\n"; }
				#warn "Define hap: 1 = 4, 2 = 3\n";
				$hap_match = 2;
				$switch_for_errors = 2;
				$TP++;
				$TP_in_sub_haplotype++;
			}
			# FN
			elsif($hap_base1 eq $hap_base2) {
				$FN++;
			}
			# FP / incorrect base (correct position) from VCF
			else {
				$false_positive++;
			}
			next POSITIONS;
		}

		# Compare
		
		##### Switch Error Rate and QAN50
		if(($hap_base1 eq $hap_base3) && ($hap_base2 eq $hap_base4)) {
			#warn "1=3, 2=4 (hap match $hap_match, switch for errors $switch_for_errors)\n";

			# correct (but previously wrong!)
			if(($hap_match eq 1) && ($switch_for_errors eq 2)) {
				#warn "switch error: correct (but previously wrong)\n";
				$switch_errors++;
				$switch_for_errors = 1;

				# QAN50
				my $span_length = ($hap_position - $span_start);
				my $adj_span_length = $span_length * ($TP_in_sub_haplotype / $positions_seen_for_sub_haplotype_with_no_switch_errors);
				$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_in_sub_haplotype}++;
				# reset
				$positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
				$span_start = $hap_position;
				$TP_in_sub_haplotype = 0;
			}

			# wrong (but previously correct!)
			elsif(($hap_match eq 2) && ($switch_for_errors eq 2)) {
				#warn "switch error: wrong (and previously correct)\n";
				$switch_errors++;
				$switch_for_errors = 1;
				
				# QAN50
				my $span_length = ($hap_position - $span_start);
				my $adj_span_length = $span_length * ($TP_in_sub_haplotype / $positions_seen_for_sub_haplotype_with_no_switch_errors);
				$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_in_sub_haplotype}++;
				# reset
				$positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
				$span_start = $hap_position;
				$TP_in_sub_haplotype = 0;
			}

			else {}
		}
		elsif(($hap_base1 eq $hap_base4) && ($hap_base2 eq $hap_base3)) { 
			#warn "1=4, 2=3 (hap match $hap_match, switch for errors $switch_for_errors)\n";

			# correct (but previously wrong!)
			if(($hap_match eq 2) && ($switch_for_errors eq 1)) {
				#warn "switch error: correct (but previously wrong)\n";
				$switch_errors++;
				$switch_for_errors = 2;
				
				# QAN50
				my $span_length = ($hap_position - $span_start);
				my $adj_span_length = $span_length * ($TP_in_sub_haplotype / $positions_seen_for_sub_haplotype_with_no_switch_errors);
				$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_in_sub_haplotype}++;
				# reset
				$positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
				$span_start = $hap_position;
				$TP_in_sub_haplotype = 0;
			}

			# wrong (but previously correct!)
			if(($hap_match eq 1) && ($switch_for_errors eq 1)) {
				#warn "switch error: wrong (and previously correct)\n";
				$switch_errors++;
				$switch_for_errors = 2;
				
				# QAN50
				my $span_length = ($hap_position - $span_start);
				my $adj_span_length = $span_length * ($TP_in_sub_haplotype / $positions_seen_for_sub_haplotype_with_no_switch_errors);
				$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_in_sub_haplotype}++;
				# reset
				$positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
				$span_start = $hap_position;
				$TP_in_sub_haplotype = 0;
			}
		}

		# and these for switch errors (although cannot change switch_for_errors as its neither hap_match1 or hap_match2)
		# FN/Homozygous
		elsif($hap_base1 eq $hap_base2) { }

		# FP / incorrect base (correct position) from VCF
		else {
			# e.g. T ne T or C && A ne T or C
			$switch_errors++;

			# QAN50
			my $span_length = ($hap_position - $span_start);
			my $adj_span_length = $span_length * ($TP / $positions_seen_for_sub_haplotype_with_no_switch_errors);
			$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP}++;
			# reset
			$positions_seen_for_sub_haplotype_with_no_switch_errors = 0;
			$span_start = $hap_position;
		}


		##### True positives, False positives, False negatives

		# TP / correct base
		if(($hap_base1 eq $hap_base3) && ($hap_base2 eq $hap_base4)) { 
			if($opt_v eq 'y') { warn "1 = 3, 2 = 4\n"; }
			if($hap_match eq 1) { 
				$TP++; 
				$TP_in_sub_haplotype++;
			} 
			else { 
				if($opt_v eq 'y') { warn "false_positive_phase\n"; }
				$FP_phase++; 
			}
		}
		elsif(($hap_base1 eq $hap_base4) && ($hap_base2 eq $hap_base3)) { 
			if($opt_v eq 'y') { warn "1 = 4, 2 = 3\n"; }
			if($hap_match eq 2) { 
				$TP++; 
				$TP_in_sub_haplotype++;
			} 
			else { 
				if($opt_v eq 'y') { warn "false_positive_phase\n"; }
				$FP_phase++; 
			}
		}

		# FN/Homozygous
		elsif($hap_base1 eq $hap_base2) {
			if($opt_v eq 'y') { warn "Homozygous: FN\n"; }
			$FN++;
		}

		# FP / incorrect base (correct position) from VCF
		else {
			# e.g. T ne T or C && A ne T or C
			if($opt_v eq 'y') { warn "$hap_base1 ne $hap_base3 or $hap_base4 && $hap_base2 ne $hap_base3 or $hap_base4\n"; }
			$false_positive++;
		}
		#warn "TP = $TP, $FP = $FP, correct base wrong hap = $FP_phase\n";
	}

	# QAN50 for either the last block with no switch errors, or the full block, if no switch error found
	my $span_length = ($span_last_seen - $span_start);
	if($span_length > 1) {
		my $adj_span_length = $span_length * ($TP_in_sub_haplotype / $positions_seen_for_sub_haplotype_with_no_switch_errors);
		$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_in_sub_haplotype}++;
	}

	# Switch error rate (-1 because the first position cannot be a switch error)
	my $switch_error_rate = 0;
	if(($switch_errors > 0) && ($positions_seen_for_sub_haplotype_with_no_switch_errors_total_for_haplotype > 2)) {
		$switch_error_rate = ($switch_errors / ($positions_seen_for_sub_haplotype_with_no_switch_errors_total_for_haplotype - 1)); 
	}
	$switch_errors_total += $switch_errors;

	# No hap match - means all mutations are not introduced. Assign hap_match=1 as it doesn't matter either way
	if(!defined $hap_match) { $hap_match = 1; }

	# Compare pairwise to identify FP (those non IRMS)
	my ($FP1, $FP2) = (0, 0);
	if($hap_match eq 1) {
		#warn "1-3\n";
		#my ($p1_3_length, $p1_3_same, $p1_3_gaps, $p1_3_dif, $p1_3_avg) = &compare_sequences($seq1, $seq3_haplotype, $id, "$contig $start - $end", \%IRMS_positions, $truth_set_start); 
		$FP1 = &compare_sequences($seq1, $seq3_haplotype, $id, "$contig $start - $end", \%IRMS_positions, $truth_set_start); 
		#warn "2-4\n";
		#my ($p2_4_length, $p2_4_same, $p2_4_gaps, $p2_4_dif, $p2_4_avg) = &compare_sequences($seq2, $seq4_haplotype, $id, "$contig_hom $start - $end", \%IRMS_positions, $truth_set_start); 
		$FP2 = &compare_sequences($seq2, $seq4_haplotype, $id, "$contig_hom $start - $end", \%IRMS_positions, $truth_set_start); 
	}
	else {
		#warn "1-4\n";
		#my ($p1_4_length, $p1_4_same, $p1_4_gaps, $p1_4_dif, $p1_4_avg) = &compare_sequences($seq1, $seq4_haplotype, $id, "$contig_hom $start - $end", \%IRMS_positions, $truth_set_start); 
		$FP1 = &compare_sequences($seq1, $seq4_haplotype, $id, "$contig_hom $start - $end", \%IRMS_positions, $truth_set_start); 
		#warn "2-3\n";
		#my ($p2_3_length, $p2_3_same, $p2_3_gaps, $p2_3_dif, $p2_3_avg) = &compare_sequences($seq2, $seq3_haplotype, $id, "$contig $start - $end", \%IRMS_positions, $truth_set_start); 
		$FP2 = &compare_sequences($seq2, $seq3_haplotype, $id, "$contig $start - $end", \%IRMS_positions, $truth_set_start); 
	}
	my $FP = ($FP1 + $FP2 + $false_positive);

	# Summary
	my $haplotype_accuracy_line = "$id\t$num_differences_in_hap\t$TP\t$FP_phase\t$FP\t$FN\t$switch_errors\t$switch_error_rate\n";
	if($opt_p eq 's') { print $haplotype_accuracy_line; }
	else { print $ofh $haplotype_accuracy_line; }

	$count++;
}

# Summary stats

# switch error rate
my $switch_error_rate = 0;
if(($switch_errors_total > 0) && ($positions_seen_total > 2)) {
	$switch_error_rate = ($switch_errors_total / ($positions_seen_total - 1)); 
}
if($opt_p eq 's') { warn "Switch Error Rate (overall) = $switch_error_rate\n"; }
else { print $ofh2 "Switch Error Rate (overall) = $switch_error_rate\n"; }

# QAN50
my $accumulative_adj_span_length = 0;
my $accumulative_TP_count = 0;
QAN50: foreach my $adj_span_length(sort { $b <=> $a } keys %adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives) {
	#warn "big to small adj span length: $adj_span_length\n";
	foreach my $TP_count(keys %{$adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}}) {
		my $tally = $adjusted_span_for_sub_blocks_with_no_switch_error_to_true_positives{$adj_span_length}{$TP_count};
		for(my $i=0; $i<$tally; $i++) {
			$accumulative_adj_span_length += $adj_span_length;
			$accumulative_TP_count += $TP_count;

			# > 1/2 total number of variants?
			if($accumulative_TP_count > ($total_number_variants / 2)) {
				if($opt_p eq 's') { warn "QAN50 = $adj_span_length\n"; }
				else { print $ofh2 "QAN50 = $adj_span_length\n"; }
				last QAN50;
			}
		}
	}
}

sub compare_sequences {
	my ($sequence1, $sequence2, $id1, $id2, $irms_pos, $truth_set_start) = @_;
	if(length($sequence1) ne length($sequence2)) { die "Sequences $id1 and $id2 not the same length... bad comparison.\n"; }
	my @s1 = split //, $sequence1;
	my @s2 = split //, $sequence2;
	my $length = length($sequence1);

	#print "irms pos = \n";
	#print Dumper($irms_pos);

	my ($same, $gaps, $different) = (0,0,0);
	my $FP = 0;
	for(my $i=0; $i<scalar(@s1); $i++) {
		my $base_one = $s1[$i];
		my $base_two = $s2[$i];
		if($base_one eq $base_two) { $same++; }
		elsif(($base_one eq '-') || ($base_two eq '-')) { $gaps++; }
		else { 
			#warn "different at pos $i ($base_one and $base_two)\n";
			my $hap_position = ($truth_set_start + $i);
			#warn "hap position = $hap_position\n";
			if(defined $$irms_pos{$hap_position}) { $different++; }
			else { $FP++; }
		}
	}
	return $FP;
	#my $average_gaps_pc = ($gaps / $length)*100;
	#return ($length, $same, $gaps, $different, $average_gaps_pc);
}

sub save_number_mutations_in_region {
	my ($irms_hash, $contig, $start, $end) = @_;
	my $num_mutations = 0;
	my %positions;
	foreach my $position(sort { $a <=> $b } keys %{$$irms_hash{$contig}}) {
		next if($position < $start);
		last if($position > $end);
		
		foreach my $nt(keys %{$$irms_hash{$contig}{$position}}) {
			my $irms = $$irms_hash{$contig}{$position}{$nt};
			#warn "$contig $position $nt $irms\n";
			$positions{$position}{$nt} = $irms;
		}
		$num_mutations++;
	}
	return (\%positions, $num_mutations);
}


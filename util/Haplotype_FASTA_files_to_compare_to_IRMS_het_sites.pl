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
Optional:\n";
our($opt_a, $opt_b, $opt_c, $opt_d);
getopt('abcd');
die $usage unless ($opt_a && $opt_b && $opt_c && $opt_d);

# Save sequences
my $fasta1 = fastafile::fasta_to_struct($opt_a);
my $fasta2 = fastafile::fasta_to_struct($opt_b);
my $fasta3 = fastafile::fasta_to_struct($opt_c);

# Save tab file
my $irms_details = tabfile::save_columns_to_column_hash($opt_d, 0, 1, 2, 3);

# For each pair or sequences, how closely do they match the contig and homologous contig (e.g. supercont1.1 or supercont1.1_homologous)
my $count = 0;
print "id\tTP\tcorrect_position_not_match_truth_set\tcorrect_base_wrong_phase\tFP\n";
ID: foreach my $id(sort keys %{$$fasta1{'seq'}}) {
	#warn "$id\n";
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

	# Pull out truth set
	my $contig_hom = ($contig . "_homologous");
	die "Cannot find $contig in $opt_c\n" if(!defined $$fasta3{'seq'}{$contig});
	die "Cannot find $contig_hom in $opt_c\n" if(!defined $$fasta3{'seq'}{$contig_hom});
	my $seq3 = $$fasta3{'seq'}{$contig};
	my $seq4 = $$fasta3{'seq'}{$contig_hom};
	my $truth_set_start = ($start - 1);
	my $seq3_haplotype = substr $seq3, $truth_set_start, ($end - $truth_set_start);
	my $seq4_haplotype = substr $seq4, $truth_set_start, ($end - $truth_set_start);

	#warn "contig = $contig\n";
	#warn "start = $start\n";
	#warn "end = $end\n";
	#warn "seq 1 =\n$seq1\n";
	#warn "seq 2 =\n$seq2\n";
	#warn "truth set seq3 =\n$seq3_haplotype\n";
	#warn "truth set seq4 =\n$seq4_haplotype\n";

	# Pull IRMS details from those regions (pos -> nt -> introduced het)
	my ($IRMS_positions, $num_mutations) = &save_number_mutations_in_region($irms_details, $contig, $truth_set_start, $end);
	#print "irms positions\n";
	#print Dumper($IRMS_positions);
	#warn "num mutation = $num_mutations\n";

	my ($IRMS_positions_homologous, $num_mutations_homologous) = &save_number_mutations_in_region($irms_details, $contig_hom, $truth_set_start, $end);
	#print "irms hom positions\n";
	#print Dumper($IRMS_positions_homologous);
	#warn "num mutation homolgous = $num_mutations_homologous\n";
	next ID if(($num_mutations + $num_mutations_homologous) eq 0);

	my %IRMS_positions = (%{$IRMS_positions}, %{$IRMS_positions_homologous});
	#print "\ncombined\n";
	#print Dumper(%IRMS_positions);

	# Go through each position that has been changed
	my $TP = 0;
	my $correct_position_not_match_truth_set = 0;
	my $hap_match;
	my $correct_base_wrong_haplotype = 0;
	foreach my $positions(sort { $a <=> $b } keys %IRMS_positions) {

		# Pull out bases of all four haplotypes (2 real, 2 truth)
		#warn "\nIRMS position = $positions\n";
		my $hap_position = ($positions - $truth_set_start);
		#warn "Hap position = $hap_position\n";
		my $hap_base1 = substr $seq1, $hap_position, 1;
		my $hap_base2 = substr $seq2, $hap_position, 1;
		my $hap_base3 = substr $seq3_haplotype, $hap_position, 1;
		my $hap_base4 = substr $seq4_haplotype, $hap_position, 1;
		#warn "Hap bases = \n";
		#warn "\t1 = $hap_base1\n"; 
		#warn "\t2 = $hap_base2\n";
		#warn "\t3 = $hap_base3\n"; 
		#warn "\t4 = $hap_base4\n";

		# Compare

		# TP / correct base
		if(($hap_base1 eq $hap_base3) && ($hap_base2 eq $hap_base4)) { 
			#warn "1 = 3, 2 = 4\n";
			if(!defined $hap_match) {
				$hap_match = 1;
				$TP++;
			} 
			elsif($hap_match eq 1) { $TP++; } 
			else { 
				#warn "correct_base_wrong_haplotype\n";
				$correct_base_wrong_haplotype++; 
			}
		}
		elsif(($hap_base1 eq $hap_base4) && ($hap_base2 eq $hap_base3)) { 
			#warn "1 = 4, 2 = 3\n";
			if(!defined $hap_match) {
				$hap_match = 2;
				$TP++;
			}
			elsif($hap_match eq 2) { $TP++; } 
			else { 
				#warn "correct_base_wrong_haplotype\n";
				$correct_base_wrong_haplotype++; 
			}
		}

		# FP / incorrect base
		else {
			#warn "$hap_base1 ne $hap_base3 or $hap_base4 ?\n";
			$correct_position_not_match_truth_set++;
		}
		#warn "TP = $TP, $FP = $FP, correct base wrong hap = $correct_base_wrong_haplotype\n";
	}

	# No hap match - means all mutations are not introduced. Assign hap_match=1 as it doesn't matter either way
	if(!defined $hap_match) {
		
		#warn "Error: No hap match?\n";
		#warn "contig = $contig\n";
		#warn "start = $start\n";
		#warn "end = $end\n";
		#warn "seq 1 = $seq1\n";
		#warn "seq 2 = $seq2\n";
		#warn "truth set seq3 = $seq3_haplotype\n";
		#warn "truth set seq4 = $seq4_haplotype\n";
		#print Dumper(%IRMS_positions);
		
		#die "end here\n";
		$hap_match = 1;
	}

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
	my $FP = ($FP1 + $FP2);

	# Summary
	print "$id\t$TP\t$correct_position_not_match_truth_set\t$correct_base_wrong_haplotype\t$FP\n";

	$count++;
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


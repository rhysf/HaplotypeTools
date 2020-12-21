#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: $0 -f <FASTA> > output_summary\n
Optional: -o\tOutput (p=pairwise within opt_f, s=summary of pairwise within opt_f, a=additional FASTA: 1:1 order) [p]
          -g\tExclude if gaps percent is greater than this [100]
	  -a\tAdditional FASTA in same order []\n";
our($opt_f, $opt_o, $opt_g, $opt_a);
getopt('foga');
die $usage if(!defined $opt_f);
if(!defined $opt_o) { $opt_o = 'p'; }
if(!defined $opt_g) { $opt_g = 100; }
die "opt_a must be specified if opt_o eq a\n" if(($opt_o eq 'a') && (! $opt_a)); 

# Save sequences
my $fasta1 = fastafile::fasta_to_struct($opt_f);

# Compare to additional FASTA (Coded to need same ID)
if($opt_o eq 'a') {
	my $fasta2 = fastafile::fasta_to_struct($opt_a);
	print "ids\tLength\tSame\tDifferent\tGaps\n";
	foreach my $id(keys %{$$fasta1{'seq'}}) {
		die "$id in $opt_f not found in $opt_o\n" if(!defined $$fasta2{'seq'}{$id});
		my $sequence1 = $$fasta1{'seq'}{$id};
		my $sequence2 = $$fasta2{'seq'}{$id};
		my ($length, $same, $gaps, $different, $average_gaps_pc) = &compare_sequences($sequence1, $sequence2, $id, $id);
		next SEQ2 if($average_gaps_pc > $opt_g);
		# Default to pairwse
		print "$id\t$length\t$same\t$different\t$gaps\n";
	}
	die "End of 1:1 comparisons\n";
}

# Headers for pairwise
if($opt_o eq 'p') { print "id1\tid2\tLength\tSame\tDifferent\tGaps\n"; }
my ($total_length, $total_same, $max, $min) = (0,0,0,100);
my ($min_id, $max_id) = ('','');

# Pairwise comparison
my %seen;
SEQ1: foreach my $id1(keys %{$$fasta1{'seq'}}) {
	SEQ2: foreach my $id2(keys %{$$fasta1{'seq'}}) {
		next SEQ2 if($id1 eq $id2);
		next SEQ2 if(defined $seen{$id1}{$id2});
		next SEQ2 if(defined $seen{$id2}{$id1});
		$seen{$id1}{$id2} = 1;
		$seen{$id2}{$id1} = 1;

		#warn "Comparing $id1 to $id2...\n";
		my $sequence1 = $$fasta1{'seq'}{$id1};
		my $sequence2 = $$fasta1{'seq'}{$id2};
		my ($length, $same, $gaps, $different, $average_gaps_pc) = &compare_sequences($sequence1, $sequence2, $id1, $id2);
		next SEQ2 if($average_gaps_pc > $opt_g);

		# pairwise summary
		if($opt_o eq 'p') { print "$id1\t$id2\t$length\t$same\t$different\t$gaps\n"; }
		# end of file summary
		if($opt_o eq 's') {
			$total_length += $length;
			$total_same += $same;
			my $average = ($same / $length)*100;
			if($average > $max) { 
				$max = $average; 
				$max_id = "$id1-$id2";
			}
			if($average < $min) {
				$min = $average;
				$min_id = "$id1-$id2";
			}
		}
	}
}

# Summary
if($opt_o eq 's') {
	my $mean = ($total_same / $total_length)*100;
	print "Average = $mean\n";
	print "Max = $max ($max_id)\n";
	print "Min = $min ($min_id)\n";
}

sub compare_sequences {
	my ($sequence1, $sequence2, $id1, $id2) = @_;
	if(length($sequence1) ne length($sequence2)) { die "Sequences $id1 and $id2 not the same length... bad comparison.\n"; }
	my @s1 = split //, $sequence1;
	my @s2 = split //, $sequence2;
	my $length = length($sequence1);

	my ($same, $gaps, $different) = (0,0,0);
	for(my $i=0; $i<scalar(@s1); $i++) {
		my $base_one = $s1[$i];
		my $base_two = $s2[$i];
		if($base_one eq $base_two) { $same++; }
		elsif(($base_one eq '-') || ($base_two eq '-')) { $gaps++; }
		else { $different++; }
	}
	my $average_gaps_pc = ($gaps / $length)*100;
	return ($length, $same, $gaps, $different, $average_gaps_pc);
}

#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_FASTA;
use read_Tab;
use VCF_phase;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF file> -l <PIA file (contig tab start tab stop)> -r <Reference FASTA>\n
Optional: -p\tPrinting option (o=outfile, s=split to opt_l_1 and _2) [s]
          -e\tExclude printing if haplotypes are identical (y/n) [n]
	  -m\tmin length [10]
	  -i\tInclude indels (y/n) [n]\n";
our($opt_v, $opt_l, $opt_r, $opt_p, $opt_e, $opt_m, $opt_i);
getopt('vlrpemi');
die $usage unless (($opt_v) && ($opt_l) && (($opt_r)));
if(!defined $opt_p) { $opt_p = 's'; }
if(!defined $opt_e) { $opt_e = 'n'; }
if(!defined $opt_m) { $opt_m = 10; }
if(!defined $opt_i) { $opt_i = 'n'; }

# Save reference FASTA
my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($opt_r);

# Save haplotypes
my $haplotypes = tabfile::save_columns_to_column_hash($opt_l, 0, 1, 2);

# Save phased variants from VCF
my $polymorphisms = vcfphase::VCF_phased_to_contig_pos_haps_to_variant($opt_v, $opt_i);

# Outfiles
my ($ofh1, $ofh2);
if($opt_p ne 'o') { 
	my $out1 = $opt_l . "_1.fasta";
	my $out2 = $opt_l . "_2.fasta";
	open $ofh1, '>', $out1 or die "Cannot open $out1";
	open $ofh2, '>', $out2 or die "Cannot open $out2";
}

# Go through haplotypes
my ($worked, $haps_seen, $haps_same) = (0, 0, 0);
foreach my $supercontig(sort keys %{$haplotypes}) {
	#warn "Haplotypes -> supercontig = $supercontig\n";
	HAP: foreach my $start(sort keys %{$$haplotypes{$supercontig}}) {
		my $stop = $$haplotypes{$supercontig}{$start};
		#warn "$start -> $stop\n";

		# Save initial haplotype sequence from FASTA
		#warn "\nRefs = $ref_seq ($group)\n";
   		my $ref_seq = substr $$sequences{$supercontig}, ($start - 1), (($stop - $start) + 1);
   		my ($hap1, $hap2) = ($ref_seq, $ref_seq);
   		$haps_seen++;

		# Run through all phased positions on this supercontig, and check if it falls within this haplotype
   		VARIANTPOS: foreach my $variant_pos(keys %{$$polymorphisms{$supercontig}}) {
			#warn "Looking for $variant_pos in $supercontig $start - $stop\n";
			next VARIANTPOS unless(($variant_pos >= $start) && ($variant_pos <= $stop));
			#warn "FOUND.\n";

			# Check the reference bases from both hapltoypes and reference FASTA are concordant
			my $test_ref = substr $$sequences{$supercontig}, ($variant_pos - 1), 1;
   			my $test_base_in_test_seq = substr $ref_seq, ($variant_pos - $start), 1;
   			if($test_ref ne $test_base_in_test_seq) { 
   				warn "$test_ref != $test_base_in_test_seq at pos $variant_pos of \n"; 
   				next VARIANTPOS;
   			}
   			else { $worked++; }

   			# Replace bases in haplotypes
   			my $variant1 = $$polymorphisms{$supercontig}{$variant_pos}{'hap1'};
   			my $variant2 = $$polymorphisms{$supercontig}{$variant_pos}{'hap2'};
			$hap1 = &replace_bases_in_haplotype($hap1, $variant_pos, $start, $stop, length($ref_seq), $variant1);
			$hap2 = &replace_bases_in_haplotype($hap2, $variant_pos, $start, $stop, length($ref_seq), $variant2);
   		}

		# dont print if same and specified
		if($hap1 eq $hap2) {
			$haps_same++;
			next HAP if($opt_e eq 'y');
		}
		# dont print if length < $opt_m
		next HAP if((length($hap1) < $opt_m) || (length($hap2) < $opt_m));

		if($opt_p eq 'o') { 
			print ">hap_1_$supercontig\_$start\_$stop\n$hap1\n";
   			print ">hap_2_$supercontig\_$start\_$stop\n$hap2\n";
		} else {
			print $ofh1 ">$supercontig\_$start\_$stop\n$hap1\n";
			print $ofh2 ">$supercontig\_$start\_$stop\n$hap2\n";
		}
	}
}

# Outfiles
if($opt_p ne 'o') { 
	close $ofh1;
	close $ofh2;
}

# Summary
warn "$worked worked according to finding correct ref. base in sequences and haps\n";
warn "$haps_seen haps seen X2\n";
warn "$haps_same are both the same\n";

sub replace_bases_in_haplotype {
	my ($hap, $variant_pos, $start, $stop, $length, $variant) = @_;
	my $start_seq = substr $hap, 0, ($variant_pos - $start); 
	my $end_seq = '';
   	if(length($start_seq) eq $length) { warn "Same length - no end part\n"; }
   	else { $end_seq   = substr $hap, (($variant_pos - $start) + 1), $stop; }
	$hap = ($start_seq . $variant . $end_seq);
	return $hap;
}
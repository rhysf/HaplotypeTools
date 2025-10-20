#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF> -r <reference FASTA>\n
Optional: -u\tRestrict to sample name []
          -i\tInclude homozygous indels (y/n) [n]
          -h\tInclude bi-allelic heterozygous positions if ref base not included (y/n) [y]
          -a\tFor heterozygous positions, use ambiguity code (a) or first position for haploid consensus (f) [f]
	  -n\tCharacter for ambiguous [N]
	  -s\tRestrict to only this supercontig [n]
	  -t\tRestrict to only this isolate [n]
Notes: Prints to opt_v-isolate-name-consensus.fasta\n";
our($opt_a, $opt_h, $opt_i, $opt_n, $opt_r, $opt_s, $opt_t, $opt_u, $opt_v);
getopt('ahinrstuv');
die $usage unless (($opt_v) && ($opt_r));
if(!defined $opt_a) { $opt_a = 'f'; }
if(!defined $opt_i) { $opt_i = 'n'; }
if(!defined $opt_h) { $opt_h = 'y'; }
if(!defined $opt_n) { $opt_n = 'N'; }
if(!defined $opt_s) { $opt_s = 'n'; }
if(!defined $opt_t) { $opt_t = 'n'; }
die "opt_a needs to be a or f: $opt_a\n" if ($opt_a !~ m/[af]/);
die "opt_i needs to be y or n: $opt_i\n" if ($opt_i !~ m/[yn]/);
die "opt_h needs to be y or n: $opt_h\n" if ($opt_h !~ m/[yn]/);

# Go through VCF saving only polymorphic sites ($variants{$isolate_name}{$sc}{$pos}{$ref} = $cons;)
my $variants = &save_polymorphic_sites_from_VCF($opt_v, $opt_s, $opt_u); 

# Save consensus reference FASTA
my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($opt_r);
if($opt_s ne 'n') {
	warn "Restricting analysis to $opt_s...\n";
	foreach my $id(keys %{$sequences}) {
		delete $$sequences{$id} if($id ne $opt_s);
	}
}

# Create consensus reference FASTA
foreach my $isolate_name(sort keys %{$variants}) {

	# Restrict to isolate?
	if($opt_t ne 'n') {
		next if($opt_t ne $isolate_name);
	}

	# Go through FASTA sequence changing variant sites
	my $sequences_consensus = $sequences;
	CONTIG: foreach my $id(sort keys %{$sequences}) {
		my @nts = split //, $$sequences_consensus{$id};
		warn scalar(keys(%{$$variants{$isolate_name}{$id}})) . " variants in contig $id for $isolate_name...\n";	
		foreach my $variant_positions(keys %{$$variants{$isolate_name}{$id}}) {
			foreach my $variant_ref(keys %{$$variants{$isolate_name}{$id}{$variant_positions}}) {
				my $variant_consensus = $$variants{$isolate_name}{$id}{$variant_positions}{$variant_ref};
				$nts[($variant_positions - 1)] = $variant_consensus;
			}
			$$sequences_consensus{$id} = join('',@nts);
		}
	}

	# Print modified genome
	#fastafile::fasta_hash_print($sequences, $opt_r, 'fasta', 'dna', $order, 'n');
	my $output_name = "$opt_v-$isolate_name-s-$opt_s-i-$opt_i-n-$opt_n-consensus.fasta";
	if($opt_s eq 'n') { fastafile::fasta_hash_print_simple_outfile($sequences_consensus, $output_name); }
	else { fastafile::fasta_hash_print_simple_outfile($sequences_consensus, $output_name); }
}

sub save_polymorphic_sites_from_VCF {
	my ($input, $restrict_to_supercontig_or_n, $restrict_to_sample) = @_;

	my (%variants, %isolate_names);
	my ($snps, $hets, $ambiguous) = (0,0,0);
	warn "Reading $input...\n";
	open my $fh, '<', $input or die "Cannot open $input\n";
	VCF: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		if($$VCF_line{'isolate_names'}) { %isolate_names = %{$$VCF_line{'isolate_names'}}; }

		# Restrict to genome or contig
		next VCF if($$VCF_line{'next'} eq 1);
		my ($sc, $pos, $ref, $cons) = ($$VCF_line{'supercontig'}, $$VCF_line{'position'}, $$VCF_line{'reference_VCF_format'}, $$VCF_line{'consensus_VCF_format'});
		if($restrict_to_supercontig_or_n ne 'n') { next unless($restrict_to_supercontig_or_n eq $sc); }
		if(defined $restrict_to_sample) {
			die "$restrict_to_sample not found in $input\n" if(!defined $isolate_names{$restrict_to_sample});
		}

		# Multi VCF stuff
		my $number_of_isolates = scalar(keys(%isolate_names));
		ISOLATES: for(my $i=0; $i<$number_of_isolates; $i++) {
			my $isolate_name = $isolate_names{$i};
			if(defined $restrict_to_sample) { next ISOLATES unless($isolate_name eq $restrict_to_sample); }
			my $base1 = $$VCF_line{($i . 'base1')};
			my $base_type = $$VCF_line{('base_type' . $i)};
			my $amb_char = $$VCF_line{('amb_char' . $i)};

			# Ignore ref
			next VCF if($base_type eq 'reference'); 

			# Indels
			if($opt_i eq 'y') {
				if($base_type eq 'insertion') { $variants{$isolate_name}{$sc}{$pos}{$ref} = $cons; }
				if($base_type eq 'deletion') {
					my $deletion_length = (length($ref) - length($cons));
					for(my $i=($pos + 1); $i<=($pos + $deletion_length); $i++) { $variants{$isolate_name}{$sc}{$i}{$ref} = ''; }
				}
			}
			# Homozygous SNPs
			if($base_type eq 'snp') {
				$snps++;
				$variants{$isolate_name}{$sc}{$pos}{$ref} = $base1;
			}

			# Heterozygous positions 
			if($base_type eq 'heterozygous') {
				$hets++;

				# take ambiguity code
				if($opt_a eq 'a') { $variants{$isolate_name}{$sc}{$pos}{$ref} = $amb_char; } 
				else { $variants{$isolate_name}{$sc}{$pos}{$ref} = $base1; }
			}

			# ambiguous positions
			if($base_type eq 'ambiguous') {
				$ambiguous++;
				$variants{$isolate_name}{$sc}{$pos}{$ref} = $opt_n;
			}
		}
	}
	close $fh;
	warn "Restrict to supercontig or all (n): $restrict_to_supercontig_or_n\n";
	warn "Found in $input: $hets heterozygous sites, $ambiguous ambiguous sites, and $snps SNPs across " . scalar(keys(%isolate_names)) . " isolates\n";
	return \%variants;
}

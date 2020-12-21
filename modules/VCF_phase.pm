package vcfphase;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;

### rfarrer@broadinstitute.org

sub return_lines_from_BioDBHTS_array {
	my ($alignments) = $_[0];
	my $lines;
	SAMSEQS: for my $a(@{$alignments}) {
		my $start = $a->start;
		my $end = $a->end;
		my $query_dna = $a->query->dna; 
		my $query_start = $a->query->start;     
		my $query_end   = $a->query->end;

		# Ignore reads specifying/aligning over indels
		my $read_length = ($query_end - $query_start);
		my $ref_seq_length = ($end - $start);
		next SAMSEQS if ($read_length ne $ref_seq_length);

		# Save
		my $line = join "\t", $start, $end, $query_dna, $query_start, $query_end;
		$lines .= "$line\n";
	}
	return $lines;
}

sub phase_reads {
	my ($alignment_file, $VCF_file, $contig) = @_;

	# Save VCF to memory (contig -> position -> line)
	warn "phase_reads_new: Finding phase for sites in $VCF_file using $alignment_file...\n";
	my $VCF_hash = vcflines::VCF_to_position_to_line_hash($VCF_file); 

	# Second phase it
	my $read_count = 0;
	open my $fh, '<', $alignment_file or die "Cannot open $alignment_file : $!\n";
	SAMSEQS: while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($start, $end, $query_dna, $query_start, $query_end) = @bits;

		# Find reads overlapping 2 or more VCF positions 
		my $saved_vcf_positions = '';
		my $count_overlapping_vcf_positions = 0;
		READPOSITIONS: for(my $i=$start; $i<=$end; $i++) {
			if(defined $$VCF_hash{$i}) {
				$count_overlapping_vcf_positions++;
				$saved_vcf_positions .= ($$VCF_hash{$i} . "\n");
			}
		}
		next SAMSEQS if ($count_overlapping_vcf_positions < 2);

		# Go through each VCF Position associated with this read and 
		# (1) phase or if alredy phased (2) check
		my $phase_position;
		my @VCF_positions = split /\n/, $saved_vcf_positions;
		VCFPOSITIONS: foreach my $line(@VCF_positions) {
			my ($VCF_line) = vcflines::read_VCF_lines($line);
			my $base_type = $$VCF_line{'base_type0'};
			my $position = $$VCF_line{'position'};
			my $ref_base = $$VCF_line{'reference_VCF_format'};
			my $consensus = $$VCF_line{'consensus_VCF_format'};
			my $id = $$VCF_line{'id'};
			my $cons_qual = $$VCF_line{'cons_qual'};
			my $filter = $$VCF_line{'filter'};
			my $info = $$VCF_line{'info'};
			my $format = $$VCF_line{'format'};
			my $sample = $$VCF_line{'sample_info0'};
			my $genotype = $$VCF_line{'GT0'};

			# Save read count
			$id .= "$read_count;";

			# Only phase heterozygous here
			next VCFPOSITIONS if($base_type ne 'heterozygous');

			# Pull out the polymorphic base
			my $polymorphic_position = ($position - $start);
			my $base = substr $query_dna, $polymorphic_position, 1;

			# Check if new phase position is the same as old. 
			# If different, switch ref and alt allele (if not already phased from previous read)
			# was $$VCF_line{'alleles'} -> $$VCF_line{'GT0'}
			if(defined $phase_position) {
				#warn "Phase position has already been defined as $phase_position\n"; 
				($id, $ref_base, $consensus) = &switch_ref_and_consensus_if_not_in_phase($genotype, $base, $ref_base, $consensus, $id, $phase_position);
			}

			# Define a phase position for the allele
			if(!defined $phase_position) {
				#warn "defining a phase position for a read over VCF: $_ \n";	
				($id, $phase_position) = &define_phase_position($genotype, $base, $ref_base, $consensus, $id)
			}

			# Replace the phasing
			#$genotype =~ s/\//|/g;
			$sample =~ s/\//|/g;

			# Replace the entry
			my $newline = join "\t", $contig, $position, $id, $ref_base, $consensus, $cons_qual, $filter, $info, $format, $sample;
			$$VCF_hash{$position} = $newline;
		}
		$read_count++;
	}
	#print Dumper($VCF_hash);
	return ($VCF_hash);
}

# Phase position is defined. Now switch ref or consensus based on haplotype info
sub switch_ref_and_consensus_if_not_in_phase {
	my ($genotype, $base, $ref_base, $consensus, $id, $phase_position) = @_;

	if($genotype =~ m/0(\/\|)1/) {
		if((($base eq $ref_base) && ($phase_position eq 1)) || (($base eq $consensus) && ($phase_position eq 0))) { 
			if($genotype =~ m/\//) {
				# Switch REF and ALT 
				my $old_ref = $ref_base;
				$ref_base = $consensus;
				$consensus = $old_ref;
			}
			elsif($genotype =~ m/\|/) { $id =~ s/;$/-PP;/; }
		}
		elsif(($base ne $ref_base) && ($base ne $consensus)) {
			$id =~ s/;$/-NT-/;
			$id .= ($base . ';');
		}
	}
	elsif($genotype =~ m/1(\/\|)2/) {
		# Switch bases if unphased
		my @alternative_consensus = split /,/, $consensus;
		if((($base eq $alternative_consensus[0]) && ($phase_position eq 1)) || (($base eq $alternative_consensus[1]) && ($phase_position eq 0))) { 
			if($genotype =~ m/\//) {
				# Switch ALT Consensus
				$consensus = ($alternative_consensus[1] . ',' . $alternative_consensus[0]);
			}
			elsif($genotype =~ m/\|/) { $id =~ s/;$/-PP;/; }
		}
		elsif(($base ne $alternative_consensus[0]) && ($base ne $alternative_consensus[1])) {
			$id =~ s/;$/-NT-/;
			$id .= ($base . ';');
		}
	}
	return ($id, $ref_base, $consensus);
}

# Phase position (0=base is ref, 1=base is consensus)
sub define_phase_position {
	my ($genotype, $base, $ref_base, $consensus, $id) = @_;

	my $phase_position;
	if($genotype =~ m/0(\/\|)1/) {
		if($base eq $ref_base) { $phase_position = 0; }
		elsif($base eq $consensus) { $phase_position = 1; }
		else { 
			$id =~ s/;$/-NT-/;
			$id .= ($base . ';');
		}
	}
	elsif($genotype =~ m/1(\/\|)2/) {
		my @alternative_consensus = split /,/, $consensus;
		if($base eq $alternative_consensus[0]) { $phase_position = 0; }
		elsif($base eq $alternative_consensus[1]) { $phase_position = 1; }
		else { 
			$id =~ s/;$/-NT-/;
			$id .= ($base . ';');
		}
	}
	return ($id, $phase_position);
}

sub replace_read_number_with_phase {
	my ($VCFline, $phase_block, $start_window, $increment, $phase_change) = @_;	
	my @VCF_parts = split /\t/, $VCFline;
	my ($contig, $position, $id, $ref_base, $consensus, $cons_qual, $filter, $info, $format, $sample1) = @VCF_parts;

	# No overlap
	$id = '.';
	#if($phase_change eq 0) { $id = '.'; }
	if($phase_change eq 0) { }

	# End of or found overlap
	elsif($phase_change eq 1) {
		if($increment eq 0) { 
			# Has not already been confirmed to overlap
			#if($format !~ m/PID/) { $id = '.'; }
			$phase_block++; 
		}
		if($increment eq 1) {
			if($format !~ m/PID/) {
				$format .= ":PID";
				$sample1 .= ":$contig-$start_window-$phase_block;";
			}
		}
	}
	# Phase below cutoff percent
	elsif($phase_change eq 2) {
		if($info !~ m/PHASEINFO/) { $info .= ',PHASEINFO=phase_below_cutoff;'; }
		$sample1 =~ s/\|/\//;
	}
	# Base below cutoff percent
	elsif($phase_change eq 3) {
		if($info !~ m/PHASEINFO/) { $info .= ',PHASEINFO=base_below_cutoff;'; }
		$sample1 =~ s/\|/\//;
	}
	# Depth below cutoff percent
	elsif($phase_change eq 4) {
		if($info !~ m/PHASEINFO/) { $info .= ',PHASEINFO=depth_confident_reads_below_cutoff;'; }
		$sample1 =~ s/\|/\//;
	}
	# No overlap
	elsif($phase_change eq 5) {
		if($format !~ m/PID/) { 
			$info .= ',PHASEINFO=No_het_overlap;';
			$sample1 =~ s/\|/\//;
		}
	}
	
	my $newline = join "\t", $contig, $position, $id, $ref_base, $consensus, $cons_qual, $filter, $info, $format, $sample1;
	return ($newline, $phase_block, $position);
}

sub VCF_phased_id_read_numbers_to_phase_pc {
	my $id = $_[0];
	my @read_numbers = split /\;/, $id;
	my ($read_count, $phase_agree, $phase_change) = (0, 0, 0);
	my ($A, $C, $T, $G) = (0, 0, 0, 0);
	READNUMS: foreach my $read_nums(@read_numbers) {

		# Phase agree
		if($read_nums =~ m/^(\d+)$/) { $phase_agree++; }

		# Different phase
		elsif($read_nums =~ m/(\d+)-(PP)/) {
			$read_numbers[$read_count] = $1;
			$phase_change++;
		}

		# Different nucleotide
		elsif($read_nums =~ m/(\d+)-(NT-)([ACTGactg])/) {
			my $rn = $1;
			my $nt = $3;
			if($nt =~ m/A|a/) { $A++; }
			elsif($nt =~ m/C|c/) { $C++; }
			elsif($nt =~ m/T|t/) { $T++; }
			elsif($nt =~ m/G|g/) { $G++; }
			$read_numbers[$read_count] = $rn;
		}
		$read_count++;
	}
	my $new_id = join(";",@read_numbers);
	$new_id .= ';';
	my $other_bases = ($A + $C + $T + $G);
	return ($new_id, $read_count, $phase_agree, $phase_change, $other_bases);
}

sub save_two_bases {
	my ($VCF_line) = $_[0];
	my ($polymorphism1, $polymorphism2);
   	if($$VCF_line{'base_type0'} eq 'reference') {
		$polymorphism1 = $$VCF_line{'reference_VCF_format'};
		$polymorphism2 = $$VCF_line{'reference_VCF_format'};
	}
	elsif($$VCF_line{'base_type0'} eq 'snp') {
		$polymorphism1 = $$VCF_line{'consensus_VCF_format'};
		$polymorphism2 = $$VCF_line{'consensus_VCF_format'};
	}
	elsif($$VCF_line{'base_type0'} =~ m/heterozygous|het_deletion|het_insertion|insertion|deletion/) {
		$polymorphism1 = $$VCF_line{'0base1'};
		$polymorphism2 = $$VCF_line{'0base2'};
	} 
	else { 
		print Dumper($VCF_line);
		die "what is this line. Did not anticipate phasing $$VCF_line{'base_type0'}\n"; 
	}
	if((!defined $polymorphism1) || (!defined $polymorphism2)) { die "Something strange over $$VCF_line{'base_type0'}\n"; }
	return ($polymorphism1, $polymorphism2);
}


sub VCF_phased_to_phase_group_contig_pos_to_bases {
	my ($file, $ofh) = @_;
	my %phased_hets;
	my $count = 0;
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		next VCF1 if($$VCF_line{'next'} eq 1);   	
		next VCF1 if(!defined $$VCF_line{'phased'});
		my ($reference, $consensus, $saved_data) = ($$VCF_line{'reference_VCF_format'}, $$VCF_line{'consensus_VCF_format'}, '-');

		# for phased homs
		if($$VCF_line{'base_type0'} eq 'reference') { $saved_data = "$reference|$reference"; }
		if($$VCF_line{'base_type0'} eq 'snp') { $saved_data = "$consensus|$consensus"; }

		# for phased hets
		if($$VCF_line{'base_type0'} eq 'heterozygous') { 
			#print Dumper($VCF_line);
			die "Error: base1 not saved by VCF_line $line\n" if(!defined $$VCF_line{'0base1'});
			die "Error: base2 not saved by VCF_line: $line\n" if(!defined $$VCF_line{'0base2'});
			$saved_data = ($$VCF_line{'0base1'} . '|' . $$VCF_line{'0base2'});
		}

		$phased_hets{$$VCF_line{'phase_group'}}{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}} = $saved_data;
		$count++;
	}
	close $fh;

	# Add info to summary
	my $summary = "$file: Found " . scalar(keys(%phased_hets)) . " phase groups over $count positions";
	warn "summary\n";
	print $ofh "$summary\n";

	return(\%phased_hets);
}

sub VCF_phased_to_contig_pos_haps_to_variant {
	my ($file, $indels) = @_;
	my %polymorphisms;
	open my $fh, '<', $file or die "Cannot open $file\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);

		# Only save phased variants
		next VCF1 if($$VCF_line{'next'} eq 1);
		next VCF1 if(!defined $$VCF_line{'phased'}); 
		next VCF1 if($$VCF_line{'base_type0'} eq 'ambiguous'); 
		if($indels eq 'n') { next VCF1 if($$VCF_line{'base_type0'} =~ m/het_deletion|het_insertion|insertion|deletion/); }

		# Catalog and save variants
		my ($polymorphism1, $polymorphism2) = &save_two_bases($VCF_line);
		$polymorphisms{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}}{'hap1'} = $polymorphism1;
		$polymorphisms{$$VCF_line{'supercontig'}}{$$VCF_line{'position'}}{'hap2'} = $polymorphism2;
	}
	close $fh;
	return \%polymorphisms;
}

1;

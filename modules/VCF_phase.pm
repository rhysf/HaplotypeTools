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
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;
use read_BAM;

### rfarrer@broadinstitute.org

sub prepare_HT_data {
	my ($platform, $bams, $cutoff, $ref_fasta, $run_grid, $min_depth, $out_folder, $out_VCF, $queue, $max_phase_length, $steps, $sample_names, $VCF_input, $summary) = @_;
	my %HT_data;
	die "Error: prepare_all_io: -a is not UGER, LSF or GridEngine: $platform\n" if($platform !~ m/^(UGER|LSF|GridEngine)$/);
	$HT_data{'platform'} = $platform;
	die "Error: prepare_all_io: -c not between 1-100: $cutoff" if(($cutoff < 1) && ($cutoff > 100));
	$HT_data{'cutoff'} = $cutoff;
	die "Error: prepare_all_io: -f not valid: $ref_fasta\n" unless(-e $ref_fasta);
	$HT_data{'ref_fasta'} = $ref_fasta;
	die "Error: prepare_all_io: -g is not n, N, y or Y: $run_grid\n" if($run_grid !~ m/^(n|N|y|Y)$/);
	$HT_data{'run_grid'} = $run_grid;
	$HT_data{'min_depth'} = $min_depth;
	$HT_data{'out_folder'} = $out_folder;
	$HT_data{'out_VCF'} = $out_VCF;
	$HT_data{'queue'} = $queue;
	$HT_data{'max_phase_length'} = $max_phase_length;
	$HT_data{'steps'} = $steps;
	die "Error: prepare_all_io: -v not valid: $VCF_input\n" unless(-e $VCF_input);
	$HT_data{'VCF_file'} = $VCF_input;
	$HT_data{'VCF_filename'} = fileparse($VCF_input);
	$HT_data{'out_summary'} = $summary;

	# BAMs and sample names
	my @bam_sep = split /,/, $bams;
	my @sample_name_sep = split /,/, $sample_names;
	my %samples_found;
	die "Error: prepare_all_io: different numbers of bams ($bams) and sample names ($sample_names)\n" if(scalar(@bam_sep) ne scalar(@sample_name_sep));
	for(my $i=0; $i<scalar(@bam_sep); $i++) {
		my $bam = $bam_sep[$i];
		my $sample_name = $sample_name_sep[$i];
		die "Error: prepare_all_io: Bam $bam not valid" unless(-e $bam);
		die "Error: prepare_all_io: Sample name $sample_name repeated. Check -u\n" if(defined $samples_found{$sample_name});
		$samples_found{$sample_name} = 1;  

		# Find sample numbers
		my ($sample_number, $sample_found) = vcflines::VCF_and_sample_name_to_sample_number($VCF_input, $sample_name);
		die "Error: prepare_all_io: Sample $sample_name not found in $VCF_input\n" if($sample_found eq 0);

		# Save files and sample numbers
		$HT_data{'BAM_file'}{$sample_number} = $bam_sep[$i];
		$HT_data{'BAM_filename'}{$sample_number} = fileparse($bam_sep[$i]);
		$HT_data{'sample_name'}{$sample_number} = $sample_name;
	}
	return \%HT_data;
}

sub phase_reads {
	my ($alignment_file, $VCF_file, $sample_number) = @_;

	# Save VCF (for 1 contig) to memory (position -> line)
	warn "phase_reads: Finding phase for sites in $VCF_file using $alignment_file...\n";
	my $VCF_hash = vcflines::VCF_to_position_to_line_hash($VCF_file); 

	# Go through every read in the alignment file
	my %phase_positions;
	my $read_count = 0;
	open my $fh, '<', $alignment_file or die "Cannot open $alignment_file : $!\n";
	SAMSEQS: while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($start, $end, $query_dna, $query_start, $query_end, $cigar) = @bits;

		# Remove softclipping from query dna
		my $cigar_soft_clipping_start_count = bamlines::count_softclipping_start($cigar);
		if($cigar_soft_clipping_start_count > 0) {
			#warn "Remove softclipping..\n";
			#warn "line = $line\n";
			$query_dna = substr $query_dna, ($query_start - 1), $query_end;
			#warn "new: $query_dna\n";
		}
		my $read_length = length($query_dna);

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
		#my $phase_position;
		my $length_of_possible_indels_in_read = 0;
		my @VCF_positions = split /\n/, $saved_vcf_positions;
		VCFPOSITIONS: foreach my $vcfline(@VCF_positions) {
			my ($VCF_line) = vcflines::read_VCF_lines($vcfline);
			my ($base_type_id, $sample_id, $genotype_id) = ("base_type$sample_number", "sample_info$sample_number", "GT$sample_number");
			my $base1_id = ($sample_number . "base1");
			my $base2_id = ($sample_number . "base2");

			my $supercontig = $$VCF_line{'supercontig'};
			my $base_type = $$VCF_line{$base_type_id};
			my $position = $$VCF_line{'position'};
			my $ref_base = $$VCF_line{'reference_VCF_format'};
			my $consensus = $$VCF_line{'consensus_VCF_format'};
			my $id = $$VCF_line{'id'};
			my $cons_qual = $$VCF_line{'cons_qual'};
			my $filter = $$VCF_line{'filter'};
			my $info = $$VCF_line{'info'};
			my $format = $$VCF_line{'format'};
			my $sample = $$VCF_line{$sample_id};
			my $genotype = $$VCF_line{$genotype_id};
			my $number_of_samples = $$VCF_line{'number_of_samples'};
			my $polymorphic_position = ($position - $start);
			my $all_bases = "$ref_base,$consensus";
			my $base1 = $$VCF_line{$base1_id};
			my $base2 = $$VCF_line{$base2_id};

			# Only phase heterozygous here
			#next VCFPOSITIONS if($base_type ne 'heterozygous');
			next VCFPOSITIONS if($base_type !~ m/het/);

			# Undefined
			next VCFPOSITIONS if(($base1 eq 'N') || ($base2 eq 'N'));
			#warn "overlapping VCF (het) lines $line\n";

			# If indel, does it go over the end of read
			my $max_length_ref_or_consensus = length($base1);
			if(length($base1) < length($base2)) { $max_length_ref_or_consensus = length($base2); }
			$length_of_possible_indels_in_read += $max_length_ref_or_consensus;
			next VCFPOSITIONS if(($polymorphic_position + $length_of_possible_indels_in_read) > $read_length);

			# No polymorphic bases (snps/hets arent showing up in cigar)\
			#my $no_polymorphic_bases = ($read_length . 'M');
			#if($cigar eq $no_polymorphic_bases) { $base = $ref_base; }
			# But X and = don't seem to present in Samtools. Maybe tool dependent. So can't rely on this. 
			# But will rely on S, M, I and D where possible
			my ($base);

			# Pull out the polymorphic base of heterozygous positions
			if($base_type eq 'heterozygous') {
				$base = substr $query_dna, $polymorphic_position, 1;
			}

			# Pull out the polymorphic base of het_insertion and het_deletion
			elsif($base_type =~ m/insertion|deletion/) {

				# Indel listed in cigar?
				my $indel_in_cigar = bamlines::check_for_indel_in_cigar($cigar, $polymorphic_position);
				my $possible_ref_base = substr $query_dna, $polymorphic_position, length($base1);
				my $possible_consensus_base = substr $query_dna, $polymorphic_position, length($base2);

				# start = correct start
				if(length($base1) eq 1) {

					#warn "possible consensus base = substr $query_dna $polymorphic_position $length_consensus\n";
					#warn "$possible_ref_base and $possible_consensus_base matching $ref_base or $consensus?\n";
					# No match (error or not found in VCF)
					if(($possible_ref_base ne $base1) && ($possible_consensus_base ne $base2)) {
						if($indel_in_cigar eq 'n') { $base = $possible_ref_base; }
						else { $base = $possible_consensus_base; }
					}

					# Something matched
					else {
						# Might be consensus base, or something different
						if($indel_in_cigar eq 'y') { $base = $possible_consensus_base; }
						else {
							if($possible_ref_base eq $base1) { $base = $possible_ref_base; }
							elsif($possible_consensus_base eq $base2) { $base = $possible_consensus_base; }
							else {
								#warn "WARNING 1: Unrecognised base: base type $base_type genotype $genotype 1) $possible_ref_base and $base1 ($ref_base) 2) $possible_consensus_base and $base2 ($consensus) 3) $indel_in_cigar 4) $vcfline 5) $line\n";
								next VCFPOSITIONS;
							}
						}
					}
				}

				# ref base length > 1
				else {

					# No match (error or not found in VCF)
					if(($possible_ref_base ne $base1) && ($possible_consensus_base ne $base2)) {
						if(length($base2) eq 1) {
							# e.g. ref AGG consensus A. Base=A
							if($indel_in_cigar eq 'n') { $base = $possible_consensus_base; }
							else { $base = $possible_ref_base; }
						}
						else {
							#warn "WARNING 2: Unrecognised base: base type $base_type genotype $genotype 1) $possible_ref_base and $base1 ($ref_base) 2) $possible_consensus_base and $base2 ($consensus) 3) $indel_in_cigar 4) $vcfline 5) $line\n";
							next VCFPOSITIONS;
						}
					}

					# Something matched
					else {
						if($indel_in_cigar eq 'y') { $base = $possible_ref_base; }
						else {
							if($possible_ref_base eq $base1) { $base = $possible_consensus_base; }
							elsif($possible_consensus_base eq $base2) { $base = $possible_ref_base; }
							else {
								#warn "WARNING 3: Unrecognised base: base type $base_type genotype $genotype 1) $possible_ref_base and $base1 ($ref_base) 2) $possible_consensus_base and $base2 ($consensus) 3) $indel_in_cigar 4) $vcfline 5) $line\n";
								next VCFPOSITIONS;
							}
						}
					}
				}

				#warn "polymorphic position = $position - $start = $polymorphic_position\n";
				#warn "base = $base\n";
				#die "end here";
			}

			else {
				warn "WARNING: Unrecognised Base type = $base_type\n";
				next VCFPOSITIONS;
				#die "end here";
			}

			# Save phase info in the id column
			$id = &define_phase_position($base, $all_bases, $id, $read_count);

			# Line for printing (leave id blank for read counts & line missing samples columns)
			my $new_line = join "\t", $supercontig, $position, $id, $ref_base, $consensus, $cons_qual, $filter, $info, $format;

			# Samples
			for(my $i=0; $i < $number_of_samples; $i++) {
				my $sample_info_name = "sample_info$i";
				my $sample_info = $$VCF_line{$sample_info_name};
				$new_line .= "\t$sample_info";
			}

			# Replace the entry
			$$VCF_hash{$position} = $new_line;

			# Check it has worked
			my @line_parts = split /\t/, $vcfline;
			my @new_line_parts = split /\t/, $new_line;
			die "Error. New line has different number of parts to old line:\n$vcfline\n$new_line\n" if(scalar(@line_parts) ne scalar(@new_line_parts));
		}
		$read_count++;
	}
	#print Dumper($VCF_hash);
	return ($VCF_hash);
}

sub define_phase_position {
	my ($base, $all_bases, $id, $read_count) = @_;

	# Add new read count
	$id .= "$read_count";

	# Define phase position (0=ref, 1=1st consensus etc.)
	my $phase_found = 0;
	my @bases = split /,/, $all_bases;
	BASES: for(my $i=0; $i<scalar(@bases); $i++) {
		my $VCF_base = $bases[$i];
		if($base eq $VCF_base) { 
			$id .= ('-PP-' . $i . ';'); 
			$phase_found = 1;
			last BASES;
		}
	}

	# Unknown phase
	if($phase_found eq 0) { $id .= ('-NT-' . $base . ';'); }
	return $id;
}

sub VCF_phased_id_read_numbers_to_summary {
	my $id = $_[0];

	my %phase_summary;
	my @read_numbers = split /\;/, $id;
	$phase_summary{'read_count'} = scalar(@read_numbers);
	my $other_bases = 0;

	# Go through each read number and save info
	READNUMS: for(my $i=0; $i<$phase_summary{'read_count'}; $i++) {
		my $read_nums = $read_numbers[$i];

		# Phase
		if($read_nums =~ m/(\d+)-(PP)-(\d+)/) {
			$phase_summary{'phase_group'}{$3}++;
			$phase_summary{'reads_in_phase_group'}{$3}{$1} = 1;
		}

		# Different nucleotide
		elsif($read_nums =~ m/(\d+)-(NT-)([ACTGNactgn])/) {
			$phase_summary{'other_bases'}{$3}++;
			$phase_summary{'reads_in_other_bases'}{$3}{$1} = 1;
		}
		# Now have things like sample-1-phased
		else { }
		#else { die "VCF_phased_id_read_numbers_to_summary: Don't recognise $read_nums\n"; }
	}

	# Move smaller phase groups (<max_haplotypes) to other_bases and reads_in_other_bases
	my $max_haplotypes = 2;
	if(scalar(keys(%{$phase_summary{'phase_group'}})) > $max_haplotypes) {

		# Save size -> phase group
		my %read_count_to_phase_group;
		foreach my $phase_group(keys %{$phase_summary{'phase_group'}}) {
			my $read_count = $phase_summary{'phase_group'}{$phase_group};

			# if tally already defined, add a small amount to the number to make unique
			while(defined $read_count_to_phase_group{$read_count}) {
				$read_count+=0.01;
				last if(!defined $read_count_to_phase_group{$read_count});
			}
			#warn "$phase_group = $read_count\n";
			$read_count_to_phase_group{$read_count} = $phase_group;
		}

		# Move small phase groups to other_bases and reads_in_other_bases
		my $count = 0;
		foreach my $read_count(sort { $b <=> $a } keys %read_count_to_phase_group) {
			my $phase_group = $read_count_to_phase_group{$read_count};
			$count++;
			if($count > $max_haplotypes) {
				#warn "remove $read_count ($phase_group)\n";
				$phase_summary{'other_bases'}{$phase_group} = $read_count;
				delete $phase_summary{'phase_group'}{$phase_group};

				# Reads
				foreach my $read(keys %{$phase_summary{'reads_in_phase_group'}{$phase_group}}) {
					$phase_summary{'reads_in_other_bases'}{$phase_group}{$read} = 1;
				}
				delete $phase_summary{'reads_in_phase_group'}{$phase_group};
			}
		}
	}

	# Save read count to phase group
	foreach my $phase(keys %{$phase_summary{'phase_group'}}) {
		my $tally = $phase_summary{'phase_group'}{$phase};
		# if tally already defined, add a small amount to the number to make unique
		while(defined $phase_summary{'read_count_to_phase_group'}{$tally}) {
			$tally+=0.01;
			last if(!defined $phase_summary{'read_count_to_phase_group'}{$tally});
		}
		$phase_summary{'read_count_to_phase_group'}{$tally} = $phase;
	}
	foreach my $other_bases(keys %{$phase_summary{'other_bases'}}) {
		my $tally = $phase_summary{'other_bases'}{$other_bases}; 
		# if tally already defined, add a small amount to the number to make unique
		while(defined $phase_summary{'read_count_to_phase_group'}{$tally}) {
			$tally+=0.01;
			last if(!defined $phase_summary{'read_count_to_phase_group'}{$tally});
		}
		$phase_summary{'read_count_to_phase_group'}{$tally} = $other_bases;
	}

	# Save number of bases in top two haplotypes
	my $count = 0;
	foreach my $read_count(sort { $b <=> $a } keys %{$phase_summary{'read_count_to_phase_group'}}) {
		# only check for 2 haplotypes
		last if($count eq 2);
		$count++;

		# Exclude any 0.01 extra values for uniqueness
		$read_count = int($read_count);

		# save
		$phase_summary{'top_two_haplotype_depth'} += $read_count;
	}	
	return (\%phase_summary);
}

sub VCF_phase_summaries_to_match {
	my ($summary1, $summary2) = @_;

	# Find overlap
	my %overlap;
	foreach my $PG1(keys %{$$summary1{'reads_in_phase_group'}}) {
		foreach my $PG2(keys %{$$summary2{'reads_in_phase_group'}}) {
			foreach my $read_num(keys %{$$summary1{'reads_in_phase_group'}{$PG1}}) {
				next if(!defined $$summary2{'reads_in_phase_group'}{$PG2}{$read_num});
				$overlap{$PG1}{$PG2}++;
			}
		}
	}

	# Find highest read count match (what do phase group1 bases join to in phase group2)
	my %phase_match;
	my %tmp_matches;
	foreach my $PG1(keys %overlap) {
		my $PG2_count = 0;
		my $PG2_match;
		foreach my $PG2(keys %{$overlap{$PG1}}) {
			my $tally = $overlap{$PG1}{$PG2};
			if($tally > $PG2_count) {
				$PG2_count = $tally;
				$PG2_match = $PG2;
			}
		}
		$phase_match{$PG1} = $PG2_match;
		$tmp_matches{$PG2_match} = 1;
	}

	# Find if both alleles are present in the highest read count matches (e.g. 0->1 and 1->0 or 0->0 and 1->1)
	my $phase_matches = scalar(keys(%tmp_matches));
	return (\%phase_match, $phase_matches);
}

sub save_two_bases {
	my ($VCF_line, $sample_number) = @_;
	my ($polymorphism1, $polymorphism2);
	my $base_type_id = ('base_type' . $sample_number);
	die "save_two_bases: Unable to find sample number $sample_number from VCF line. Check VCF\n" if(!defined $$VCF_line{$base_type_id});
	my $base_type = $$VCF_line{$base_type_id};

	# Save bases
   	if($base_type eq 'reference') {
		$polymorphism1 = $$VCF_line{'reference_VCF_format'};
		$polymorphism2 = $$VCF_line{'reference_VCF_format'};
	}
	elsif($base_type eq 'snp') {
		$polymorphism1 = $$VCF_line{'consensus_VCF_format'};
		$polymorphism2 = $$VCF_line{'consensus_VCF_format'};
	}
	elsif($base_type =~ m/heterozygous|het_deletion|het_insertion|insertion|deletion/) {
		my $base1_id = ($sample_number . 'base1');
		my $base2_id = ($sample_number . 'base2');
		$polymorphism1 = $$VCF_line{$base1_id};
		$polymorphism2 = $$VCF_line{$base2_id};
	} 
	else { 
		print Dumper($VCF_line);
		die "save_two_bases: What is this line. Did not expect $base_type\n"; 
	}
	if((!defined $polymorphism1) || (!defined $polymorphism2)) { die "Something strange over $$VCF_line{'base_type0'}\n"; }
	return ($polymorphism1, $polymorphism2);
}

sub VCF_phased_to_phase_group_contig_pos_to_bases {
	my ($file, $ofh, $sample_number) = @_;
	my %phased_hets;
	my $count = 0;
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		my $contig = $$VCF_line{'supercontig'};
		my $position = $$VCF_line{'position'};
		my $reference = $$VCF_line{'reference_VCF_format'};
		my $consensus = $$VCF_line{'consensus_VCF_format'};
		my $saved_data = '-';
		my $phased_id = ($sample_number . 'phased');
		my $base_type_id = ('base_type' . $sample_number);
		my $base_type = $$VCF_line{$base_type_id};
		my $base1_id = ($sample_number . 'base1');
		my $base2_id = ($sample_number . 'base2');

		# Check it's phased
		next VCF1 if($$VCF_line{'next'} eq 1);
		next VCF1 if(!defined $$VCF_line{$phased_id});

		# for phased homs
		if($base_type eq 'reference') { $saved_data = "$reference|$reference"; }
		if($base_type eq 'snp') { $saved_data = "$consensus|$consensus"; }

		# for phased hets
		if($base_type eq 'heterozygous') { 
			#print Dumper($VCF_line);
			die "Error: base1 not saved by VCF_line $line\n" if(!defined $$VCF_line{$base1_id});
			die "Error: base2 not saved by VCF_line: $line\n" if(!defined $$VCF_line{$base2_id});
			$saved_data = ($$VCF_line{$base1_id} . '|' . $$VCF_line{$base2_id});
		}

		$phased_hets{$$VCF_line{'phase_group'}}{$contig}{$position} = $saved_data;
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
	my ($file, $indels, $sample_name) = @_;
	my %polymorphisms;

	# Find sample number
	my ($sample_number, $sample_found) = vcflines::VCF_and_sample_name_to_sample_number($file, $sample_name);

	open my $fh, '<', $file or die "Cannot open $file\n";
	warn "Reading $file...\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = vcflines::read_VCF_lines($line);
		my $phased_id = ($sample_number . 'phased');
		my $base_type_id = ('base_type' . $sample_number);
		my $supercontig = $$VCF_line{'supercontig'};
		my $position = $$VCF_line{'position'};

		# Only save phased variants
		next VCF1 if($$VCF_line{'next'} eq 1);
		next VCF1 if(!defined $$VCF_line{$phased_id}); 
		next VCF1 if($$VCF_line{$base_type_id} eq 'ambiguous'); 
		if($indels eq 'n') { next VCF1 if($$VCF_line{$base_type_id} =~ m/het_deletion|het_insertion|insertion|deletion/); }

		# Catalog and save variants
		my ($polymorphism1, $polymorphism2) = &save_two_bases($VCF_line, $sample_number);
		$polymorphisms{$supercontig}{$position}{'hap1'} = $polymorphism1;
		$polymorphisms{$supercontig}{$position}{'hap2'} = $polymorphism2;
	}
	close $fh;
	return \%polymorphisms;
}


1;

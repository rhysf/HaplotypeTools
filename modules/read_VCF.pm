package vcflines;
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

### rfarrer@broadinstitute.org

sub read_VCF_lines {
	my $VCF_line = $_[0];
	my @bits = split /\t/, $VCF_line;
	my %VCF_info;

	# Save header info including isolate names
	if($VCF_line =~ m/^\#/) { 
		my $VCF_struct = &VCF_header_to_struct($VCF_line, \%VCF_info); 
		return $VCF_struct; 
	} else {
		$VCF_info{'next'}=0;
		$VCF_info{'header'}='N';
	}

	# Initial quality check
	if(@bits < 9) {
		warn "$0: Bad VCF with < 9 columns: $VCF_line\n";
		$VCF_info{'next'}=1; 
		return \%VCF_info;
	}

	# Multi VCF
	if(@bits > 10) {
		for(my $i=9; $i < scalar(@bits); $i++) {
			my $sample_info_id = ('sample_info' . ($i - 9));
			#warn "$sample_info_id = $i = $bits[$i]\n";
			$VCF_info{$sample_info_id} = $bits[$i];
		}
	}

	# Parts continued
	$VCF_info{'supercontig'}          = $bits[0];
	$VCF_info{'position'}             = $bits[1];
	$VCF_info{'id'}                   = $bits[2];
	$VCF_info{'reference_VCF_format'} = $bits[3];
	$VCF_info{'consensus_VCF_format'} = $bits[4];
	$VCF_info{'cons_qual'}            = $bits[5];
	$VCF_info{'filter'}               = $bits[6];
	$VCF_info{'info'}                 = $bits[7];
	$VCF_info{'format'}               = $bits[8];

	# Split format parts, and save available names 
	my @format_parts = split /:/, $VCF_info{'format'};
	my %format_part_ids_available;
	foreach(@format_parts) { $format_part_ids_available{$_} = 1; }
	$VCF_info{'number_of_samples'} = (scalar(@bits) - 9);

	# Sample_info deals with Multi VCF as well
	for(my $i=9; $i < scalar(@bits); $i++) {

		# Save sample_info[isolate number]
		my $isolate_number = ($i - 9);
		my $sample_info_id = ('sample_info' . $isolate_number);
		$VCF_info{$sample_info_id} = $bits[$i];
		my @sample_info_parts = split /:/, $VCF_info{$sample_info_id};

		# Check format parts matches sample info parts. Error if less format parts than sample info parts
		die "$VCF_line: scalar(@sample_info_parts) < scalar(@format_parts) (sample info < format parts). Have not coded for this eventuality!\n" if(scalar(@format_parts) < scalar(@sample_info_parts));
		# Subset sample_info_parts by format_parts
		if(scalar(@format_parts) > scalar(@sample_info_parts)) {
			my @reduced_sample_info_parts;
			for(my $i=0; $i<scalar(@format_parts); $i++) {
				push @reduced_sample_info_parts, $sample_info_parts[$i];
			}
			@sample_info_parts = @reduced_sample_info_parts;
		}
		die "Format and Sample Info do not match: $VCF_line\n" if(scalar(@format_parts) ne scalar(@sample_info_parts));

		# Save genotype, depth etc.
		my ($GT_id, $DP_id, $base_type_ID, $amb_char_ID, $pid) = ("GT$isolate_number", "DP$isolate_number", "base_type$isolate_number", "amb_char$isolate_number", "PID$isolate_number");
		for(my $f=0; $f<scalar(@format_parts); $f++) {
			my $format_part = $format_parts[$f];
			my $sample_part = $sample_info_parts[$f];
			my $format_part_id = ($format_part . $isolate_number);
			die "$format_part_id overwrites known format id on $VCF_line\n" if(defined $format_part_ids_available{$format_part_id});
			$VCF_info{$format_part_id} = $sample_part;
		}
		die "Unable to find genotype ($GT_id): $VCF_line\n" if(!defined $VCF_info{$GT_id});
		if(!defined $VCF_info{$DP_id}) { $VCF_info{$DP_id} = '?'; }

		# Check alleles are different
		if($VCF_info{$GT_id} =~ m/\/|\|/) {
			my @allele_parts = split /\/|\|/, $VCF_info{$GT_id};
			die "Have not coded for multiple alleles for $VCF_line in check alleles are different\n" if(scalar(@allele_parts) ne 2);
			if($allele_parts[0] eq $allele_parts[1]) { 
				$VCF_info{$GT_id} = $allele_parts[0]; 
			}
		}

		# Determine base (base2 only if diploid)
		my ($base1, $base2, $base_type) = &VCF_struct_determine_bases_and_base_type(\%VCF_info, $GT_id);
		$VCF_info{$base_type_ID}= $base_type;
		$VCF_info{($isolate_number . 'base1')} = $base1;
		$VCF_info{($isolate_number . 'base2')} = $base2;
		
		# Ambiguity character
		if($VCF_info{$base_type_ID} eq 'heterozygous') { $VCF_info{$amb_char_ID} = &get_ambiguity_char($base1, $base2); }

		# Phasing
		if(defined $VCF_info{$pid}) {
			$VCF_info{'phased'} = 1;
			$VCF_info{'phase_group'} = $VCF_info{$pid}; # Not necessary. It's already saved. but for now, keep it so it works with my phasing
			#die "new code. I've found PID: $VCF_info{$pid} . All good?\n";
		}
	}	

	# Return
	return \%VCF_info;
}

# Phase metrics from my phasing
sub get_phase_metrics {
	my ($phase_info, $current_phase_group, $current_position, $current_supercontig) = @_;

	# Entend the previous phase group
	if(($$phase_info{'previous_group'}) && ($$phase_info{'previous_group'} eq $current_phase_group)) {
		$$phase_info{'previous_last_position'} = $current_position;

		# Remove save
		if(defined $$phase_info{'save'}) { delete $$phase_info{'save'}; }
	}

	# New phase group
	else {
		# Save
		if(defined $$phase_info{'previous_group'}) {
			$$phase_info{'save'}{'group'} = $$phase_info{'previous_group'};
			$$phase_info{'save'}{'supercontig'} = $$phase_info{'previous_supercontig'};
			$$phase_info{'save'}{'first_position'} = $$phase_info{'previous_first_position'};
			$$phase_info{'save'}{'last_position'} = $$phase_info{'previous_last_position'};
		}

		# New
		$$phase_info{'previous_group'} = $current_phase_group;
		$$phase_info{'previous_supercontig'} = $current_supercontig;
		$$phase_info{'previous_first_position'} = $current_position;
		$$phase_info{'previous_last_position'} = $current_position;
	}
	return $phase_info;
}

sub summarise_phases {
	my $phase_sizes = $_[0];
	my ($min, $max, $mean, $total, $num) = (0,0,0,0,0);
	foreach my $length(sort { $a <=> $b } keys %{$phase_sizes}) {
		my $tally = $$phase_sizes{$length};
		if(!defined $min) { $min = $length; }
		if($length > $max) { $max = $length; }
		$num += $tally;
		$total += ($length * $tally);
		#warn "in summarise_phases with $length and $tally\n";
	}
	if($num eq 0) { $mean = 0; }
	else { $mean = ($total / $num); }
	return ($min, $max, $mean, $total, $num);
}

sub VCF_split_for_phasing {
	my ($file, $outfolder, $haplotype_length, $contig_lengths) = @_;

	# Save VCF
	warn "VCF_split_for_phasing: $file...\n";
	my (%VCF_header);
	my $ofh;
	my $contig;
	my ($start_window, $stop_window) = (0, $haplotype_length);
	open my $fh, '<', $file or die "Cannot open VCF: $file: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Save the header of VCF. Ignore ambiguous sites
		if($$VCF_line{'header'} eq 'Y') {
   			if ($line =~ m/^\#/) {
				if($line =~ m/=/) {
					my @header_parts = split /=/, $line;
					# ignore previous PID or PGT
					next VCF1 if($line =~ m/ID=PID|ID=PGT/);

					$VCF_header{$header_parts[0]} .= "$line\n";
				} else { 
					$VCF_header{'rest'} .= "$line\n";
				}
			}
			next VCF1;
		}
		next VCF1 if($$VCF_line{'next'} eq 1);
		next VCF1 if($$VCF_line{'base_type0'} =~ m/insertion|deletion/);
		#next VCF1 if($$VCF_line{'GT0'} eq '.');
		die "VCF_split_for_phasing: Multisample VCF detected. Currently only compatible with single sample VCFs" if($$VCF_line{'number_of_samples'} > 1);
		die "VCF_split_for_phasing: Cannot find genotype: $line" if(!defined $$VCF_line{'GT0'});

		# Init outfile
		if(!defined $ofh) {
			$contig = $$VCF_line{'supercontig'};
			warn "\tProcessing $contig...\n";
			my $outfile = "$outfolder/$file-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}
		if($$VCF_line{'supercontig'} ne $contig) {
			close $ofh;
			$contig = $$VCF_line{'supercontig'};
			($start_window, $stop_window) = (0, $haplotype_length);
			warn "\tProcessing $contig...\n";
			my $outfile = "$outfolder/$file-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}

		# Check for new window
		if($$VCF_line{'position'} > $stop_window) {
			close $ofh;
			$start_window += $haplotype_length;
			$stop_window += $haplotype_length; 
			my $outfile = "$outfolder/$file-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}

		# Check for previous phasing and remove
		if($$VCF_line{'format'} =~ m/PID|PGT/) {
			my $new_format;
			my $new_sample_info;
			my @format_parts = split /:/, $$VCF_line{'format'};
			my @sample_parts = split /:/, $$VCF_line{'sample_info0'};
			for(my $i=0; $i<scalar(@format_parts); $i++) {
				my $format_part = $format_parts[$i];
				my $sample_part = $sample_parts[$i];
				if(($format_part ne 'PID') && ($format_part ne 'PGT')) {
					$new_format .= "$format_part:";
					$new_sample_info .= "$sample_part:";
				}
			}
			$new_format =~ s/:$//;
			$new_sample_info =~ s/:$//;
			$$VCF_line{'format'} = $new_format;
			$$VCF_line{'sample_info0'} = $new_sample_info;
		}
		if($$VCF_line{'sample_info0'} =~ m/\|/) { $$VCF_line{'sample_info0'} =~ s/\|/\//g; }

		# Print (leave id blank for read counts
		my $new_line = join "\t", $$VCF_line{'supercontig'}, $$VCF_line{'position'}, '', $$VCF_line{'reference_VCF_format'}, $$VCF_line{'consensus_VCF_format'}, $$VCF_line{'cons_qual'}, $$VCF_line{'filter'}, $$VCF_line{'info'}, $$VCF_line{'format'}, $$VCF_line{'sample_info0'};
		print $ofh "$new_line\n";
	}
	close $fh;
	close $ofh;

	# Make new Phase block info
	$VCF_header{'##FORMAT'} .= "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">\n";

	# Print header
	warn "VCF_split_for_phasing: Printing header...\n";
	my $header_output = "$outfolder/$file-header";
	open my $ofh2, '>', $header_output or die "Cannot open $header_output";
	if(defined $VCF_header{'##fileformat'}) { print $ofh2 "$VCF_header{'##fileformat'}"; }
	foreach my $header_part(keys %VCF_header) {
		next if($header_part eq '##fileformat');
		next if($header_part eq 'rest');
		print $ofh2 "$VCF_header{$header_part}";
	}
	print $ofh2 "$VCF_header{'rest'}";
	close $ofh2;
	return;
}

sub VCF_to_position_to_line_hash {
	my ($file) = @_;
	#warn "Saving VCF file: $file\n";
	my (%VCF);
	open my $fh, '<', $file or die "Cannot open VCF: $file: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);
		$VCF{$$VCF_line{'position'}} = $line;
	}
	close $fh;
	return (\%VCF);
}

sub get_ambiguity_char {
	my ($base1, $base2) = @_;
	my $ambiguity_char;
				
	# K
	if(($base1 eq 'T') && ($base2 eq 'G')) { $ambiguity_char = 'K'; }
	elsif(($base1 eq 'G') && ($base2 eq 'T')) { $ambiguity_char = 'K'; }
							
	# M
	elsif(($base1 eq 'A') && ($base2 eq 'C')) { $ambiguity_char = 'M'; }
	elsif(($base1 eq 'C') && ($base2 eq 'A')) { $ambiguity_char = 'M'; }
										
	# R
	elsif(($base1 eq 'A') && ($base2 eq 'G')) { $ambiguity_char = 'R'; }
	elsif(($base1 eq 'G') && ($base2 eq 'A')) { $ambiguity_char = 'R'; }
													
	# Y
	elsif(($base1 eq 'T') && ($base2 eq 'C')) { $ambiguity_char = 'Y'; }
	elsif(($base1 eq 'C') && ($base2 eq 'T')) { $ambiguity_char = 'Y'; }
																
	# S
	elsif(($base1 eq 'G') && ($base2 eq 'C')) { $ambiguity_char = 'S'; }
	elsif(($base1 eq 'C') && ($base2 eq 'G')) { $ambiguity_char = 'S'; }
																			
	# W
	elsif(($base1 eq 'A') && ($base2 eq 'T')) { $ambiguity_char = 'W'; }
	elsif(($base1 eq 'T') && ($base2 eq 'A')) { $ambiguity_char = 'W'; }
	return $ambiguity_char;
}

############### Local subroutines
sub VCF_header_to_struct {
	my ($VCF_line, $VCF_struct) = @_;
	my @bits = split /\t/, $VCF_line;
	$$VCF_struct{'next'}=1; 
	$$VCF_struct{'header'}='Y';
	if($VCF_line =~ m/^\#CHROM\tPOS\tID\tREF/) {
		for(my $i=9; $i < scalar(@bits); $i++) {
			$$VCF_struct{'isolate_names'}{($i - 9)} = $bits[$i];
		}
	}
	return $VCF_struct;
}

sub VCF_struct_determine_bases_and_base_type {
	my($VCF_struct, $GT_id) = @_;
	my ($base1, $base2, $base_type);
	$base2 = 'None';

	# ambiguous
	if(($$VCF_struct{'reference_VCF_format'} eq 'N') || ($$VCF_struct{'consensus_VCF_format'} eq 'N') || ($$VCF_struct{$GT_id} eq '.')) { 
		$base1 = 'N';
		$base_type = 'ambiguous';
		return ($base1, $base2, $base_type);
	}

	# Homozygous ref-calls
	if($$VCF_struct{$GT_id} eq 0) { 
		$base1 = $$VCF_struct{'reference_VCF_format'};
		$base_type = 'reference';
		return ($base1, $base2, $base_type);
	}

	# QC that GT matches a base for homozygous variants
	my @bases = split /,/, $$VCF_struct{'consensus_VCF_format'};
	if(($$VCF_struct{$GT_id} !~ m/(\d)([\/\|])(\d)/) && (!defined $bases[($$VCF_struct{$GT_id} - 1)])) {
		warn "Nothing found for this VCF entry:\n";
		warn Dumper($VCF_struct);
		$base1 = 'N';
		$base_type = 'ambiguous';
		return ($base1, $base2, $base_type);
	}

	# Not heterozygous
	if(($$VCF_struct{$GT_id} ne 0) && ($$VCF_struct{$GT_id} !~ m/(\d)([\/\|])(\d)/)) {
		my $consensus = $bases[($$VCF_struct{$GT_id} - 1)]; # won't be defined for heterozygous positions
	
		# Homozygous SNP
		if(length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) { 

			# A SNP
			if((length($$VCF_struct{'reference_VCF_format'}) eq 1) && (length($consensus) eq 1)) { 
				$base1 = $consensus;
				$base_type = 'snp';
				return ($base1, $base2, $base_type);
			}
	
			# SNP(s) disguised as an indel
			if((length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) && ($consensus !~ m/\./)) {
				my @bases_reference = split //, $$VCF_struct{'reference_VCF_format'};
				my @bases_consensus = split //, $consensus;
				my $snp_count = 0;
				for(my $i=0; $i<scalar(@bases_reference); $i++) {
					my $ref_base = $bases_reference[$i];
					my $cons_base = $bases_consensus[$i];
					if($ref_base ne $cons_base) { $snp_count++; }
				}
				if($snp_count ne 0) {
					$base1 = $consensus;
					$base_type = ('snp_multi' . $snp_count);
					return ($base1, $base2, $base_type);
				}
			}
	
			# Ambiguous
			warn "Nothing found for this apparant homozygous snp:\n";
			warn Dumper($VCF_struct);
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}

		# Homozygous indel
		if(length($$VCF_struct{'reference_VCF_format'}) ne length($consensus)) {

			# Deletion (maybe with snps in there too!)
			if(length($$VCF_struct{'reference_VCF_format'}) > length($consensus)) { 
				$base1 = $consensus;
				$base_type = 'deletion';
				return ($base1, $base2, $base_type);
			}
			if((length($$VCF_struct{'reference_VCF_format'}) eq length($consensus)) && ($consensus =~ m/^\./)) { 
				$base1 = $consensus;
				$base_type = 'deletion';
				return ($base1, $base2, $base_type);
			}	
	
			# Insertion (maybe with snps in there too!)
			if(length($$VCF_struct{'reference_VCF_format'}) < length($consensus)) { 
				$base1 = $consensus;
				$base_type = 'insertion';
				return ($base1, $base2, $base_type);
			}
			
			# Ambiguous
			warn "Nothing found for this apparent homozygous indel:\n";
			warn Dumper($VCF_struct);
			$base1 = 'N';
			$base_type = 'ambiguous';
			return ($base1, $base2, $base_type);
		}
	}

	# Bi-allelic heterozygous positions & indels
	if($$VCF_struct{$GT_id} =~ m/(\d)([\/\|])(\d)/) { 
		$base_type = 'heterozygous';
		my @bases_het;
		if($$VCF_struct{'consensus_VCF_format'} =~ m/\,/) {
			@bases = split /,/, $$VCF_struct{'consensus_VCF_format'};
			foreach(@bases) {
				if(length($_) > length($$VCF_struct{'reference_VCF_format'})) { $base_type = 'het_insertion'; }
				if(length($_) < length($$VCF_struct{'reference_VCF_format'})) { $base_type = 'het_deletion'; }
			}
		} else { 
			push @bases, $$VCF_struct{'reference_VCF_format'};
			push @bases, $$VCF_struct{'consensus_VCF_format'};
			if(length($bases[1]) > length($bases[0])) { $base_type = 'het_insertion'; }
			if(length($bases[1]) < length($bases[0])) { $base_type = 'het_deletion'; }
		}
		$base1 = $bases[0];
		$base2 = $bases[1];
		return ($base1, $base2, $base_type);
	}
	#return ($base1, $base2, $base_type);
}

1;

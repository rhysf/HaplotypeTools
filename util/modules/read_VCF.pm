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
use File::Basename;

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
		die "$VCF_line: scalar(@format_parts) < scalar(@sample_info_parts) (format parts < sample info parts). Have not coded for this eventuality!\n" if(scalar(@format_parts) < scalar(@sample_info_parts));

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
		my ($GT_id, $DP_id, $base_type_ID, $amb_char_ID, $pid, $ps_id) = ("GT$isolate_number", "DP$isolate_number", "base_type$isolate_number", "amb_char$isolate_number", "PID$isolate_number", "PS$isolate_number");
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
		if(!defined $VCF_info{$base_type_ID}) { die "Error: read_VCF not identified base type ($base_type_ID) from $VCF_line\n"; } 
		if($VCF_info{$base_type_ID} eq 'heterozygous') { $VCF_info{$amb_char_ID} = &get_ambiguity_char($base1, $base2); }

		# Phasing
		if(defined $VCF_info{$pid}) {
			$VCF_info{($isolate_number . 'phased')} = 1;
			$VCF_info{($isolate_number . 'phase_group')} = $VCF_info{$pid}; # Not necessary. It's already saved. but for now, keep it so it works with my phasing
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
	my ($HT_data, $contig_lengths) = @_;

	# Input
	my $file = $$HT_data{'VCF_file'};
	my $outfolder = $$HT_data{'out_folder'};
	my $haplotype_length = $$HT_data{'max_phase_length'};
	foreach($$HT_data{'VCF_file'}, $$HT_data{'out_folder'}, $$HT_data{'max_phase_length'}) { die "VCF_split_for_phasing: HT_data not properly initialised: $!\n" if(!defined $_); }

	# Save VCF
	warn "VCF_split_for_phasing: $file...\n";
	my $file_no_dir = fileparse($file);
	my (%VCF_header, %isolate_names);
	my $ofh;
	my $contig;
	my $sample_number = 0;
	my ($start_window, $stop_window) = (0, $haplotype_length);
	open my $fh, '<', $file or die "Cannot open VCF: $file: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Save the VCF header (except any previous PID or PGT)
		if($$VCF_line{'header'} eq 'Y') {
			if($line =~ m/=/) {
				my @header_parts = split /=/, $line;
				# ignore previous PID or PGT
				next VCF1 if($line =~ m/ID=PID|ID=PGT/);
				$VCF_header{$header_parts[0]} .= "$line\n";
			} else { $VCF_header{'rest'} .= "$line\n"; }
			next VCF1;
		}
		next VCF1 if($$VCF_line{'next'} eq 1);

		# Ignore ambiguous sites
		my ($GT, $base_type) = ("GT$sample_number", "base_type$sample_number");
		die "VCF_split_for_phasing: Cannot find genotype: $line" if(!defined $$VCF_line{$GT});

		# Init outfile
		if(!defined $ofh) {
			$contig = $$VCF_line{'supercontig'};
			warn "\tProcessing $contig...\n";
			my $outfile = "$outfolder/$file_no_dir-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}
		if($$VCF_line{'supercontig'} ne $contig) {
			close $ofh;
			$contig = $$VCF_line{'supercontig'};
			($start_window, $stop_window) = (0, $haplotype_length);
			warn "\tProcessing $contig...\n";
			my $outfile = "$outfolder/$file_no_dir-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}

		# Check for new window
		if($$VCF_line{'position'} > $stop_window) {
			close $ofh;
			$start_window += $haplotype_length;
			$stop_window += $haplotype_length; 
			my $outfile = "$outfolder/$file_no_dir-$contig-$start_window-$stop_window";
			open $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		}

		# Remove PID|PGT from format if present
		my $new_format = &remove_phase_from_format($$VCF_line{'format'});

		# Line for printing (leave id blank for read counts & line missing samples columns)
		my $new_line = join "\t", $$VCF_line{'supercontig'}, $$VCF_line{'position'}, '', $$VCF_line{'reference_VCF_format'}, $$VCF_line{'consensus_VCF_format'}, $$VCF_line{'cons_qual'}, $$VCF_line{'filter'}, $$VCF_line{'info'}, $new_format;

		# Remove PID|PGT from samples and add sham samples if needed (needs to be for all samples or format parts < sample info parts)
		for(my $i=0; $i < $$VCF_line{'number_of_samples'}; $i++) {
			my $sample_info_name = "sample_info$i";
			my $sample_info = $$VCF_line{$sample_info_name};
			$sample_info = &remove_phase_and_optional_add_sham_for_sample($sample_info, $$VCF_line{'format'});
			$new_line .= "\t$sample_info";
		}

		# Print
		print $ofh "$new_line\n";
	}
	close $fh;
	close $ofh;

	# Make new Phase block info
	$VCF_header{'##FORMAT'} .= "##FORMAT=<ID=PID,Number=1,Type=String,Description=\"Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group\">\n";

	# Print header
	warn "VCF_split_for_phasing: Printing header...\n";
	my $header_output = "$outfolder/$file_no_dir-header";
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

sub parse_VCF {
	my ($vcf_file, $printing_option, $feature, $summarise_per_contig, $sample_name) = @_;
	my (%counts, %types_found);
	my $phase_info;

	# Sample number
	my $sample_number = &VCF_and_sample_name_to_sample_number($vcf_file, $sample_name); 

	# Go through the VCF 
	open my $fh, '<', $vcf_file or die "Cannot open $vcf_file\n";
	warn "Reading $vcf_file...\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Print 
		# Header
		if($$VCF_line{'header'} eq 'Y') {
			&print_VCF_header($line, $printing_option);
			next VCF1;
		}
		next VCF1 if($$VCF_line{'next'} eq 1);
		last VCF1 if(($printing_option eq 'header') && ($$VCF_line{'header'} eq 'N'));

		# Print non-header VCF lines
		my ($base_type_id, $amb_char_id) = ("base_type$sample_number", "amb_char$sample_number");
		die "base type for sample number $sample_number not found for $line\n" if(!defined $$VCF_line{$base_type_id});
		my ($base_type, $info, $amb_char) = ($$VCF_line{$base_type_id}, $$VCF_line{'info'}, $$VCF_line{$amb_char_id});
		&print_VCF_lines($line, $printing_option, $feature, $$VCF_line{$base_type}, "$vcf_file-$$VCF_line{'supercontig'}");

		# Tally variants
		my $hash_contig_name = "all";
		if($summarise_per_contig eq 'y') { $hash_contig_name = $$VCF_line{'supercontig'}; }
		$counts{$vcf_file}{$sample_name}{$hash_contig_name}{$base_type}++; 
		$types_found{'Variant_types'}{$base_type} = 1;

		# Special types
		if(($base_type eq 'heterozygous') && ($info =~ m/TRIPROB|TRIALLELIC/)) { $counts{$vcf_file}{$sample_name}{$hash_contig_name}{$info}++; } 
		if(($base_type eq 'snp') && ($info eq 'phased')) { $counts{$vcf_file}{$sample_name}{$hash_contig_name}{'Phased_SNP'}++; }
		if(($base_type eq 'heterozygous') && ($info eq 'phased')) { $counts{$vcf_file}{$sample_name}{$hash_contig_name}{'Phased_HET'}++; }
		if($amb_char) { $types_found{'Type_of_Het'}{$amb_char}{$vcf_file}{$sample_name}{$hash_contig_name}++; }
		next VCF1 if($base_type eq 'ambiguous');
		next VCF1 if($$VCF_line{'consensus_VCF_format'} eq '.');

		# Phase metrics
		if($$VCF_line{'phased'}) {
			($phase_info) = &get_phase_metrics($phase_info, $$VCF_line{'phase_group'}, $$VCF_line{'position'}, $$VCF_line{'supercontig'});

			if(defined $$phase_info{'save'}) {

				# New phase group = push haplotype into array
				my $phase_group_to_save = $$phase_info{'save'}{'group'};
				my $supercontig_of_phase_group_to_save = $$phase_info{'save'}{'supercontig'};
				my $first_phase_group_position = $$phase_info{'save'}{'first_position'};
				my $last_phase_group_position = $$phase_info{'save'}{'last_position'};
				my $length_of_last_phase_group = ($last_phase_group_position - $first_phase_group_position);
				$types_found{'Phase_group_lengths'}{$vcf_file}{$sample_name}{$hash_contig_name}{$length_of_last_phase_group}++;
			}
		}
	}
	close $fh;
	return (\%counts, \%types_found);
}

sub VCF_to_position_to_line_hash {
	my ($file) = @_;
	#warn "Saving VCF file: $file\n";
	my (%VCF);
	open my $fh, '<', $file or die "Cannot open VCF: $file: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my @bits = split /\t/, $line;

		# ignore headers
		next if(!defined $bits[1]); 
		$VCF{$bits[1]} = $line;
	}
	close $fh;
	return (\%VCF);
}

sub VCF_line_make {
	my $VCF_line = $_[0];
	my $number_of_samples = $$VCF_line{'number_of_samples'};
	my $line = join "\t", 
	$$VCF_line{'supercontig'}, 
	$$VCF_line{'position'}, 
	$$VCF_line{'id'}, 
	$$VCF_line{'reference_VCF_format'}, 
	$$VCF_line{'consensus_VCF_format'}, 
	$$VCF_line{'cons_qual'}, 
	$$VCF_line{'filter'}, 
	$$VCF_line{'info'}, 
	$$VCF_line{'format'};

	# Samples
	for(my $i=0; $i < $number_of_samples; $i++) {
		my $sample_info_name = "sample_info$i";
		my $sample_info = $$VCF_line{$sample_info_name};
		$line .= "\t$sample_info";
	}
	return $line;
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

sub VCF_and_sample_name_to_sample_number {
	my ($VCF, $sample_name) = @_;
	die "VCF_and_sample_name_to_sample_number: VCF doesn't look valid: $VCF : $!\n" if(! -e $VCF);
	my $sample_number = 0;
	my $found = 0;
	my $found_header = 0;
	open my $fh, '<', $VCF or die "Cannot open VCF: $VCF: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Save the isolate_names
		if(defined $$VCF_line{'isolate_names'}) {
			$found_header = 1;
			foreach my $vcf_isolate_number(keys %{$$VCF_line{'isolate_names'}}) {
				my $vcf_sample_name = $$VCF_line{'isolate_names'}{$vcf_isolate_number};
				if($vcf_sample_name eq $sample_name) {
					$sample_number = $vcf_isolate_number; 
					$found =1;
				}
			}
			last VCF1;
		}
	}
	close $fh;
	if($found_header eq 0) { warn "Did not find VCF header specifying sample names\n"; }
	if($found eq 0) { warn "Did not find sample $sample_name in VCF header, so using sample 0 by default (if this is a multi VCF and something other than the first sample is wanted - it needs to be re-run\n"; }
	if($found eq 1) { warn "VCF_and_sample_name_to_sample_number: found sample $sample_name on sample column: $sample_number\n"; }
	return ($sample_number, $found);
}

sub VCF_to_sample_names {
	my $VCF = $_[0];
	my %sample_names;
	open my $fh, '<', $VCF or die "Cannot open VCF: $VCF: $!\n";
	VCF1: while (my $line = <$fh>) {
   		chomp $line;
		my ($VCF_line) = &read_VCF_lines($line);

		# Save the isolate_names
		if(defined $$VCF_line{'isolate_names'}) {
			return $$VCF_line{'isolate_names'};
		}
	}
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

sub print_VCF_header {
	my ($line, $printing_option) = @_; 
	print "$line\n" if($printing_option eq 'vcf');
	print "$line\n" if($printing_option eq 'header');
	return 1;
}

sub print_VCF_lines {
	my ($line, $printing_option, $feature, $base_type, $outfile) = @_;
	if(($feature eq 'all') || ($feature eq $base_type)) { 
		print "$line\n" if($printing_option eq 'vcf');
		print "$line\n" if($printing_option eq 'no_header');
		if($printing_option eq 'split_by_contig') { 
			open OUT, ">>$outfile" or die "Cannot open $outfile: $!\n";
			print OUT "$line\n";
			close OUT;
		}
	}
	return 1;
}

sub remove_phase_and_optional_add_sham_for_sample {
	my ($sample, $format) = @_;

	if($sample =~ m/\|/) { $sample =~ s/\|/\//g; }
	my @sample_parts = split /:/, $sample;
	my @format_parts = split /:/, $format;

	# Check for empty samples or ambigious samples and replace with sham (explicitly ambigious!)
	my $new_sample_info = $sample;
	if(scalar(@sample_parts) < scalar(@format_parts)) { $new_sample_info = &make_sham_sample_from_format($format); }

	# Check for previous phasing from old format and remove from samples
	if($format =~ m/PID|PGT/) {
		my @sample_parts = split /:/, $new_sample_info;
		$new_sample_info = '';
		for(my $i=0; $i<scalar(@format_parts); $i++) {
			my $format_part = $format_parts[$i];
			next if($format_part =~ m/PID|PGT/);
			die "No sample part found from $sample $format . Check VCF\n" if(!defined $sample_parts[$i]);
			my $sample_part = $sample_parts[$i];
			$new_sample_info .= "$sample_part:";
		}
		$new_sample_info =~ s/:$//;
	}
	return $new_sample_info;
}

sub make_sham_sample_from_format {
	my $format = $_[0];
	my @format_parts = split /:/, $format;
	my $sham_sample = ".:" x scalar(@format_parts);
	$sham_sample =~ s/:$//;
	return $sham_sample;
}

sub remove_phase_from_format {
	my $format = $_[0];
	my $unphased_format;

	if($format =~ m/PID|PGT/) {
		my @format_parts = split /:/, $format;
		for(my $i=0; $i<scalar(@format_parts); $i++) {
			my $format_part = $format_parts[$i];
			next if($format_part =~ m/PID|PGT/);
			$unphased_format .= "$format_part:";
		}
		$unphased_format =~ s/:$//;
	} else { $unphased_format = $format; }
	return $unphased_format;
}

sub VCF_struct_determine_bases_and_base_type {
	my($VCF_struct, $GT_id) = @_;
	my ($base1, $base2, $base_type);
	my $ref_base = $$VCF_struct{'reference_VCF_format'};
	my $consensus = $$VCF_struct{'consensus_VCF_format'};
	my $GT = $$VCF_struct{$GT_id};

	# Pull all bases
	my @all_bases = ($ref_base, split(/,/, $consensus));
	my @gt_indices = split(/[\/|]/, $GT);
	my @bases = map { $all_bases[$_] // 'N' } @gt_indices;

	# Determine base type

	my %seen; @seen{@gt_indices} = ();
    my @unique = keys %seen;

    #my $base_type;
    #if ($GT eq '.' or grep { $_ eq 'N' } @bases) {
    #    $base_type = 'ambiguous';
    #}

	# ambiguous
	if(($base1 eq 'N') || ($GT eq '.')) { 
		$base_type = 'ambiguous';
		return ($base1, $base2, $base_type);
	}

	# Homozygous reference
	if($GT eq 0) { 
		$base_type = 'reference'; 
		return ($base1, $base2, $base_type);
	}

	# Homozygous variant
	if(($GT ne 0) && ($GT =~ m/^(\d)$/)) {

		# Homozygous SNP
		if((length($ref_base) eq length($base1)) && (length($ref_base) eq 1)) { 
			$base_type = 'snp'; 
			return ($base1, $base2, $base_type);
		}

		# SNP(s) disguised as an indel
		elsif(((length($ref_base) eq length($base1)) && (length($ref_base) ne 1)) && ($base1 !~ m/\./)) {
			my @bases_reference = split //, $ref_base;
			my @bases_consensus = split //, $base1;
			my $snp_count = 0;
			my ($ref_base_saved, $cons_base_saved);
			for(my $i=0; $i<scalar(@bases_reference); $i++) {
				my $ref_base = $bases_reference[$i];
				my $cons_base = $bases_consensus[$i];
				if($ref_base ne $cons_base) { 
					$ref_base_saved = $ref_base;
					$cons_base_saved = $cons_base;
					$snp_count++; 
				}
			}

			# no changes found
			if($snp_count eq 0) { $base_type = 'reference'; }

			# update bases with just those that are variant
			if($snp_count eq 1) {
				$base1 = $ref_base_saved;
				$base2 = $cons_base_saved;
				$base_type = 'snp';
				return ($base1, $base2, $base_type);
			}

			# multiple snps found
			else { $base_type = ('snp_multi' . $snp_count); }
			return ($base1, $base2, $base_type);
		}

		# Homozygous indel
		elsif(length($ref_base) ne length($base1)) {

			# Deletion (maybe with snps in there too!)
			if(length($ref_base) > length($base1)) { $base_type = 'deletion'; }
			if((length($ref_base) eq length($base1)) && ($base1 =~ m/^\./)) { $base_type = 'deletion'; }

			# Insertion (maybe with snps in there too!)
			if(length($ref_base) < length($base1)) { $base_type = 'insertion'; }

			return ($base1, $base2, $base_type);
		}
	}

	# Bi-allelic heterozygous positions & indels
	if($GT =~ m/^(\d)([\/\|])(\d)$/) { 
		my ($GT1, $GT2) = ($1, $3);

		$base_type = 'heterozygous';
		if(length($base1) > length($base2)) { $base_type = 'het_insertion'; }
		if(length($base1) < length($base2)) { $base_type = 'het_deletion'; }
		return ($base1, $base2, $base_type);
	}

	return ($bases[0] // 'N', $bases[1] // 'N', $base_type);
}

1;

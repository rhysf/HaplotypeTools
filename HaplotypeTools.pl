#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use lib "/home/unix/rfarrer/perl5/lib/perl5/x86_64-linux-thread-multi/";
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_BAM;
use read_FASTA;
use read_VCF;
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF> -b <sorted BAMs (separated by comma)> -u <VCF sample names in order of input BAM files (separated by comma)> -f <reference fasta>
Optional: -c\tCut-off percent reads supporting phase group [90]
          -m\tMinimum read depth overlapping two heterozygous positions [4]
	  -r\tMax phase length [10000]
	  -s\tSteps (1=process VCF, 2=process BAM, 3=assign read info to VCF, 4=validate and assign phase groups, 5=concatenate) [12345]\n
Parallel: -g\tRun commands on the grid (y/n) [n]
	  -a\tPlatform (UGER, LSF, GridEngine) [UGER]
	  -q\tQueue name [short]\n
Outputs:  -o\tOutput folder for tmp files [opt_v-HaplotypeTools-phased-r-opt_r]
          -p\tPhased VCF [opt_v-Phased-m-opt_m-c-opt_c-r-opt_r.vcf]
	  -y\tPhased summary [opt_v-Phased-m-opt_m-c-opt_c-r-opt_r-summary.tab]\n";
our($opt_a, $opt_b, $opt_c, $opt_f, $opt_g, $opt_m, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_u, $opt_v, $opt_y);
getopt('abcfgmopqrsuvy');
die $usage unless ($opt_v && $opt_b && $opt_u && $opt_f);
if(!defined $opt_a) { $opt_a = 'UGER'; }
if(!defined $opt_c) { $opt_c = 90; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_r) { $opt_r = 10000; }
if(!defined $opt_s) { $opt_s = 12345; }
if(!defined $opt_q) { $opt_q = 'short'; }
if(!defined $opt_o) { $opt_o = $opt_v . '-HaplotypeTools-phased-r-' . $opt_r; }
if(!defined $opt_p) { $opt_p = "$opt_v-Phased-m-$opt_m-c-$opt_c-r-$opt_r.vcf"; }
if(!defined $opt_y) { $opt_y = "$opt_v-Phased-m-$opt_m-c-$opt_c-r-$opt_r-summary.tab"; }

# Prepare all input and output variables
my $HT_data = vcfphase::prepare_HT_data($opt_a, $opt_b, $opt_c, $opt_f, $opt_g, $opt_m, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_u, $opt_v, $opt_y); 
print Dumper($HT_data);

# Output folder
if(! -d $$HT_data{'out_folder'}) { my $cmd1 = `mkdir $$HT_data{'out_folder'}`; }

# Dependencies
my $VCF_phase_script = "$Bin/util/VCF_phase.pl";
my $VCF_phase_validate_and_assign_script = "$Bin/util/VCF_phased_validate_and_assign_phase_groups.pl";
my $Run_Commands_python = "$Bin/util/Run_cmds_on_grid.py";
foreach($VCF_phase_script, $VCF_phase_validate_and_assign_script, $Run_Commands_python) { die "Cannot find dependency script $_ : $!\n" if(! -e $_); }

# Save names and length of reference sequence
my $fasta_lengths = fastafile::fasta_id_to_seq_length_hash($$HT_data{'ref_fasta'});

# Step 1: Split VCF
if($$HT_data{'steps'} =~ m/1/) {
	warn "Step 1: Split VCF...\n";
	my $VCF_split = vcflines::VCF_split_for_phasing($HT_data, $fasta_lengths);
}

# Step 2: Split BAMs
if($$HT_data{'steps'} =~ m/2/) {
	warn "Step 2: Split BAMs...\n";
	my $BAM_split = bamlines::BAM_split_for_phasing($HT_data, $fasta_lengths);
}

# VCF header
my @header_files = <$$HT_data{'out_folder'}/*-header>;
my $header_file = $header_files[0];

# Step 3: Assign read info to VCFs
if($$HT_data{'steps'} =~ m/3/) {
	warn "Step 3: Assign read info to VCFs...\n";

	# Foreach BAM
	my @cmds;
	foreach my $sample_number(sort keys %{$$HT_data{'BAM_filename'}}) {
		my $BAM_filename = $$HT_data{'BAM_filename'}{$sample_number};
		#warn "\t$BAM_filename and $VCF_filename ($sample_number) ->\n";

		# Foreach contig
		CONTIGS: foreach my $contig(sort keys %{$fasta_lengths}) {
			my $last_position_in_loop = 1;
			warn "\tProcessing $contig...\n";

			# Foreach window
			WINDOW: for(my $i=0; $i<$$fasta_lengths{$contig}; $i+=$$HT_data{'max_phase_length'}) {
				my $stop_window = ($i + $$HT_data{'max_phase_length'});

				# Check both alignment and VCF files have been made
				my $VCF_file = "$$HT_data{'out_folder'}/$$HT_data{'VCF_filename'}-$contig-$i-$stop_window";
				my $alignment_file = "$$HT_data{'out_folder'}/$BAM_filename-$contig-$i-$stop_window";

				# Sometimes no positions in VCF or reads in BAM, so don't try and phase (but give warning)
				if(! -e $VCF_file) {
					warn "VCF file not found: $VCF_file. Either no positions in VCF in that window or need to re-run step 1\n";
					next WINDOW;
				}
				if(! -e $alignment_file) {
					warn "Alignment file not found: $alignment_file. Either no reads aligning in that window or need to re-run step 2\n";
					next WINDOW;
				}

				# Phase (not on grid)
				my $cmd = "perl $VCF_phase_script -a $alignment_file -v $VCF_file -s $sample_number\n";
				if($$HT_data{'run_grid'} eq 'n') { system($cmd); } 
				# Phase in parallel
				else {
					push @cmds, $cmd;
				}
				#die "end here\n";
			}
		}
	}
	# Run paralel cmds:
	if($$HT_data{'run_grid'} ne 'n') { &run_cmds_on_grid(\@cmds, 'cmds-HaplotypeTools-for-grid.txt', $$HT_data{'platform'}, $$HT_data{'queue'}); }
}

# Step 4: Validate and assign phase groups
if($$HT_data{'steps'} =~ m/4/) {
	warn "Step 4: Validate and assign phase groups...\n";
	
	my @cmds;
	foreach my $sample_number(keys %{$$HT_data{'BAM_filename'}}) {
		my @phased_files = <$$HT_data{'out_folder'}/*-phased-$sample_number>;
		foreach my $file(@phased_files) {
			warn "\t$file...\n";

			# Validate (not on grid)
			my $cmd = "perl $VCF_phase_validate_and_assign_script -v $file -c $$HT_data{'cutoff'} -m $$HT_data{'min_depth'}\n";
			if($$HT_data{'run_grid'} eq 'n') { system($cmd); }
			# Phase in parallel
			else {
				push @cmds, $cmd;
			}
		}
		# Run paralel cmds:
		if($$HT_data{'run_grid'} ne 'n') { &run_cmds_on_grid(\@cmds, 'cmds2-HaplotypeTools-for-grid.txt', $$HT_data{'platform'}, $$HT_data{'queue'}); }
	}
}

# Step 5: Concatenate
if($$HT_data{'steps'} =~ m/5/) {
	warn "Step 5: Concatenate...\n";

	# Header
	warn "\tcat $header_file...\n";
	my $cat1 = `cp $header_file $$HT_data{'out_VCF'}`;

	# Foreach contig
	warn "\tMerging individually phased samples...\n";
	CONTIGS: foreach my $contig(sort keys %{$fasta_lengths}) {
		my $last_position_in_loop = 1;
		warn "\tProcessing $contig...\n";

		# Foreach window
		WINDOW: for(my $i=0; $i<$$fasta_lengths{$contig}; $i+=$$HT_data{'max_phase_length'}) {
			my $stop_window = ($i + $$HT_data{'max_phase_length'});

			# VCF_file for most of the VCF info
			my $VCF_file = "$$HT_data{'out_folder'}/$$HT_data{'VCF_filename'}-$contig-$i-$stop_window";

			# Sometimes no positions in VCF so don't try and process (but give warning)
			if(! -e $VCF_file) {
				warn "VCF file not found: $VCF_file. Either no positions in VCF in that window or need to re-run step 1\n";
				next WINDOW;
			}

			# Save VCF to memory
			my $VCF_positions_to_line = vcflines::VCF_to_position_to_line_hash($VCF_file);

			# Replace ID and sample from individual phased files
			my %phased_lines;
			my @phased_files = <$VCF_file-phased-*-and-assigned.tab>;
			foreach my $file(@phased_files) {
				my ($contig, $start, $stop, $sample_number) = $file =~ /-(.+)-(\d+)-(\d+)-phased-(\d+)-and-assigned.tab$/;
				#my $sample_info = "sample_info$sample_number";
				#warn "phased file: $file ($contig, $start, $stop, $sample_number)\n";

				# Open each phased file updating ID and samples for saved VCF
				open my $fh, '<', $file or die "Cannot open $file : $!\n";
				while(my $line=<$fh>) {
					chomp $line;
					my @bits = split /\t/, $line;
					die "Error: $file not formatted correctly. Re-run steps 34\n" if(!defined $bits[3]);
					my ($contig, $pos, $id, $sample) = @bits;
					die "Error: pos $pos not found in $file. Re-run steps 34\n" if(!defined $$VCF_positions_to_line{$pos});
					if($id =~ m/phased/) { $phased_lines{$pos} = 1; }
					my $saved_line = $$VCF_positions_to_line{$pos};

					# Save only "sample-" or "phase_info" info (not read numbers)
					my @id_parts = split /\;/, $id;
					my $new_id = '';
					foreach my $id_part(@id_parts) {
						if($id_part =~ m/^sample|^phase_info/) { $new_id .= "$id_part;"; }
					}
					$id = $new_id;

					# Replace ID and sample
					my @bits2 = split /\t/, $saved_line;
					$bits2[2] .= $id;
					$bits2[($sample_number + 9)] = $sample;

					# Save
					my $new_line = join "\t", @bits2;
					#my $new_line = vcflines::VCF_line_make($VCF_line);
					$$VCF_positions_to_line{$pos} = $new_line;
				}
			}

			# update Format and non-phased samples on phased lines
			foreach my $position(keys %phased_lines) {
				my $VCF_line = $$VCF_positions_to_line{$position};
				#warn "pos $position and line $VCF_line (trying to get format with PID and all samples correct)\n";
				my @bits = split /\t/, $VCF_line;

				# Save samples that are phased
				my %phased_samples;
				my @id_parts = split /;/, $bits[2];
				foreach my $id_part(@id_parts) {
					next if($id_part !~ m/phased/);
					my @id_part_parts = split /-/, $id_part;
					$phased_samples{$id_part_parts[1]} = 1;
					#warn "sample phased = $id_part_parts[1]\n";
				}

				# Update format
				$bits[8] .= ":PID";
		
				# Update every sample
				my $number_of_samples = (scalar(@bits) - 9);
				for(my $i=9; $i<scalar(@bits); $i++) {
					my $sample_number = ($i - 9);
					if(!defined $phased_samples{$sample_number}) {
						$bits[$i] .= ":.";
						#my @sample_id_parts = split /:/, $bits[$i];
					}
				}

				# Save new line
				my $new_line = join "\t", @bits;
				#warn "new line = $new_line\n";
				$$VCF_positions_to_line{$position} = $new_line;
			}

			# Print new VCFs VCF_position_to_line
			my $VCF_windows_output = $VCF_file . '-phased';
			open my $ofh, '>', $VCF_windows_output or die "Cannot open $VCF_windows_output : $!";
			VCFPOSITIONS: foreach my $positions(sort { $a <=> $b } keys %{$VCF_positions_to_line}) {
				my $line = $$VCF_positions_to_line{$positions};
				print $ofh "$line\n";
			}
			close $ofh;
		}
	}

	# Merge VCF-phased files into single file
	warn "\tMerging individually phased files...\n";
	my %phased_file_order;
	my @phased_files = <$$HT_data{'out_folder'}/*-phased>;
	foreach my $file(sort { number_strip($a) <=> number_strip($b) } @phased_files) {
		my ($contig, $start, $stop) = $file =~ /-(.+)-(\d+)-(\d+)-phased$/;
		$phased_file_order{$contig} .= "$file\n";
	}
	foreach my $contig(sort keys %phased_file_order) {
		my $files = $phased_file_order{$contig};
		my @file = split /\n/, $files;
		foreach my $f(@file) {
			warn "\tcat $f...\n";
			my $cat2 = `cat $f >> $$HT_data{'out_VCF'}`;
		}
	}
	#warn "$header_file\n";

	warn "\tMerging summary files...\n";
	my %percent_agree_to_tally_summary;

	# Save content of every summary
	my @summary_files = <$$HT_data{'out_folder'}/*-pc_cutoff.tab>;
	foreach my $file(@summary_files) {

		# Open each summary file
		open my $fh, '<', $file or die "Cannot open $file : $!\n";
		while(my $line=<$fh>) {
			chomp $line;
			next if($line =~ m/^Percent/);
			my @bits = split /\t/, $line;
			#die "Error: $file not formatted correctly. Re-run steps 34\n" if(!defined $bits[1]);
			if(!defined $bits[1]) { warn "$file not formatted correctly with line $line.\n"; }
			$percent_agree_to_tally_summary{$bits[0]} += $bits[1];
		}
	}

	# Print
	open my $ofh, '>', $$HT_data{'out_summary'} or die "Cannot open $$HT_data{'out_summary'} : $!";
	print $ofh "Percent_reads_agreeing_haplotype_phase\tTally\n";
	foreach my $percent(sort { $a <=> $b } keys %percent_agree_to_tally_summary) {
		my $tally = $percent_agree_to_tally_summary{$percent};
		print $ofh "$percent\t$tally\n";
	}
	close $ofh;
}
warn "Finished\n";

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	die "Error, cmd $cmd died with return $ret\n" if($ret);
	return 1;
}

sub number_strip {
	my $line = shift;
	my ($start, $stop) = $line =~ /(\d+)-(\d+)-phased$/;
	return $start;
}

sub run_cmds_on_grid {
	my ($cmds, $cmds_file, $platform, $queue) = @_;
	open my $ofh, '>', $cmds_file or die "Cannot open $cmds_file : $!\n";
	print $ofh @{$cmds};
	close $ofh;
	warn "Running commands in parallel...\n";
	my $cmd_parallel = "$Run_Commands_python --platform $platform --queue $queue --mem 2 --throttle_nodes 95 --cmds_per_node 10 $cmds_file 2>&1 | tee $cmds_file.log";
	process_cmd($cmd_parallel);	
}

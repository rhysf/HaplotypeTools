#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use lib "/home/unix/rfarrer/perl5/lib/perl5/x86_64-linux-thread-multi/";
use Bio::SeqIO;
use Bio::DB::HTS;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_VCF;
use read_FASTA;
use VCF_phase;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -v <VCF> -b <BAM (sorted)> -f <reference fasta>
Optional: -c\tCut-off percent reads supporting phase group [90]
          -m\tMinimum read depth overlapping two heterozygous SNPs [4]
	  -r\tMax phase length (longer=slower) [10000]
	  -s\tSteps (1=process VCF, 2=process BAM, 3=phase, 4=validate and assign phase groups, 5=concatenate) [12345]\n
Parallel: -g\tRun commands on the grid (y/n) [n]
	  -a\tPlatform (UGER, LSF, GridEngine) [UGER]
	  -q\tQueue name [short]\n
Outputs:  -o\tOutput folder for tmp files [opt_v-HaplotypeTools-phased-r-opt_r]
	  -p\tPhased VCF [opt_v-Phased-m-opt_m-c-opt_c-r-opt_r.vcf]\n";
our($opt_a, $opt_b, $opt_c, $opt_f, $opt_g, $opt_l, $opt_m, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_v);
getopt('abcfglmopqrsv');
die $usage unless ($opt_v && $opt_b && $opt_f);
if(!defined $opt_a) { $opt_a = 'UGER'; }
if(!defined $opt_c) { $opt_c = 90; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_m) { $opt_m = 4; }
if(!defined $opt_r) { $opt_r = 10000; }
if(!defined $opt_s) { $opt_s = 12345; }
if(!defined $opt_q) { $opt_q = 'short'; }
if(!defined $opt_o) { $opt_o = $opt_v . '-HaplotypeTools-phased-r-' . $opt_r; }
if(!defined $opt_p) { $opt_p = "$opt_v-Phased-m-$opt_m-c-$opt_c-r-$opt_r.vcf"; }
if(!defined $opt_l) { $opt_l = $opt_v .'-Tally-agree-disagree-' . $opt_m . '-MD-' . $opt_c . '-pc_cutoff.tab'; }
die "-g is not n, N, y or Y: $opt_g\n" if($opt_g !~ m/^(n|N|y|Y)$/);
die "-a is not UGER, LSF or GridEngine: $opt_a\n" if($opt_a !~ m/^(UGER|LSF|GridEngine)$/);
$opt_c /= 100;

# Output folder
if(! -d $opt_o) { my $cmd1 = `mkdir $opt_o`; }

# Dependencies
my $VCF_phase_script = "$Bin/util/VCF_phase.pl";
my $VCF_phase_validate_and_assign_script = "$Bin/util/VCF_phased_validate_and_assign_phase_groups.pl";
my $Run_Commands_python = "$Bin/util/Run_cmds_on_grid.py";
foreach($VCF_phase_script, $VCF_phase_validate_and_assign_script, $Run_Commands_python) { die "Cannot find dependency script $_ : $!\n" if(! -e $_); }

# Save names and length of reference sequence
my ($fasta_lengths) = fastafile::fasta_id_to_seq_length_hash($opt_f);

# Step 1: Split VCF
my ($VCF, $VCF_header);
if($opt_s =~ m/1/) {
	warn "Step 1: Split VCF...\n";
	my $VCF_split = vcflines::VCF_split_for_phasing($opt_v, $opt_o, $opt_r, $fasta_lengths);
}

# Step 2: Split BAM
if($opt_s =~ m/2/) {
	warn "Step 2: Split BAM...\n";
	my $sam = Bio::DB::HTS->new(-bam => $opt_b, -fasta=> $opt_f);
	CONTIGS: foreach my $contig(sort keys %{$fasta_lengths}) {

		warn "\tProcessing $contig...\n";
		my $last_position_in_loop = 1;
		for(my $i=$opt_r; $i<$$fasta_lengths{$contig}; $i+=$opt_r) {
			my $loop_stop = ($i - 1);

			# For print windows
			my $start_window = $last_position_in_loop;
			if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }

			my $stop_window = ($i + $opt_r);
			#warn "Processing $contig $last_position_in_loop - $loop_stop\n";
			my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $loop_stop);
			my $lines = vcfphase::return_lines_from_BioDBHTS_array(\@alignments);

			# Print
			my $outfile = "$opt_o/$opt_b-$contig-$start_window-$i";
			open my $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
			if(defined $lines) { print $ofh $lines; }
			close $ofh;
			$last_position_in_loop = $i;
		}

		# Last group of reads for contig
		#warn "Processing $contig $last_position_in_loop - $$fasta_lengths{$contig}\n";
		my $start_window = $last_position_in_loop;
		if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }
		my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $$fasta_lengths{$contig});
		my $lines = vcfphase::return_lines_from_BioDBHTS_array(\@alignments);

		# Print
		my $stop_window = ($start_window + $opt_r);
		my $outfile = "$opt_o/$opt_b-$contig-$start_window-$stop_window";
		open my $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
		if(defined $lines) { print $ofh $lines; }
		close $ofh;
	}
}

# Step 3: Phase
if($opt_s =~ m/3/) {
	warn "Step 3: Phase...\n";
	my @cmds;
	CONTIGS: foreach my $contig(sort keys %{$fasta_lengths}) {
		my $last_position_in_loop = 1;
		warn "\tProcessing $contig...\n";
		#for(my $i=$opt_r; $i<$$fasta_lengths{$contig}; $i+=$opt_r) {
		for(my $i=0; $i<$$fasta_lengths{$contig}; $i+=$opt_r) {
			my $stop_window = ($i + $opt_r);

			# Check both alignment and VCF files have been made
			my $VCF_file = "$opt_o/$opt_v-$contig-$i-$stop_window";
			my $alignment_file = "$opt_o/$opt_b-$contig-$i-$stop_window";
			die "VCF file not found: $VCF_file. Re-run Step 1 : $!" if(! -e $VCF_file);
			die "Alignment file not found: $alignment_file. Re-run Step 2 : $!" if(! -e $alignment_file);

			# Phase (not on grid)
			my $cmd = "perl $VCF_phase_script -a $alignment_file -v $VCF_file -c $contig\n";
			if($opt_g eq 'n') { system($cmd); } 
			# Phase in parallel
			else {
				push @cmds, $cmd;
			}
		}
	}
	# Run paralel cmds:
	if($opt_g ne 'n') { &run_cmds_on_grid(\@cmds, 'cmds-HaplotypeTools-for-grid.txt'); }
}

# Step 4: Validate and assign phase groups
if($opt_s =~ m/4/) {
	warn "Step 4: Finalise...\n";
	my @cmds;
	my @phased_files = <$opt_o/*-phased>;
	foreach my $file(@phased_files) {
		warn "\t$file...\n";

		# Validate (not on grid)
		my $cmd = "perl $VCF_phase_validate_and_assign_script -v $file -c $opt_c -m $opt_m\n";
		if($opt_g eq 'n') { system($cmd); }
		# Phase in parallel
		else {
			push @cmds, $cmd;
		}
	}
	# Run paralel cmds:
	if($opt_g ne 'n') { &run_cmds_on_grid(\@cmds, 'cmds2-HaplotypeTools-for-grid.txt'); }
}

# Step 5: Concatenate
if($opt_s =~ m/5/) {
	warn "Step 5: Concatenate...\n";

	# Header
	my @header_files = <$opt_o/*-header>;
	my $header_file = $header_files[0];
	warn "\tcat $header_file...\n";
	my $cat1 = `cat $header_file > $opt_p`;

	# VCF-phased
	my %phased_file_order;
	my @phased_files = <$opt_o/*-assigned.VCF>;
	foreach my $file(sort { number_strip($a) <=> number_strip($b) } @phased_files) {
		my ($contig, $start, $stop) = $file =~ /-(.+)-(\d+)-(\d+)-phased-and-assigned.VCF$/;
		$phased_file_order{$contig} .= "$file\n";
	}
	foreach my $contig(sort keys %phased_file_order) {
		my $files = $phased_file_order{$contig};
		my @file = split /\n/, $files;
		foreach my $f(@file) {
			warn "\tcat $f...\n";
			my $cat2 = `cat $f >> $opt_p`;
		}
	}
	#warn "$header_file\n";
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
	my ($start, $stop) = $line =~ /(\d+)-(\d+)-phased-and-assigned.VCF$/;
	return $start;
}

sub run_cmds_on_grid {
	my ($cmds, $cmds_file) = @_;
	open my $ofh, '>', $cmds_file or die "Cannot open $cmds_file : $!\n";
	print $ofh @{$cmds};
	close $ofh;
	warn "Running commands in parallel...\n";
	my $cmd_parallel = "$Run_Commands_python --platform $opt_a --queue $opt_q --mem 2 --throttle_nodes 95 --cmds_per_node 10 $cmds_file 2>&1 | tee $cmds_file.log";
	process_cmd($cmd_parallel);	
}

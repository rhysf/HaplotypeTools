#!/usr/bin/perl -w
use strict;
use lib "/home/unix/rfarrer/perl5/lib/perl5/";
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_Tab;
use read_FASTA;
use read_Tree;
use sliding_windows;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -p <phased in any/all; PIA> -a <Haplotype file 1 FASTA> -b <Haplotype file 2 FASTA> -n <Name tab Type tab Location for consensus genomes>\n
Optional: -f Folder for output FASTA and Trees [HaplotypeTools_output]
          -c Clades/Lineages/Metadata to cluster isolates by. Format (Name tab Clades/Lineages/Metadata) []
          -t Location of FastTree [/seq/annotation/bio_tools/FastTree/current/FastTreeMP]
          -m Mimimum length of haplotypes [100]
          -v Verbose (y/n) [n]\n
Windows:  -d Calculate windows (y/n) [n]
          -l Reference FASTA []
          -z Window length [10000]\n
Output:   -f Output folder for genomic FASTA files and Trees [HaplotypeTools_output]
          -s Output details [HaplotypeTools_details]
          -u Output relatives summary [HaplotypeTools_relatives_summary]
          -e Output clade summary, if -c used [HaplotypeTools_clade_summary]
          -w Output windows [HaplotypeTools_windows]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_i, $opt_l, $opt_m, $opt_n, $opt_p, $opt_s, $opt_t, $opt_u, $opt_v, $opt_w, $opt_z);
getopt('abcdefilmnpstuvwz');
die $usage unless ($opt_a && $opt_b && $opt_n && $opt_p);
foreach($opt_a, $opt_b, $opt_n, $opt_p) { die "Cannot open $_ : $!\n" unless (-e $_); }
if(!defined $opt_d) { $opt_d = 'n'; }
if(!defined $opt_f) { $opt_f = "HaplotypeTools_output"; }
if(!defined $opt_i) { $opt_i = "name"; }
if(!defined $opt_m) { $opt_m = 100; }
if(!defined $opt_t) { $opt_t = "/seq/annotation/bio_tools/FastTree/current/FastTreeMP"; }
if(!defined $opt_s) { $opt_s = "HaplotypeTools_details"; }
if(!defined $opt_u) { $opt_u = "HaplotypeTools_relatives_summary"; }
if(!defined $opt_e) { $opt_e = "HaplotypeTools_clade_summary"; }
if(!defined $opt_w) { $opt_w = "HaplotypeTools_windows"; }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_z) { $opt_z = 10000; }
if($opt_d eq 'y') { 
	if(!defined $opt_l) { die "opt_d requires opt_l reference_fasta\n" unless(-e $opt_l); }
}

# Save Name tab Type Location
warn "Saving consensus genome locations...\n";
my $genomes_name_to_type_to_location = tabfile::save_columns_to_column_hash($opt_n, 0, 1, 2);
my $opt_c_data;
if($opt_c) { 
	die "Cannot open $opt_c : $!\n" unless(-e $opt_c); 
	warn "Saving clades/lineages/metadata to cluster by...\n";
	$opt_c_data = tabfile::save_columns_to_column_hash($opt_c, 0, 1);
}
#print Dumper($genomes_name_to_type_to_location);

# Save haplotypes
warn "Saving haplotypes...\n";
my $haplotypes1_struct = fastafile::fasta_to_struct($opt_a);
my $haplotypes2_struct = fastafile::fasta_to_struct($opt_b);

# Check they match
warn "Check haplotype names match in $opt_a and $opt_b...\n";
for(my $i=0; $i< scalar(@{$$haplotypes1_struct{'order'}}); $i++) {
	my $id1 = $$haplotypes1_struct{'order'}[$i];
	my $id2 = $$haplotypes2_struct{'order'}[$i];
	die "id $id1 in $opt_a does not have a match or in the same order in $opt_b ($id2). Consult README\n" if($id1 ne $id2);
}

# Go through PIA and create FASTA file
warn "Go through PIA and create folder of FASTA files in $opt_f...\n";
my $cmd1 = "mkdir $opt_f";
system($cmd1);
my %closest_relatives_to_haplotypes;
open my $fh, '<', $opt_p or die "Cannot open $opt_p : $!\n";
HAP: while(my $line=<$fh>) {
	chomp $line;
	my @bits = split /\t/, $line;
	my ($contig, $start, $stop) = @bits;

	# Check I have corresponding haplotypes
	my $haplotype_name = "$contig\_$start\_$stop";
	die "Cannot find $haplotype_name in $opt_a (as defined by PIA). Re-run pipeline\n" if(!defined $$haplotypes1_struct{'seq'}{$haplotype_name});
	die "Cannot find $haplotype_name in $opt_b (as defined by PIA). Re-run pipeline\n" if(!defined $$haplotypes2_struct{'seq'}{$haplotype_name});
	my $haplotype1 = $$haplotypes1_struct{'seq'}{$haplotype_name};
	my $haplotype2 = $$haplotypes2_struct{'seq'}{$haplotype_name};

	# Minimum length
	my $length = ($stop - $start);
	if($length < $opt_m) {
		if($opt_v eq 'y') { warn "Skipping $haplotype_name. Length $length < $opt_m.\n"; }
		next HAP;
	}

	# Outfile
	my $outfile = "$opt_f/$haplotype_name";
	my $treefile = "$outfile.tree";
	open my $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
	print $ofh ">haplotype1\n$haplotype1\n";
	print $ofh ">haplotype2\n$haplotype2\n";

	# Go through each genome extracting sequence
	foreach my $isolate(sort keys %{$genomes_name_to_type_to_location}) {
		die "no FASTA found" if(!defined $$genomes_name_to_type_to_location{$isolate}{'FASTA'});
		my $genomic_seq = &save_genomic_sequence($$genomes_name_to_type_to_location{$isolate}{'FASTA'}, $contig, $start, $stop, $opt_v);
		print $ofh ">$isolate\n$genomic_seq\n";
	
		# Check type
		if($opt_c) {
			die "No isolate $isolate found in $opt_c\n" if(!defined $$opt_c_data{$isolate});
		}
	}

	# Construct Tree
	my $cmd2 .= "$opt_t -nt ";
	eval {
		$cmd2 .= " $outfile > $treefile"; # sometimes this step fails if it's too large.
		&process_cmd($cmd2);
	};
	if ($@) { print "Error:  $@"; }

	# Save tree files
	my $tree_lengths = treefile::save_pairwise_distances_from_tree($treefile);
	#print Dumper($tree_lengths);

	# Find closest relative
	if($opt_v eq 'y') { warn "Finding closest relatives...\n"; }
	foreach my $haplotype(keys %{$tree_lengths}) {
		# e.g. chrD_64_1185.tree
		foreach my $isolate1(keys %{$$tree_lengths{$haplotype}}) {
			next if(($isolate1 eq 'haplotype1') || ($isolate1 eq 'haplotype2'));
			foreach my $isolate2(keys %{$$tree_lengths{$haplotype}{$isolate1}}) {
				next if(($isolate2 ne 'haplotype1') && ($isolate2 ne 'haplotype2'));

				# Check if smaller distance. Then save
				my $length = $$tree_lengths{$haplotype}{$isolate1}{$isolate2};
				if(!defined $closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'length'}) { 
					$closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'length'} = $length; 
					$closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'relative'} = $isolate1;
				} else {
					if($length < $closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'length'}) {
						$closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'length'} = $length;
						$closest_relatives_to_haplotypes{$haplotype}{$isolate2}{'relative'} = $isolate1;
					}
				}
			}
		}
	}
}
close $fh;

# Print summary
my %summary;
warn "print summary...\n";
open my $ofh1, '>', $opt_s or die "Cannot open $opt_s : $!\n";
if(defined $opt_c) { print $ofh1 "HaplotypeName\tContig\tStart\tStop\tHap1_relative\tHap2_relative\tHap1_clade_or_other\tHap2_clade_or_other\n"; }
else { print $ofh1 "HaplotypeName\tContig\tStart\tStop\tHap1_relative\tHap2_relative\n"; }
foreach my $haplotype(keys %closest_relatives_to_haplotypes) {

	# Contig, start, stop
	my @hap_parts = split /\_/, $haplotype;
	my $contig = '';
	for(my $i=0; $i<(scalar(@hap_parts) - 2); $i++) { $contig .= "$hap_parts[$i]_"; }
	$contig =~ s/\_$//;
	my $start = $hap_parts[(scalar(@hap_parts) - 2)];
	my $stop = $hap_parts[(scalar(@hap_parts) - 1)];
	$stop =~ s/\.tree//;

	# Relatives
	my $hap1_relative = $closest_relatives_to_haplotypes{$haplotype}{'haplotype1'}{'relative'};
	my $hap2_relative = $closest_relatives_to_haplotypes{$haplotype}{'haplotype2'}{'relative'};
	my ($hap1_clade, $hap2_clade) = ('NA', 'NA');
	if(defined $opt_c) {
		$hap1_clade = $$opt_c_data{$hap1_relative};
		$hap2_clade = $$opt_c_data{$hap2_relative};
	} 

	# summary
	$summary{'all'}{$hap1_relative} = $hap1_clade;
	$summary{'all'}{$hap2_relative} = $hap2_clade;
	$summary{'haplotype_relative'}{$hap1_relative} += ($stop - $start);
	$summary{'haplotype_relative'}{$hap2_relative} += ($stop - $start);
	$summary{'haplotype_clades'}{$hap1_clade} += ($stop - $start);
	$summary{'haplotype_clades'}{$hap2_clade} += ($stop - $start);
	$summary{'total'} += ($stop - $start) *2;

	# Print
	if(defined $opt_c) { print $ofh1 "$haplotype\t$contig\t$start\t$stop\t$hap1_relative\t$hap2_relative\t$hap1_clade\t$hap2_clade\n"; }
	else { print $ofh1 "$haplotype\t$contig\t$start\t$stop\t$hap1_relative\t$hap2_relative\n"; }
}
close $ofh1;
die "Finished without calculating windows (-d $opt_d)\n" if($opt_d eq 'n');

# Print summary of relatives
warn "print relatives summary...\n";
open my $ofh2, '>', $opt_u or die "Cannot open $opt_u : $!\n";
print $ofh2 "Haplotype relative\tClade\tTotal nt\thaplotype relatives (nt)\thaplotype relatives (%)\n";
foreach my $relative(sort keys %{$summary{'all'}}) {
	my $clade = $summary{'all'}{$relative};
	my $hap_relative = $summary{'haplotype_relative'}{$relative};
	my $total = $summary{'total'};
	my $hap_relative_pc = sprintf("%.2f", (($hap_relative / $total)) * 100);
	print $ofh2 "$relative\t$clade\t$total\t$hap_relative\t$hap_relative_pc\n";
}
close $ofh2;

# Print summary of clades
if($opt_c) {
	warn "Print clade summery...\n";
	open my $ofh3, '>', $opt_e or die "Cannot open $opt_e : $!\n";
	print $ofh3 "Clade\tTotal nt\thaplotype clades (nt)\thaplotype clades (%)\n";
	foreach my $clade(sort keys %{$summary{'haplotype_clades'}}) {
		my $total = $summary{'total'};
		my $hap_clades = $summary{'haplotype_clades'}{$clade};
		my $hap_clades_pc = sprintf("%.2f", (($hap_clades / $total)) * 100);
		print $ofh3 "$clade\t$total\t$hap_clades\t$hap_clades_pc\n";
	}
	close $ofh3;
}

# Windows
warn "calculating windows...\n";
my $sequences = fastafile::fasta_to_struct($opt_l);

# Make the windows
my $windows_want = "Hap1 Hap2";
my ($windows, $total_GC) = slidingwindows::make_windows_for_seq_struct($sequences, $opt_z, $opt_i, $windows_want);
open my $fh2, '<', $opt_s or die "Cannot open $opt_s : $!\n";
while(my $line=<$fh2>) {
	chomp $line;

	# Ignore header
	next if($line =~ m/^HaplotypeName\t/);

	my @bits = split /\t/, $line;
	die "Malformed $opt_s: $line\n" if(!defined $bits[5]);
	my ($hap_name, $contig, $start, $stop, $hap1_relative, $hap2_relative) = @bits;
	if($opt_c) {
		$hap1_relative = $bits[6];
		$hap2_relative = $bits[7];
	}

	# Find windows and fill (if multiple haplotypes are found per window - it will overwrite them. Also resolution is 10 bp)
	for(my $i=$start; $i< $stop; $i+=10) {
		my ($start_window, $stop_window) = slidingwindows::find_window_for_position($$windows{$opt_i}{$contig}, $contig, $i);
		$$windows{$opt_i}{$contig}{$start_window}{$stop_window}{'Hap1'} = $hap1_relative; 
		$$windows{$opt_i}{$contig}{$start_window}{$stop_window}{'Hap2'} = $hap2_relative; 
	}
}

# Print windows
warn "printing windows...\n";
open my $ofh3, '>', $opt_w or die "Cannot open $opt_w : $!\n";
print $ofh3 "Contig\tWindow_Start\tWindow_Stop\tRunning_Position\tHap1\tHap2\n";

my $running_window_length;
foreach my $contigs(@{$$sequences{'order'}}) {
	STARTS: foreach my $starts(sort {$a <=> $b} keys %{$$windows{$opt_i}{$contigs}}) {
		STOPS: foreach my $stops(sort {$a <=> $b} keys %{$$windows{$opt_i}{$contigs}{$starts}}) {
			my $Hap1 = $$windows{$opt_i}{$contigs}{$starts}{$stops}{'Hap1'};
			my $Hap2 = $$windows{$opt_i}{$contigs}{$starts}{$stops}{'Hap2'};

			$running_window_length += ($stops - $starts);
			print $ofh3 "$contigs\t$starts\t$stops\t$running_window_length\t$Hap1\t$Hap2\n";
		}
	}
}
close $ofh3;
warn "Finished.\n";

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	if ($ret) { die "Error, cmd: $cmd died with ret: $ret"; }
	return;
}

sub save_genomic_sequence {
	my ($file, $contig, $start, $stop, $verbose) = @_;
	my $genomic_seq;
	if($verbose eq 'y') { warn "save_genomic_sequence: Extract $contig\_$start\_$stop from $file...\n"; }
	my $inseq = Bio::SeqIO->new('-file' => "<$file",'-format' => 'fasta');
	FASTA: while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;

		next FASTA if($id ne $contig);
		$genomic_seq = substr $seq, ($start - 1), (($stop - $start) + 1);
		last FASTA;
	}
	die "No genomic sequence found for $contig $start - $stop in $file\n" if(!defined $genomic_seq);
	return $genomic_seq;
}

#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use read_Tab;
use read_FASTA;
use read_Tree;
use sliding_windows;
use Data::Dumper;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -l <Reference FASTA> -s <Input details (HaplotypeTools_details)>\n
Optional: -c Clades/Lineages/Metadata to cluster isolates by. Format (Name tab Clades/Lineages/Metadata) []
          -s Input details [HaplotypeTools_details]
          -t Location of FastTree [/seq/annotation/bio_tools/FastTree/current/FastTreeMP]
          -v Verbose (y/n) [n]\n
Windows:  -z Window length [10000]\n
Output:   -f Output folder for genomic FASTA files and Trees [HaplotypeTools_output]
          -w Output windows [HaplotypeTools_windows]\n";
our($opt_c, $opt_i, $opt_l, $opt_s, $opt_v, $opt_w, $opt_z);
getopt('cilpsvwz');
die $usage unless ($opt_l && $opt_s);
foreach($opt_l) { die "Cannot open $_ : $!\n" unless (-e $_); }
if(!defined $opt_i) { $opt_i = "name"; }
if(!defined $opt_s) { $opt_s = "HaplotypeTools_details"; }
if(!defined $opt_w) { $opt_w = "HaplotypeTools_windows"; }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_z) { $opt_z = 10000; }

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

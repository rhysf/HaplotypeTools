package bamlines;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use lib "/home/unix/rfarrer/perl5/lib/perl5/x86_64-linux-thread-multi/";
use Bio::DB::HTS;
use Data::Dumper;

### rfarrer@broadinstitute.org

sub BAM_split_for_phasing {
	my ($HT_data, $fasta_lengths) = @_;

	# Input
	my $fasta = $$HT_data{'ref_fasta'};
	my $outfolder = $$HT_data{'out_folder'};
	my $haplotype_length = $$HT_data{'max_phase_length'};
	foreach($fasta, $$HT_data{'out_folder'}, $$HT_data{'max_phase_length'}) { die "VCF_split_for_phasing: HT_data not properly initialised: $!\n" if(!defined $_); }

	# Foreach BAM
	foreach my $sample_number(keys %{$$HT_data{'BAM_file'}}) {
		my $BAM = $$HT_data{'BAM_file'}{$sample_number};
		my $BAM_filename = $$HT_data{'BAM_filename'}{$sample_number};
		warn "\t$BAM ->\n";
		my $sam = Bio::DB::HTS->new(-bam => $BAM, -fasta=> $fasta);

		# Foreach contig
		CONTIGS: foreach my $contig(sort keys %{$fasta_lengths}) {
			warn "\tProcessing $contig...\n";

			# Foreach non-overlapping window
			my $last_position_in_loop = 1;
			for(my $i=$haplotype_length; $i<$$fasta_lengths{$contig}; $i+=$haplotype_length) {
				my $loop_stop = ($i - 1);

				# For print windows
				my $start_window = $last_position_in_loop;
				if($last_position_in_loop eq 1) { $start_window = ($last_position_in_loop - 1); }

				my $stop_window = ($i + $haplotype_length);
				#warn "Processing $contig $last_position_in_loop - $loop_stop\n";
				my @alignments = $sam->get_features_by_location(-seq_id =>$contig,-start => $last_position_in_loop,-end => $loop_stop);
				my $lines = &return_lines_from_BioDBHTS_array(\@alignments);

				# Print
				my $outfile = "$outfolder/$BAM_filename-$contig-$start_window-$i";
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
			my $lines = &return_lines_from_BioDBHTS_array(\@alignments);

			# Print
			my $stop_window = ($start_window + $haplotype_length);
			my $outfile = "$outfolder/$BAM_filename-$contig-$start_window-$stop_window";
			open my $ofh, '>', $outfile or die "Cannot open $outfile : $!\n";
			if(defined $lines) { print $ofh $lines; }
			close $ofh;
		}
	}
	return 1;
}

sub return_lines_from_BioDBHTS_array {
	my ($alignments) = $_[0];
	my $lines;
	SAMSEQS: for my $a(@{$alignments}) {
		# where does the alignment start in the reference sequence
		my $start = $a->start;
		my $end = $a->end;
		# query/read sequence
		my $query_dna = $a->query->dna; 
		# where does the alignment start in the query sequence
		my $query_start = $a->query->start;     
		my $query_end   = $a->query->end;
		# new for determining indels
		my $cigar = $a->cigar_str;

		# Ignore reads specifying/aligning over indels
		#my $read_length = ($query_end - $query_start);
		#my $ref_seq_length = ($end - $start);
		#next SAMSEQS if ($read_length ne $ref_seq_length);

		# Save
		my $line = join "\t", $start, $end, $query_dna, $query_start, $query_end, $cigar;
		$lines .= "$line\n";
	}
	return $lines;
}

#### local subroutine
sub check_for_indel_in_cigar {
	my ($cigar, $pos) = @_;
	my @cigar_parts = split (/(\d+)/, $cigar);
	my $cumulative_position = 0;
	my $indel_in_cigar = 'n';
	my $softclipping_start = 0;
	foreach my $cigar_part(@cigar_parts) {
		next if($cigar_part eq '');
		#warn "found $cigar_part\n";
		if($cigar_part =~ m/\d+/) { $cumulative_position += $cigar_part; }
		else {
			# found insertion of deletion? (+1 for 0 based dna coord, and +1 because the indel is listed after the number)
			if($cumulative_position eq ($pos + 2)) {
				if($cigar_part =~ m/I|D/) { $indel_in_cigar = 'y'; }
				#warn "mutation = $cigar_part\n";
			}

			# softclipping start.
			if($cigar_part eq 'S') {
				if($softclipping_start eq 0) { 
					$softclipping_start = 1;
					$cumulative_position = 0;
				}
			}
		}
		#warn "cumulative position = $cumulative_position\n";
	}
	return $indel_in_cigar;
}

sub count_softclipping_start {
	my $cigar = $_[0];
	my $count = 0;

	my @cigar_parts = split (/(\d+)/, $cigar);
	my $cumulative_position = 0;
	my $indel_in_cigar = 'n';
	foreach my $cigar_part(@cigar_parts) {
		next if($cigar_part eq '');
		#warn "found $cigar_part\n";
		if($cigar_part =~ m/\d+/) { $cumulative_position += $cigar_part; }
		else {
			if($cigar_part eq 'S') { $count = $cumulative_position; }
			last;
		}
	}
	return $count;
}

1;

#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 <haplotypes list> <full sequence file>\n";
die $usage if (@ARGV !=2);
my ($genes_of_interest, $sequence_file) = @ARGV;

# Count haplotypes per contig
my $contigs_found = &save_contigs_to_hash($genes_of_interest);

# Haplotypes per contig
warn "Haplotypes per contig:\n";
foreach my $contig(keys %{$contigs_found}) {
	my $tally = $$contigs_found{$contig};
	print "$contig\t$tally\n";
}

# Go through FASTA and if found - save info
my %found_per_contig_per_kb;
my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",'-format' => 'fasta' ) ;
while (my $seq_obj = $inseq->next_seq ) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my $desc = $seq_obj->description;
	
	if(defined $$contigs_found{$id}) { $found_per_contig_per_kb{$id} = ($$contigs_found{$id} / length($seq)); }
}
warn "greatest by kb:\n";
foreach my $contig(keys %found_per_contig_per_kb) {
	my $perkb = $found_per_contig_per_kb{$contig};
	print "$contig\t$perkb\n";
}

sub save_contigs_to_hash {
	my $file = $_[0];
	my %found;
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;	
		my @bits = split /\t/, $line;
		my ($contig, $group) = @bits;
		$found{$contig}++;
	}
	close $fh;
	return \%found;
}

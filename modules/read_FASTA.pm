package fastafile;
use strict;
use Bio::SeqIO;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub fasta_to_struct {
	my $input = $_[0];
	my %struct;
	$struct{'filename'} = $input;
	warn "fasta_to_struct: saving from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;

		# Save
		$struct{'seq'}{$id} = $seq;
		$struct{'desc'}{$id} = $desc;
		push @{$struct{'order'}}, $id;
	}
	return \%struct;
}

sub fasta_id_to_seq_hash {
	my $input = $_[0];
	my (%sequences, %descriptions);
	my @order;
	warn "fasta_id_to_seq_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;
		$sequences{$id}=$seq;
		$descriptions{$id}=$desc;
		push @order, $id;
	}
	return (\%sequences, \%descriptions, \@order);
}

sub fasta_id_to_seq_length_hash {
	my $input = $_[0];
	my (%lengths);
	warn "fasta_id_to_seq_length_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    		my $seq = $seq_obj->seq;
		my $length = length($seq);
		$lengths{$id} = $length;
	}
	return (\%lengths);
}

sub fasta_hash_print_simple_outfile {
	my ($fasta_hash, $outfile) = @_;
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!\n";
	FASTA: foreach my $id(keys %{$fasta_hash}) {
		my $seq = $$fasta_hash{$id};
		$seq =~ s/(\S{60})/$1\n/g;
		print $ofh ">$id\n$seq\n";
	}
	close $ofh;
	return;
}


1;

package genomearray;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);

### rfarrer@broadinstitute.org

sub make_genome_hash_array_from_seq_hash {
	my $seq_hash = $_[0];
	my %genome;
	warn "filling genome array...\n";
	foreach my $ids(keys %{$seq_hash}) {
		my $length = length($$seq_hash{$ids});
		my @array;
		for(my $i=0; $i<=$length; $i++) { $array[$i] = 0; }
		push @{$genome{$ids}}, @array;
	}
	return \%genome;
}

1;

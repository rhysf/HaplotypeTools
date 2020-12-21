package tabfile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub save_columns_to_one_hash {
	my $file = $_[0];
	my $column_count = (scalar(@_) - 1);
	die "save_columns_to_one_hash: less than 1 column specified\n" if($column_count < 1);

	# Save from file
	my %interest;
	my $line_count = 0;
	warn "save_columns_to_one_hash: saving $column_count columns from $file...\n";
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	LIST: while(my $line=<$fh>) {
		chomp $line;
		next LIST if($line =~ m/^\n/);
		next LIST if($line eq '');

		# Save wanted columns to array
		my @save = &save_wanted_columns_to_array($line, $line_count, \@_);

		# Save wanted columns to hash
		if($column_count eq 1)    { $interest{$save[0]} = 1; }
		elsif($column_count eq 2) { $interest{$save[0]}{$save[1]} = 1; }
		elsif($column_count eq 3) { $interest{$save[0]}{$save[1]}{$save[2]} = 1; }
		elsif($column_count eq 4) { $interest{$save[0]}{$save[1]}{$save[2]}{$save[3]} = 1; }
		elsif($column_count eq 5) { $interest{$save[0]}{$save[1]}{$save[2]}{$save[3]}{$save[4]} = 1; }
		else { die "Not coded for column_count > 5 ($column_count). This may cause problems. Consider using Deep::Hash instead\n"; }
		$line_count++;
	}
	close $fh;
	warn "$line_count lines in $file\n";
	warn scalar(keys(%interest)) . " unique whole entries saved\n";
	return(\%interest);
}

sub save_columns_to_column_hash {
	my $file = $_[0];
	my $column_count = (scalar(@_) - 1);
	die "save_columns_to_column_hash: less than 2 columns specified\n" if($column_count < 2);

	# Save from file
	my %interest;
	my $line_count = 0;
	warn "save_columns_to_column_hash: saving $column_count columns from $file...\n";
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	LIST: while(my $line=<$fh>) {
		chomp $line;

		# Save wanted columns to array
		my @save = &save_wanted_columns_to_array($line, $line_count, \@_);

		# Save wanted columns to hash
		if($column_count eq 2)    { $interest{$save[0]} = $save[1]; }
		elsif($column_count eq 3) { $interest{$save[0]}{$save[1]} = $save[2]; }
		elsif($column_count eq 4) { $interest{$save[0]}{$save[1]}{$save[2]} = $save[3]; }
		elsif($column_count eq 5) { $interest{$save[0]}{$save[1]}{$save[2]}{$save[3]} = $save[4]; }
		else { die "Not coded for column_count > 5 ($column_count). This may cause problems. Consider using Deep::Hash instead\n"; }
		$line_count++;
	}
	close $fh;
	warn "$line_count lines in $file\n";
	warn scalar(keys(%interest)) . " unique whole entries saved\n";
	return(\%interest);
}

# Local sub routine
sub save_wanted_columns_to_array {
	my ($line, $line_count, $values) = @_;
	my @bits = split /\t/, $line;
	my @wanted_columns;
	for(my $i=1; $i<scalar(@{$values}); $i++) {
		my $column_wanted = $$values[$i];
		if($line_count eq 0) { warn "\tColumn $column_wanted wanted\n"; }
		die "Undefined columns $column_wanted: $line\n" if(!defined $bits[$column_wanted]);
		push @wanted_columns, $bits[$column_wanted];
	}
	return @wanted_columns;
}

1;

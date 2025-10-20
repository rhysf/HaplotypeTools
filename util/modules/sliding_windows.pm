package slidingwindows;
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
use lib "$Bin/";
use read_Tab;
use Data::Dumper;

### rfarrer@broadinstitute.org

sub make_windows_for_seq_struct {
	my ($seq_struct, $window_length, $isolatename, $windows_interest) = @_;
	my @windows_of_interest = split /\s/, $windows_interest;
	#my @windows_of_interest = qw(anything HOMS HETS INDEL PHASE TOTALDEPTH AVERAGEDEPTH AVERAGEDEPTHEWL NORMALISEDDEPTH NORMALISEDDEPTHEWL NORMALISEDDEPTHGC NORMALISEDDEPTHGCEWL TAB);

	my (%windows, %total_GC);
	warn "make_windows_for_seq_hash: making $window_length windows according to seq hash...\n";
	foreach my $id(keys %{$$seq_struct{'seq'}}) {
		my $seq = $$seq_struct{'seq'}{$id};
		my $length = length($seq);

		#warn "Windows for $id = 1 - $length\n";
		my $previous_length = 0;

		# Along the contig
		for(my $i=$window_length; $i < $length; $i+=$window_length) {
			# warn "$previous_length - $i\n";	
			my $seq_window = substr $seq, $previous_length, ($i - $previous_length);
			my ($EWL, $GC_PC) = &get_effective_window_length_and_gcpc($seq_window);
			$total_GC{$isolatename}{$GC_PC}{'TOTALLENGTH'} += $window_length;
			$total_GC{$isolatename}{$GC_PC}{'TOTALLENGTHEWL'} += $EWL;
			$windows{$isolatename}{$id}{$previous_length}{$i}{'EWL'} = $EWL;
			$windows{$isolatename}{$id}{$previous_length}{$i}{'GC_PC'} = $GC_PC;
			foreach my $window_type(@windows_of_interest) { $windows{$isolatename}{$id}{$previous_length}{$i}{$window_type} = 0; }
			$previous_length = $i;
		}

		# End of contig	
		if($previous_length < $length) {
			my $seq_window = substr $seq, $previous_length, $length;
			my ($EWL, $GC_PC) = &get_effective_window_length_and_gcpc($seq_window);
			$total_GC{$isolatename}{$GC_PC}{'TOTALLENGTH'} += $window_length;
			$total_GC{$isolatename}{$GC_PC}{'TOTALLENGTHEWL'} += $EWL;
			$windows{$isolatename}{$id}{$previous_length}{$length}{'EWL'} = $EWL;
			$windows{$isolatename}{$id}{$previous_length}{$length}{'GC_PC'} = $GC_PC;

			foreach my $window_type(@windows_of_interest) { $windows{$isolatename}{$id}{$previous_length}{$length}{$window_type} = 0; }
		}
	}
	warn scalar(keys(%windows)) . " isolates with empty windows initalised.\n";
	return(\%windows, \%total_GC);
}

sub get_effective_window_length_and_gcpc {
	my ($seq) = $_[0];
	# EWL
	my $Ns = ($seq =~ tr/N|n//);
	my $EWL = (length($seq) - $Ns);
	# GC%
	my $CGs = ($seq =~ tr/C|c|G|g//);
	my $GC_PC = 0;
	if($EWL ne 0) { $GC_PC = ($CGs / $EWL)*100; }
	$GC_PC = sprintf("%.0f", $GC_PC);
	return($EWL, $GC_PC);
}

sub find_window_for_position {
	my ($windows_df, $supercontig, $position) = @_;
	my ($start_winds, $stop_winds);

	FINDWINDOW: foreach my $starts(sort {$a <=> $b} keys %{$windows_df}) {
		next FINDWINDOW if($position < $starts);

		# Potential stop 
		STOPS: foreach my $stop(sort {$a <=> $b} keys %{$$windows_df{$starts}}) {
			next STOPS if($position > $stop);
			$start_winds = $starts;
			$stop_winds = $stop;
			last FINDWINDOW;
		}
	}

	# Check window has been found
	if((!defined $start_winds) || (!defined $stop_winds)) {
		die "Could not find a window for position $position in $supercontig Check input files\n";
	}	
	return ($start_winds, $stop_winds);
}

sub save_dataframe_one_column {
	my ($windows_df, $feature_of_interest) = @_;
	my %info;
	my $position_of_interest = 0;	
	open WIN, "<$windows_df" or die "Cannot open $windows_df: $!\n";
	WINDOWS: while(my $line=<WIN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		if ($line =~ m/^Contig/) {
			for(my $i=0; $i<scalar(@bits); $i++) {
				if($bits[$i] eq $feature_of_interest) { $position_of_interest = $i; }
			}
			next WINDOWS;
		}
		my ($contig, $start, $stop, $running_position) = @bits;
		#$info{$isolate}{$stop} = $bits[$position_of_interest];
		$info{$contig}{$stop} = $bits[$position_of_interest];

	}
	close WIN;
	return \%info;
}

sub print_contig_borders {
	my ($cumulative_contig_lengths, $fh) = @_;
	
	# Contig borders	
	print $fh "abline(v=0, col=\"blue\", lwd=1.5, lty=2)\n";
	foreach my $contigs(keys %{$cumulative_contig_lengths}) {
		my $length = $$cumulative_contig_lengths{$contigs};
		print $fh "abline(v=$length, col=\"blue\", lwd=1.5, lty=2)\n";
	}
	return 1;
}

sub save_dataframe_specs {
	my $window_df = $_[0];

	# Save these aspects of the data
	my ($length);
	my (@contig_order);
	my (%contig_lengths, %cumulative_contig_lengths, %contig_half_way, %contig_start);

	# Save sum window length, contig lengths, cumulative lengths and contig order
	open WIN, "<$window_df" or die "Cannot open $window_df: $!\n";
	WINDOWS: while(my $line=<WIN>) {
		chomp $line;
		next WINDOWS if ($line =~ m/^Contig/);
		my @bits = split /\t/, $line;
		my ($contig, $start, $stop, $running_position) = @bits;
		$length = $running_position;
		$contig_lengths{$contig} = $stop;
		$cumulative_contig_lengths{$contig} = $running_position;

		if(!defined $contig_order[0]) { 
			$contig_order[0] = $contig; 
			$contig_start{$contig} = 0;
		}
		if(!defined $contig_start{$contig}) { $contig_start{$contig} = $running_position; }
		if($contig_order[(scalar(@contig_order) - 1)] !~ m/$contig/) { push @contig_order, $contig; }
	}
	close WIN;

	# Calculate half way points
	my $temp_length = 0;
	foreach(@contig_order) {
		my $half_contig_length = ($contig_lengths{$_} / 2);
		my $cumulative_length = $cumulative_contig_lengths{$_};
		$contig_half_way{$_} = ($cumulative_length - $half_contig_length);
	}
	#warn "end of df specs. length = $length\nContig order = @contig_order\n";
	return ($length, \@contig_order, \%cumulative_contig_lengths, \%contig_half_way, \%contig_start);
}



1;

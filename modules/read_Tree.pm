package treefile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;
use File::Glob;
use Hash::Merge qw(merge);
use List::Util qw(max);
use Bio::TreeIO;

### rfarrer@broadinstitute.org

sub newick_struct_to_pairwise_differences {
	my $newick_struct = $_[0];
	my %pairwise_differences;

	# Go through tree
	warn "newick_struct_to_pairwise_distances : saving pairwise distances...\n";
	while(my $tree = $newick_struct->next_tree) {

		# Save nodes
		my @nodes;
		foreach my $node ( $tree->get_nodes ) { push @nodes, $node; }

		# node count (1..n)
		my $node_count = scalar($tree->get_nodes);
		#print "node count = $node_count\n";

		# Pairwise
		for(my $i=0; $i<$node_count; $i++) {
			for(my $j=0; $j<$node_count; $j++) {
				next if($j <= $i);
				#warn "compare $i and $j\n";
				my $node_id1 = $nodes[$i]->id;
				my $node_id2 = $nodes[$j]->id;
				next unless((defined $node_id1) && (defined $node_id2));
				next unless(($node_id1 =~ m/^([a-zA-Z0-9]+)$/) && ($node_id2 =~ m/^([a-zA-Z0-9]+)$/)); # no . allowed in gene names!
				#warn "$node_id1 vs $node_id2\n";
				my $distance = $tree->distance(-nodes => [$nodes[$i], $nodes[$j]]);
				#warn "distance between $i) $node_id1 to $j) $node_id2 = $distance\n";
				$pairwise_differences{$node_id1}{$node_id2} = $distance;
				$pairwise_differences{$node_id2}{$node_id1} = $distance;
			}
		}
	}
	return \%pairwise_differences;
}

sub save_pairwise_distances_from_tree {
	my ($treefile) = $_[0];
	my %trees_struct;
	warn "save_pairwise_distances_from_tree: $treefile. Saving filename_prefix -> seq_ID1 -> seq_ID2 -> tree_length hash...\n";

	# Save trees and cluster numbers
	my @tree_parts = split /\//, $treefile;
	my $filename = $tree_parts[scalar(@tree_parts) - 1];
	my $cluster = $filename;

	# Save newick tree.
	warn "Save tree = $treefile, Cluster = $cluster\n";
	my $newick = new Bio::TreeIO(-file   => $treefile, -format => "newick");

	# Save Newick struct with info
	my $pairwise_distances = &newick_struct_to_pairwise_differences($newick);

	# Save pairwise distances by cluster
	my %pairwise_distances_with_cluster_name;
	$pairwise_distances_with_cluster_name{$cluster} = $pairwise_distances;
	%trees_struct = %{ merge( \%trees_struct, \%pairwise_distances_with_cluster_name ) };

	return \%trees_struct;
}

1;

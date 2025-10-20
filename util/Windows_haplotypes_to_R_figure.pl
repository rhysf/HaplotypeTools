#!/usr/bin/perl -w
use strict;
use Cwd;
my $dir = getcwd;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/modules";
use sliding_windows;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage: perl $0 -w <Windows haplotypes dataframe>\n
Optional: -s Scaling factor for width (nucleotides per mm) [20000]
          -h Scaling factor for height (genome per mm) [70]
          -r Resolution (The nominal resolution in ppi) [100]
          -c Minimum contig length to include label [1000000]
          -g Optional size of axis text [2]
          -l Legend size [2]
          -x X-min [0]\n
Outputs:  -o Output R plot [windows_plot_for_Haplotypes.R]\n";
our($opt_c, $opt_g, $opt_h, $opt_l, $opt_o, $opt_r, $opt_s, $opt_w, $opt_x);
getopt('cghlorswxz');
die $usage unless ($opt_w);
if(!defined $opt_c) { $opt_c = 1000000; }
if(!defined $opt_s) { $opt_s = 20000; }
if(!defined $opt_h) { $opt_h = 70; }
if(!defined $opt_r) { $opt_r = 100; }
if(!defined $opt_l) { $opt_l = 2; }
if(!defined $opt_g) { $opt_g = 2; }
if(!defined $opt_o) { $opt_o = 'windows_plot-for-Haplotypes.R'; }
if(!defined $opt_x) { $opt_x = 0; }
foreach($opt_s, $opt_h, $opt_r) { die "$_ needs to be entirely numerical\n" if ($_ !~ m/^\d+$/); }

# Count number of different relatives
my ($haplotype_keys, $haplotype_string) = &count_different_relatives_and_assign_numbers($opt_w);
my $number_of_different_species = scalar(keys(%{$haplotype_keys}));
warn "Number of different species = $number_of_different_species\n";

# Replace Hap1 and Hap2 column with those numbers
my $tmp_windows = ($opt_w . '-tmp');
open my $ofh_tmp, '>', $tmp_windows or die "Cannot open $tmp_windows : $!";
open my $fh, '<', $opt_w or die "Cannot open $opt_w : $!";
WINDOWS: while(my $line=<$fh>) {
	chomp $line;

	# print header
	if($line =~ m/Contig\tWindow_Start/) {
		print $ofh_tmp "$line\n";
		next WINDOWS;
	}

	# print rest
	my @bits = split /\t/, $line;
	my ($contig, $start, $stop, $running_pos, $hap1, $hap2) = @bits;
	my ($hap1_number, $hap2_number) = (0, 0);
	if(defined $$haplotype_keys{$hap1}) { $hap1_number = $$haplotype_keys{$hap1}; }
	if(defined $$haplotype_keys{$hap2}) { $hap2_number = $$haplotype_keys{$hap2}; }
	my $new_line = join "\t", $contig, $start, $stop, $running_pos, $hap1_number, $hap2_number;
	print $ofh_tmp "$new_line\n";
}
close $fh;
close $ofh_tmp;

# Width specified from 2nd column, end of 1st file.
my ($length_of_windows_dataframe, $order_of_contigs, $cumulative_contig_lengths, $contig_half_way_lengths, $contig_cumulative_starts) = slidingwindows::save_dataframe_specs($opt_w);

# Specifications for .pdf
my $width = ($length_of_windows_dataframe / $opt_s);
$width = sprintf("%.0f", $width);	
my $height = (400 + (2 * $opt_h));

# Make R file
my $R_file = ($dir . '/' . 'windows_plot-for-Haplotypes.R');
warn "Printing to $opt_o\n";
open my $ofh, '>', $opt_o or die "Cannot open $opt_o: $!\n";
print $ofh "require(graphics)\n";
print $ofh "library(plyr)\n";
#print $ofh "library(ggplot2)\n";
print $ofh "library(RColorBrewer)\n";
my $height_inches = ($height / 87);
my $width_inches = ($width / 87);
print $ofh "pdf(\"$R_file.pdf\", height=$height_inches, width=$width_inches)\n";

# par mar (margin area) and oma (outer margin area) = distance S, W, N, E of fig
# las=1 = Labels horizontal
print $ofh "par(mfrow=c(3,1), mar=c(1, 5, 0, 5), oma=c(15, 3, 7, 3), las=1, cex.axis=$opt_g)\n";

# Attach Hap1 (there is also Hap2)
my $name = 'Hap1';
my $window = $opt_w;
die "window $window (name $name) not defined\n" if (!defined $window);
print $ofh "sample <- read.table(\"$tmp_windows\", header=T, sep=\"\\t\")\n";
print $ofh "attach(sample)\n";
print $ofh "names(sample)\n";

# Colours
print $ofh "cols <- brewer.pal($number_of_different_species, \"Paired\")\n";
print $ofh "cols\n";

# legend
print $ofh "plot.new()\n";
print $ofh "legend(\"center\", c($haplotype_string), xpd=TRUE, horiz=TRUE, inset=c(0,0), bty=\"n\", cex=$opt_l, fill=cols)\n";

foreach my $haplotype(qw(Hap1 Hap2)) {

	# Colour values
	print $ofh "sample\$Col <- cols[as.numeric(cut(sample\$$haplotype,breaks = $number_of_different_species))]\n";

	# Plot
	print $ofh "plot(sample\$Running_Position, (sample\$$haplotype), type=\"h\", ylab=\"\", yaxt=\"n\", xlim=c($opt_x, $length_of_windows_dataframe), xaxt = \"n\", ylim=c(0,1), cex.lab=2, xaxs = \"r\", col=sample\$Col)\n";

	# Axes on side 
	my $position_of_name = (1 / 2);
	print $ofh "axis(2, at=$position_of_name, lwd.ticks=0, lab=c(\"$haplotype\"))\n";
	
	# Box around the current plot
	print $ofh "box(which = \"plot\", lty = \"solid\", lwd=1.5)\n";

	# Top x-axis and right y-axis	
	if($haplotype eq 'Hap1') { 
		# Gives x axis at top of graph, pos=how high up above axis (ONLY ONCE)
		my $x_axis_pos = (1 + (1 / 10));
		my $second_y_lab_pos = ($length_of_windows_dataframe + ($length_of_windows_dataframe / 20));
		print $ofh "lablist.top<-as.vector(c(\"0\", \"6000000\", \"12000000\", \"18000000\", \"24000000\", \"30000000\", \"36000000\", \"42000000\",\"48000000\",\"54000000\",\"60000000\"))\n";
		print $ofh "axis(3, at=lablist.top, lab=c(\"0\",\"6\",\"12\",\"18\",\"24\",\"30\",\"36\",\"42\",\"48\",\"54\",\"60\"), pos=$x_axis_pos, lwd=1.5, cex=2)\n";
		print $ofh "mtext(\"Position in genome (Mb)\", outer = TRUE, side=3, line=1, cex=2)\n";
		print $ofh "pointsiwant <- c(0, 1)\n";
		#print $ofh "axis(4, at = pointsiwant, lab = c(\"0\", \"1\"), tick = TRUE, pos=$second_y_lab_pos, lwd=1.5)\n";
		#print $ofh "axis(4, at=$position_of_name, lwd.ticks=0, lab=c(\"$opt_z\"), las=2, line=3)\n";
		#print $ofh "axis(4, at = 0, labels=FALSE, pos=$length_of_windows_dataframe, lwd=1)\n";
	}
	# Bottom axis
	if($haplotype eq 'Hap2') { 
		my $bottom_x_axis_pos = (0 - (1 / 10));
		my $bottom_x_axis_lab = $bottom_x_axis_pos - 0.2;
	
		# Make vector of locations for, and names of contigs
		my ($vector_start, $vector_half_way_points, $vector_labels);
		my $last_length = 0;
		foreach my $temp_vector_label(@{$order_of_contigs}) {
			my $half_way_points = $$contig_half_way_lengths{$temp_vector_label};

			# Enough space for contig label?
			#warn "Enough space? last length = $last_length... new half way = $half_way_points\n";
			my $contig_length = ($half_way_points - $last_length);
			#warn "$temp_vector_label goes on $contig_length\n";
			next if($contig_length < $opt_c);
			$last_length = $half_way_points;

			# Plot
			$vector_half_way_points .= "\"$half_way_points\",";
			if(!defined $vector_start) { $vector_start = "\"0\", \"$half_way_points\""; }
			$vector_labels .= ("\"$temp_vector_label\",");
		}
		$vector_half_way_points =~ s/\,$//;
		$vector_labels =~ s/\,$//;
		print $ofh "labpos.y<-as.vector(c($vector_half_way_points))\n";
		print $ofh "labpos.long<-as.vector(c($vector_start))\n";	
		print $ofh "axis(1, at=labpos.long, lwd.ticks=0, pos=$bottom_x_axis_pos, lwd=1.5, labels=FALSE)\n";
		print $ofh "axis(1, at=labpos.y, lwd.ticks=0, pos=$bottom_x_axis_pos, lwd=1.5, labels=FALSE)\n";
		print $ofh "axis(1, at=labpos.y, lab=c($vector_labels), pos=$bottom_x_axis_lab, lwd=0, las=2)\n";
		#print $ofh "mtext((expression(paste(\"Supercontig (\",n^o,\")\"))), outer = TRUE, side=1, line=8, cex=3)\n";
	}

	# Contig borders	
	print $ofh "abline(v=0, col=\"black\", lwd=1.5, lty=1)\n";
	foreach my $contigs(keys %{$cumulative_contig_lengths}) {
		my $length = $$cumulative_contig_lengths{$contigs};
		print $ofh "abline(v=$length, col=\"black\", lwd=1.5, lty=1)\n";
	}
}

print $ofh "dev.off()\n";
print $ofh "q()\n\n";
close $ofh;

# Run R script
my $CMD = "cat $R_file | R --vanilla";
system($CMD);	

sub count_different_relatives_and_assign_numbers {
	my $file = $_[0];
	warn "count_different_relatives: $file...\n";
	my %haplotype_keys;
	open my $fh, '<', $file or die "Cannot open $file : $!";
	while(my $line=<$fh>) {
		chomp $line;
		next if($line =~ m/Contig\tWindow_Start/);
		my @bits = split /\t/, $line;
		my ($contig, $start, $stop, $running_pos, $hap1, $hap2) = @bits;
		next if(($hap1 eq 0) || ($hap2 eq 0));
		$haplotype_keys{$hap1} = 1;
		$haplotype_keys{$hap2} = 1;
	}
	close $fh;

	# Give each relative a different number
	my %haplotype_keys2;
	my $haplotype_keys_string;
	my $count = 0;
	foreach my $hap_relatives(sort keys %haplotype_keys) {
		$haplotype_keys2{$hap_relatives} = $count;
		$count++;
		$haplotype_keys_string .= "\"$hap_relatives\",";
	}
	$haplotype_keys_string =~ s/,$//;
	return (\%haplotype_keys2, $haplotype_keys_string);
}

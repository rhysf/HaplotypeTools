
<img src="https://github.com/rhysf/HaplotypeTools/blob/main/resources/logo.png?raw=true" width="400" height="400" />


## Introduction

HaplotypeTools is a set of scripts designed to phase aligned WGS data, whereby
haplotypes are constructed from sufficient numbers of reads overlapping two or
more bi-allelic heterozygous position. Both haplotypes will need to be recovered
according to a cutoff (-c) in HaplotypeTools.pl. HaplotypeTools.pl can take a 
lond time when not submitted in parallel to a cluster, so where available it should
be run in parallel (via settings -g, -a and -q). Once phased, the VCF can be compared 
to consensus genomes belonging to relatives in order to identify possible parental
strains, or can be compared to other phased VCF's (or if multi VCF, samples within)
to look for crossovers.

Genomic regions that are phased in >=1 VCF (or multiple samples within a VCF) 
are called 'phased in any/all; PIA'. Phased in any (-p 1) is the default, while 
phased in all (-p 2) will output only regions phased in every sample.
PIA in either setting serve as an intermediate file needed to 
extract haplotypes in FASTA format, as well as HaplotypePlacer. For comparing
haplotypes to consensus genomes, phased in any (-p 1) is best, as it will maximise
the number of genomic regions, or phased in just a single sample for fewer sites
that are most informative to that sample. Some applications may require identifying
regions that are phased in 2+ samples, in which case phased in all (-p 2) should
be used.

Haplotype_placer.pl will compare haplotypes from a phased VCF and consensus
genomes for other isolates, and iteratively construct trees using FastTree to
identify the nearest relative of each haplotype, which can be useful for identifying
parents of recombinants or hybrids, as well as visualising breakpoints. The output
includes both summaries and window files that can be plotted using the script XXX.

VCF_phased_compare_to_VCF_phased.pl can be used to compare two phased VCFs and
directly identify genomic regions that contain crossovers between them.


## Documentation

All documentation for HaplotypeTools can be found at https://github.com/rhysf/HaplotypeTools

## Support

For issues, questions, comments or feature requests, please check or post to the issues tab on github: https://github.com/rhysf/HaplotypeTools/issues

## Prerequisites:

* R (plyr, RColorBrewer)
* Perl
* BioPerl
* Bio::DB::HTS
* Hash::Merge
* Samtools v0.1.10 or higher (samtools.sourceforge.net)
* FastTree (http://www.microbesonline.org/fasttree/)

## Getting started / example pipeline for phasing individual sample

```bash
git clone git@github.com:rhysf/HaplotypeTools.git
cd HaplotypeTools/
perl HaplotypeTools.pl -v <vcf> -b <sorted BAM> -u <VCF sample name> -f <reference.fasta>
perl util/VCF_phased_to_PIA.pl \
	-v <vcf-Phased-m-4-c-90-r-10000.vcf> \
	-f <reference.fasta>
perl util/VCF_phased_and_PIA_to_FASTA.pl \
	-v <vcf-Phased-m-4-c-90-r-10000.vcf> \
	-l <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab> \
	-r <reference.fasta>
perl util/FASTA_compare_sequences.pl \
	-f <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_1.fasta> \
	-a <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_2.fasta> \ 
	-o a > summary
```

## Example pipeline for phasing multi sample VCF

```bash
perl HaplotypeTools.pl -v <multi sample vcf> \
	-b <sorted BAMs (separated by comma)> \
	-u <VCF sample names in order of input BAM files (separated by comma)> \
	-f <reference.fasta>
perl util/VCF_phased_to_PIA.pl \
	-v <vcf-Phased-m-4-c-90-r-10000.vcf> \
	-f <reference.fasta>
	[if PIA is only for subset of samples] -s <VCF sample names to restrict analysis to>
perl util/VCF_phased_and_PIA_to_FASTA.pl \
	-v <vcf-Phased-m-4-c-90-r-10000.vcf> \
	-l <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab> \
	-r <reference.fasta>
	-u <sample name in VCF to pull haplotypes from>
perl util/FASTA_compare_sequences.pl \
	-f <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-<opt_u>.tab_1.fasta> \
	-a <vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-<opt_u>.tab_2.fasta> \ 
	-o a > summary
```

## Example pipeline for comparing haplotypes to other genomes

```bash
perl util/VCF_and_FASTA_to_consensus_FASTA.pl \
	-v <sample1.vcf> \
	-r <reference.fasta>
```

Will make your consensus genomes. Next, we need to make a tabular file in the 
following format (Name_Type_Location.tab) E.g.:

Sample1	FASTA	/dir/with/your/consensus/FASTA/sample1.vcf-WGS-s-n-i-n-n-N-consensus.fasta
Sample2	FASTA	/dir/with/your/consensus/FASTA/sample2.vcf-WGS-s-n-i-n-n-N-consensus.fasta
Sample3	FASTA	/dir/with/your/consensus/FASTA/sample3.vcf-WGS-s-n-i-n-n-N-consensus.fasta
...

And an optional tabular file for clades/lineages E.g.:
Sample1	clade1
Sample2	clade1
Sample3	clade2

```bash
perl util/VCF_phased_to_PIA.pl \
	-v <phased_sample1.vcf>,<phased_sample2.vcf>,<phased_sample3.vcf> \
	-f <reference.fasta>
perl util/VCF_phased_and_PIA_to_FASTA.pl \
	-v <phased_sample1.vcf> \
	-l <phased_sample1.vcf-plus_other_VCFs-PIA-p-1-c-10-s-all.tab> \
	-r <reference.fasta> 
perl util/Haplotype_placer.pl \
	-p <phased_sample1-plus_other_VCFs-PIA-p-1-c-10-s-all.tab> \
	-a <phased_sample1.vcf-plus_other_VCFs-PIA-p-1-c-10-s-all.tab_1.fasta> \
	-b <phased_sample1.vcf-plus_other_VCFs-PIA-p-1-c-10-s-all.tab_2.fasta> \
	-n <Name_Type_Location.tab>
```

## For window plots as well

```bash
perl util/Haplotype_placer.pl \
	-p <phased_sample1-plus_other_VCFs-PIA-p-1-c-10-s-all.tab> \
	-a <phased_sample1.vcf-plus_other_VCFs-PIA-p-1-c-10-s-all.tab_1.fasta> \
	-b <phased_sample1.vcf-plus_other_VCFs-PIA-p-1-c-10-s-all.tab_2.fasta> \
	-n <Name_Type_Location.tab> \
	-d y \
	-l <reference.fasta>
perl util/Windows_haplotypes_to_R_figure.pl -w HaplotypeTools_windows
```

##Example pipeline to identify crossovers

```bash
perl util/VCF_phased_compare_to_VCF_phased.pl \
	-a <phased_sample1.vcf> \
	-b <phased_sample2.vcf> \
	-o sample1_vs_sample2
```

or

```bash
perl util/VCF_phased_compare_to_VCF_phased.pl \
	-a <phased.vcf> \
	-c <sample ID 1 from VCF>
	-d <sample ID 2 from VCF>
	-o sample1_vs_sample2
```

## Example pipeline using test data (chr1 Bd) to phase individual sample and plot haploytpe windows

```bash
cd HaplotypeTools/example
perl ../HaplotypeTools.pl \
	-v Hybrid-SA.vcf-chr1.vcf \
	-b Hybrid-SA.vcf-chr1.bam \
	-u Hybrid-SA-EC3 \
	-f Hybrid-SA.vcf-chr1.fasta
```

* [note] if you can submit jobs to a cluster, use options -g, -a and -q to speed this step up.

```bash
perl ../util/VCF_phased_to_PIA.pl \
	-v Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf \
	-f Hybrid-SA.vcf-chr1.fasta

perl ../util/VCF_phased_and_PIA_to_FASTA.pl \
	-v Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf \
	-l Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab \
	-r Hybrid-SA.vcf-chr1.fasta

perl ../util/FASTA_compare_sequences.pl \
	-f Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_1.fasta \
	-a Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_2.fasta \
	-o a > summary

perl ../util/Haplotype_placer.pl \
	-p Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab \
	-a Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_1.fasta \
	-b Hybrid-SA.vcf-chr1.vcf-Phased-m-4-c-90-r-10000.vcf-PIA-p-1-c-10-s-all.tab_2.fasta \
	-n Name_Type_Location_consensus_genomes.tab \
	-d y \
	-l Hybrid-SA.vcf-chr1.fasta

perl util/Windows_haplotypes_to_R_figure.pl \
	-w HaplotypeTools_windows -s 7000
```

* [note] This should produce windows_plot-for-Haplotypes.R.pdf. For reference, an example copy is saved as Example-windows_plot-for-Haplotypes.R.pdf


## Individual script details:

* HaplotypeTools parameters are shown below, followed by their default settings given in [].
* Files are highlighted in <>.

```bash
HaplotypeTools.pl 
Parameters: -v <VCF> 
            -b <BAM (sorted) separated by comma if phasing multiple samples in a multi VCF> 
            -f <reference fasta>
	    -u VCF sample names in order of input BAM files (separated by comma if phasing multiple samples)
Optional:   -c Cut-off percent reads supporting phase group [90]
            -m Minimum read depth overlapping two heterozygous SNPs [4]
	    -r Max phase length [10000]
	    -s Steps (1=process VCF, 2=process BAM, 3=assign read info to VCF, 4=validate and assign phase groups, 5=concatenate) [12345]
Parallel:   -g Run commands on the grid (y/n) [n]
            -a Platform (UGER, LSF, GridEngine) [UGER]
	    -q Queue name [short]
Outputs:    -o Output folder for tmp files [opt_v-HaplotypeTools-phased-r-opt_r]
	    -p Phased VCF [opt_v-Phased-m-opt_m-c-opt_c-r-opt_r.vcf]
            -y Phased summary [opt_v-Phased-m-opt_m-c-opt_c-r-opt_r-summary.tab]

FASTA_compare_sequences.pl
Parameters: -f <FASTA>
Optional:   -o Output (p=pairwise within opt_f, s=summary of pairwise within opt_f, a=additional FASTA: 1:1 order) [p]
            -g Exclude if gaps percent is greater than this [100]
	    -a <Additional FASTA in same order> []

Haplotype_list_count.pl <haplotypes list> <full sequence file>

Haplotype_placer.pl 
Parameters: -p <phased in any/all; PIA>
            -a <Haplotype file 1 FASTA>
            -b <Haplotype file 2 FASTA>
            -n <Name tab Type tab Location for consensus genomes>
Optional:   -f Folder for output FASTA and Trees [HaplotypeTools_output]
            -c Clades/Lineages/Metadata to cluster isolates by. Format (Name tab Clades/Lineages/Metadata) []
            -t Location of FastTree [/seq/annotation/bio_tools/FastTree/current/FastTreeMP]
            -m Mimimum length of haplotypes [100]
            -v Verbose (y/n) [n]
Windows:    -d Calculate windows (y/n) [n]
            -l <Reference FASTA> []
            -z Window length [10000]
Output:     -f Output folder for genomic FASTA files and Trees [HaplotypeTools_output]
            -s Output summary [HaplotypeTools_summary]
	    -w Output windows [HaplotypeTools_windows]

Haplotype_FASTA_files_to_compare_to_IRMS_het_sites.pl
Parameters: -a <PIA_1 FASTA> 
            -b <PIA_2 FASTA> 
            -c <IRMS HET.fasta> 
            -d <IRMS HET.details>
Optional:   -p print to stdout (s) or (f) [s]
            -o File to print to if opt_p=f <opt_a-and-2.accuracy>
            -z Verbose for testing (y/n) [n]


VCF_and_FASTA_to_consensus_FASTA.pl 
Parameters: -v <VCF> 
            -r <reference FASTA>
Optional:   -i Include homozygous indels (y/n) [n]
            -h Include bi-allelic heterozygous positions if ref base not included (y/n) [y]
            -a For heterozygous positions, use ambiguity code (a) or first position for haploid consensus (f) [f]
	    -n Character for ambiguous [N]
	    -s Restrict to only this supercontig [n]
	    -t Restrict to only this isolate [n]
Notes: Prints to opt_v-isolate-name-consensus.fasta

VCF_phase.pl
Parameters: -v <VCF from step 1 of HaplotypeTools>
            -a <alignment file from step 2 of HaplotypeTools>

Optional:   -s Sample number
Output:     -o Output [opt_v-phased-opt_s]

VCF_phase_haploid.pl <VCF file 1> > outfile.vcf

VCF_phased_and_PIA_to_FASTA.pl 
Parameters: -v <VCF file> 
            -l <PIA file (contig tab start tab stop)> 
            -r <Reference FASTA>
Optional:   -u Sample name in VCF to pull haplotypes from [WGS]
            -p Printing option (o=outfile, s=split to opt_l_1 and _2) [s]
            -e Exclude printing if haplotypes are identical (y/n) [n]
	    -m min length [10]
	    -i Include indels (y/n) [n]
	    -z Verbose for testing (y/n) [n]

VCF_phased_append_name_to_phase_group.pl <VCF file 1> > outfile.vcf

VCF_phased_calculate_haplotype_lengths.pl <Phased VCF file (with optional additional VCFs separated by comma)>
Optional:  -s Sample name (with optional additional sample names separated by comma)
           -e Output results for every sample (y/n) [n]

VCF_phased_compare_to_VCF_phased.pl
Parameters: -a <VCF file 1>
            -b <VCF file 2>
	    or
	    -a <VCF>
	    -c Sample ID 1
	    -d Sample ID 2
Optional:   -o Outfile extension [HaplotypeTools_VCF1_vs_VCF2]

VCF_phased_to_PIA.pl
Parameters: -v <VCF's (separated by comma)> 
            -f <reference FASTA>
Optional:   -c Length cut-off for minimum haplotype to be considered [10]
            -p Phased in any (1), Phased in all (2) [1]
            -s Sample names to restrict analysis too (separated by comma)
            -t Phase Tag (HaplotypeTools uses PID, Whatshap uses PS or HP etc) [PID]
Output:     -u Summary [opt-v-PIA-p-Opt_p-c-Opt_c-s-Opt_s.summary]
            -o Output [opt-v-PIA-p-Opt_p-c-Opt_c-s-Opt_s.tab]

VCF_phased_validate_and_assign_phase_groups.pl
Parameters: -v <VCF-phased-Sample_Number (from VCF_phase.pl)>
Optional:   -c Cut-off percent reads supporting phase group [90]
            -m Minimum read depth overlapping two heterozygous SNPs [4]
            -z Verbose for error checking (y/n) [n]
Output:     -o Output VCF [opt_v-and-assigned.tab]
            -l Tallies for percent of reads agreeing with phases [opt_v-Tally-agree-disagree-opt_m-MD-opt_c-pc_cutoff.tab]

Windows_haplotypes_to_R_figure.pl
Parameters: -w <Windows haplotypes dataframe>
Optional    -s Scaling factor for width (nucleotides per mm) [20000]
            -h Scaling factor for height (genome per mm) [70]
            -r Resolution (The nominal resolution in ppi) [100]
            -c Minimum contig length to include label [1000000]
            -g Optional size of axis text [2]
            -l Legend size [2]
            -x X-min [0]
```
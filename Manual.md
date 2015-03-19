

# Requirements #

  * Perl 5.8 or higher
  * R 2.7 or higher
  * a short read aligner (e.g.: bwa)
  * samtools
  * Perl libraries:
    * Text::NSP::Measures::2D::Fisher::twotailed

# Experimental Design #

PoPoolation2 is designed for comparing allele frequencies between two or more populations and should thus solely be used with pooled genomic DNA. Most importantly the amount of DNA per individual in a single pool should be constant. Furthermore the number of individuals per pool should be similar. Sequencing data should be available for two or more populations. At least a partial assembly of the sequenced species should be available.

We furthermore recommend that the pool size (number of individuals in the pool) should be larger than the coverage. This minimizes re-sampling the same allele from a single individual several times.

The accuracy of the estimated differences in allele frequencies between populations will critically depend on the pool size and on the sequencing coverage. In principle it should thus be aimed to maximize both the pool size and the coverage.
However ultimately the required coverage and pool size will depend on the effect size that needs to be reliably detected. For example to detect significant (Fisher's exact test; p=0.01) allele frequency differences between two populations of 30%, a coverage (and pool size) of 50 will be sufficient. If however significant (Fisher's exact test; p=0.01) allele frequency differences of 10% need to be detected the coverage as well as the pool size should be about 400.


# Download PoPoolation2 #

PoPoolation2 may be obtained directly from the subversion repository. Just enter the command line and move to the directory where you want to install PoPoolation2 and type the following command:
```
svn checkout http://popoolation2.googlecode.com/svn/trunk/ popoolation2
```

To update your copy of PoPoolation2 with the latest improvements enter your PoPoolation2 directory and enter the command:
```
svn update
```

Alternatively PoPoolation2 may be downloaded from the project main page: http://code.google.com/p/popoolation2/

However we recommend to use subversion as bugfixes will be immediately available in the repository.

# Benchmarks for main scripts #
This section provides some estimates about computation time and RAM consumption.

| script | time (min:sec) | RAM (MB)|
|:-------|:---------------|:--------|
| mpileup2sync.pl |11:56 | 5 |
| mpileup2sync.pl |00:16 | 317 |
| snp-frequency-diff.pl | 2:07 | 143 |
| fst-sliding.pl | 1:09 | 174 |
| fisher-test.pl | 3:00 | 174 |
| cmh-test.pl | 1:36 |201|


`mpileup2sync.pl` has been tested with a mpileup file for 2 populations having 1 million entries and a coverage of 100 per population. `snp-frequency-diff.pl` has been tested with a synchronized file (see below) for 2 populations having 1 million entries. `fst-sliding.pl` has been tested with a synchronized file for 2 populations having 1 million entries. We used a window size of 1000 bp, the memory consumption will increase with the window size. For example, when increasing the window size to 1,000,000 bp the RAM demand will increase to 1357 MB. `fisher-test.pl` has been tested with a synchronized file for 2 populations having 1 million entries. We used a window size of 1 (calculating a p-value for every SNP). Memory demand will increase with the window size as described for `fst-sliding.pl`. `cmh-test.pl` was tested using a synchronized file for 4 populations having 1 million entries.

All benchmarks have been done on a Mac OS X 10.5.8, 2x2.26 GHz Quad-Core Intel Xeon, with 16GB of RAM using a single thread (CPU).


# PoPoolation2 scripts #

## General note ##

Help for all PoPoolation2 scripts can be requested using the '--help' option.
For example
```
perl fst-sliding.pl --help
```

All main scripts of PoPoolation2 contain Unit tests which can be run using the '--test' option
For example:
```
perl trim-fastq.pl --test
```
We generally recommend to run this tests before using a script.

## mpileup2sync.pl ##

Converts a mpileup file into a so called 'synchronized' format. mpileup files can be obtained with samtools and the 'synchronized' file format is the standard format used with all PoPoolation2 scripts. The synchronized file contain the allele frequencies for all bases in the reference genome and for all populations being analyzed.
Following an example of two populations:

```
2L	79	G	0:0:0:15:0:0	0:0:0:38:0:0
2L	80	A	12:0:0:0:0:0	38:0:0:0:0:0
2L	81	A	14:0:0:0:0:0	43:0:0:0:0:0
2L	82	A	14:0:0:0:0:0	42:0:0:0:0:0
```

  * column 1: reference contig
  * column 2: position in the reference contig
  * column 3: refernce character
  * column >3: allele frequencies for all populations in the form A-count:T-count:C-count:G-count:N-count:deletion-count

**Note** the synchronized file may contain entries for two ore more populations
**Note** allele frequencies are provided after discarding low quality bases (minimum quality)

## mpileup2sync.jar ##

The same as above only implemented in Java-multithreading. Synchronizing the mpileup file is a major bottleneck in data analysis with PoPoolation2. To speed up data processing we implemented this tool in Java multi-threading which speeds up the analysis about 78x (with 8 CPUs).
The input, output, and commands are identical to the implementation in Perl.

## synchronize-pileup.pl ##

Allows to create a synchronized file (see above) when starting with several pileup files. Before the mpileup file format was introduced samtools only supported the pileup file format. We thus only recommend to use this script when you have old legacy pileup files and no mpileup file can be obtained!

This script was necessary as the (outdated) pileup files may start at any position with respect to the reference contig (eg. file 1 starts at position 12 on chromosome 2L whereas file 2 starts at position 100 of chromosome 2L). It was thus necessary to put these files 'in sync' where entries refering to the same positions (eg.: 2L:20,001) are grouped into the same column.

## snp-frequency-diff.pl ##

Provides detailed statistics about the major and minor alleles for all SNPs in the given populations. Furthermore calculates for every SNP, the exact allele frequency differences for all pairwise comparisons of the given populations.

## fst-test.pl ##
Allows to calculate the pairwise Fst using a sliding window approach. All possible pairwise comparisons will be calculated. Fst values may either be calculated in the classical way (http://www.sinauer.com/detail.php?id=3082)  or by an equation suggested by Karlsson et al. (http://www.nature.com/ng/journal/v39/n11/full/ng.2007.10.html). Furthermore Fst values may be calculated for windows or single bases. Following an example output
```
2L      68500   360     1.000   62.1    1:2=0.01873725  1:3=0.02131245  2:3=0.01521177
2L      69500   118     1.000   71.9    1:2=0.00969479  1:3=0.00116059  2:3=0.00905794
2L      70500   269     1.000   63.6    1:2=0.01955417  1:3=0.01547995  2:3=0.01300569
```

  * column 1: reference contig
  * column 2: position of the window (middle)
  * column 3: SNPs identified in the window
  * column 4: fraction of window having sufficient coverage
  * column 5: average minimum coverage
  * column >5: Fst values for all pairwise comparisons; For example "2:3=0.02" states that the Fst for comparing population 2 with population 3 (using the order of the synchronized file) is 0.2

## fisher-test.pl ##
Allows to calculate the pairwise Fisher's exact test for all populations being present in a given synchronized file. The Fisher's exact test may be used to estimate the significance of allele frequency differences between populations. The output format is identical to `fst-sliding.pl`. P-values are provided in -log10(p-value)

## cmh-test.pl ##
Allows to perform a Cochran-Mantel-Haenszel test for repeated tests of independence for every SNP. This may be useful for experiments where several biological replicates are used. This test thus identifies consistent changes in allele frequencies for several replicates.

To elucidate the use of the cmh-test following a rather contrived example in Arabidopsis: Imagine you want to identify genes responsible for height. In this case you may separate all small and all large plants into groups and sequence these groups independently. Now you could measure the allele frequency differences between these groups. However to increase the resolution and reduce the number of false positives you may repeat this several times, maybe using Arabidopsis populations from very distinctive geographical locations. In this case the cmh-test may help to identify consistent changes in allele frequencies between small and large plants.

**NOTE**: with the feature population it is essential that every population is only tested once. ie: 1-4,2-5,3-6 is correct whereas 1-2,1-3,1-4 is false and violates the assumptions of the cmh-test

**NOTE**: The cmh test will check the min-count, min-coverage, max-coverage of ALL populations in the sync-file, not just for those being compared in the test.  If these restrictions should however only be applied to the tested populations, a sync-file containing only the relevant populations may be created with awk (cat inputfile.sync | awk 'BEGIN{OFS="\t"}{print $1,$2,$....} > outputfile.sync')

## subsample-synchronized.pl ##

This script reduces the coverage in a synchronized file for every population to the targeted coverage by random sampling of bases.

As genomic reads are randomly distributed in the genome the coverage of mapped reads shows marked fluctuations along chromosomes. There are also well known biases like the GC bias, in this case regions having a high GC content also have elevated coverages.
Statistical test's, like the Fisher's exact test or the CMH-test, more readily identify allele frequency differences between populations in regions having high coverages. This may result in artefactual results as for example a higher density of significant allele frequency differences in regions having a high GC content. This script allows to subsample bases to a uniform coverage, which should thus eliminate artefactual results that are caused by coverage fluctuations.

However, several methods for subsampling have been implemented (with replacement, without replacement, exact fraction).

## create-genewise-sync.pl ##
This script is a central inovation of PoPoolation2 and enables the analysis of genes, exons, introns or any other feature present in a gtf file. It converts a synchronized file which has a coordinate system that is typically based on reference chromosomes into a synchronized file with a coordinate system based on the features present in a given `gtf-file` (http://mblab.wustl.edu/GTF2.html). This conversion allows to use all analysis of PoPoolation2 (Fst, fisher's exact test, CMH-test, allele frequency differences) for the given features.
Example before conversion:
```
2L	5088	T	0:5:0:0:0:0	0:6:0:0:0:0
2L	5089	T	0:4:0:0:0:0	0:5:0:0:0:0
...
2R	1001	T	0:5:0:0:0:0	0:6:0:0:0:0
2R	1002	T	0:4:0:0:0:0	0:5:0:0:0:0
...
```

Example after conversion:
```
CG40127-RA      1       T       0:95:0:0:0:0    0:95:0:0:0:0
CG40127-RA      2       G       0:0:0:97:0:0    0:0:0:97:0:0
CG40127-RA      3       A       98:0:0:0:0:0    98:0:0:0:0:0
```

Internally all features having the same `gene_id` are grouped and this `gene_id` will act as the new id of the reference contigs. This grouping of features having the same `gene_id` allows to deal with exon intron structure. Note that the `transcript_id` of a gtf file is ignored. The position of the entries will be provided with respect to this `gene_id`, where intronic sequences will be ignored. Features on the minus strand will be reverse complemented (including the reference character and the allele-counts).
The script also deals with overlapping features where two cases need to be distinguished. If the overlapping features have the same `gene_id` the script calculates the union (e.g.: many exons for every gene). On the other hand, if they have a different `gene_id` they are treated as separate entries and output will be provided for every individual `gene_id`. After conversion the main analysis of PoPoolation2 may be performed either with a sliding window or for individuals SNPs. If for example Fst should be calculated for the whole genes than the window size should be set to a number being longer than the length of any gene (eg.: 1000000).


## export/pwc2igv.pl ##
Converts the results of scripts conducting pairwise comparisions (fiser-test.pl fst-sliding.pl) into the `.igv` format. This may be directly loaded into the IGV (http://www.broadinstitute.org/software/igv/home). Note that the `.igv` file format allows to load multiple tracks simultanously.

## export/cmh2gwas.pl ##
Converts the results of `cmh-test.pl` into the `.gwas` format, which may be directly imported into the IGV (http://www.broadinstitute.org/software/igv/home). This may be useful to visualize results of genome wide association studies

## export/subsample\_sync2fasta.pl ##
Converts a synchronized file into a multiple fasta file. This may be useful to use third party software for Population Genetics like Mega5 (DnaSP). For every population present in the synchronized file the same number (`--target-coverage`) of fasta entries  will be created.
As conversion into fasta file format requires a uniform coverage a random sub-sampling step is performed internally.

**Warning** Pooling of individuals results in a loss of haplotype information, the haplotype information present in the fasta file is thus merely an artefact. Thus, do not conduct any analysis with the resulting fasta files that use haplotype information.

## export/subsmaple\_sync2GenePop.pl ##
Converts a synchronized file into a GenePop file. This file may be used with GenePop (http://genepop.curtin.edu.au/) and Arlequin (http://cmpg.unibe.ch/software/arlequin35/). GenePop files may be directly imported with Arlequin or converted into the Arlequin file format (.arp) at the GenePop webpage.
As conversion into GenePop files requires a uniform coverage a random sub-sampling step is performed internally.
Following an example of the GenePop file:
```
Microsatelite loci on Chiracus radioactivus, a pest species
Loc1, Loc2, Loc3
Pop
AA8, 04 07 03
AA9, 04 06 02
A10, 02 06 01
A11, 04 06 01
Pop
AF, 00 00 00
AF, 02 03 01
AF, 02 03 02
AF, 02 03 00
```

**Warning** Pooling of individuals results in a loss of haplotype information, the haplotype information present in the fasta file is thus merely an artefact. Thus, do not conduct any analysis with the resulting fasta files that use haplotype information.

## indel\_filtering/identify-indel-regions.pl ##
Identifies regions containing an indel from a mpileup file and converts the location of regions surrounding the indels into a `.gtf` file format (http://mblab.wustl.edu/GTF2.html)
May be used to remove regions around indels from synchronized files.

## indel\_filtering/filter-sync-by-gtf.pl ##

Remove the regions specified in a `.gtf` file (see above) from a synchronized file. Maybe used to filter regions around indels.

# After PoPoolation2 #
What to do after you successfully used PoPoolation2,  having identified several hundreds (or thousands) SNPs being associated with a trait of interest (or environmental difference). You may use Gowinda to test whether your SNPs show an enrichment for any GO category (or any gene set): http://code.google.com/p/gowinda/
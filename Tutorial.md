

# Introduction #

PoPoolation2 facilitates comparison of allele frequencies between different populations from pooled next generation sequencing data. In this walkthrough the main functionality of PoPoolation2 will be demonstrated by comparing the allele frequencies of two populations.

The sample data set used in this walkthrough consists of simulated `fastq` reads derived from the first 1.4 Mbp of chromosome 2R of _D.melanogaster_. Two pooled populations have been simulated each having a coverage of 100. In total 10.000 SNPs, located in unambiguous regions, have been introduced.  We simulated different allele frequencies for these SNPs in the two populations, with differences ranging from 0.1 to 0.9. We also used a minimum distance of 100 bp between two successive SNPs. The simulated reads have a length of 75 bp and a sequencing error rate of 1%.

# Requirements #

  * PoPoolation2
  * Perl http://www.perl.org/
  * R http://www.r-project.org/
  * bwa http://bio-bwa.sourceforge.net/
  * samtools http://samtools.sourceforge.net/
  * IGV http://www.broadinstitute.org/igv/


# Data #

The sample data can be downloaded here: http://popoolation2.googlecode.com/files/demo-data.zip

Subsequently unzip the archive.


# Walkthrough #

## Prepare the reference genome ##
Open the command line and change directory to the unziped folder containing the sample data of PoPoolation2
```
mkdir ref
mv 2R.chr ref
bwa index ref/2R.chr
```

## Map the reads to the reference genome ##
```
mkdir map
bwa aln -n 0.01 -l 100 -o 1 -d 12 -e 12 -t 8 ref/2R.chr pop1.fastq > map/pop1.sai
bwa aln -n 0.01 -l 100 -o 1 -d 12 -e 12 -t 8 ref/2R.chr pop2.fastq > map/pop2.sai
bwa samse ref/2R.chr map/pop1.sai pop1.fastq > map/pop1.sam
bwa samse ref/2R.chr map/pop2.sai pop2.fastq > map/pop2.sam
```

## Remove ambiguously mapped reads ##
```
samtools view -q 20 -bS map/pop1.sam | samtools sort - map/pop1
samtools view -q 20 -bS map/pop2.sam | samtools sort - map/pop2
```

## Create a synchronized file ##

Synchronized files are the main input files for PoPoolation2. They basically contain the allele frequencies for every population at every base in the reference genome in a concise format. Note that the synchronized file format contains the allele frequencies after filtering for base quality.

```
samtools mpileup -B map/pop1.bam map/pop2.bam > p1_p2.mpileup
perl <popoolation2-path>/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input p1_p2.mpileup --output p1_p2.sync
```

Synchronizing the mpileup file is quite time consuming. To remove this bottleneck we implemented 'mpileup2sync' in Java multi-threading which is about 78x faster as the implementation in perl.
```
java -ea -Xmx7g -jar <popoolation2-path>/mpileup2sync.jar --input p1_p2.mpileup --output p1_p2_java.sync --fastq-type sanger --min-qual 20 --threads 8
```


Sample of a synchronized file:

```
2R	2302	N	0:7:0:0:0:0	0:7:0:0:0:0
2R	2303	N	0:8:0:0:0:0	0:8:0:0:0:0
2R	2304	N	0:0:9:0:0:0	0:0:9:0:0:0
2R	2305	N	1:0:9:0:0:0	0:0:9:1:0:0
```

  * col1: reference contig
  * col2: position within the refernce contig
  * col3: reference character
  * col4: allele frequencies of population number 1
  * col5: allele frequencies of population number 2
  * coln: allele frequencies of population number n

The allele frequencies are in the format A:T:C:G:N:del, i.e: count of bases 'A', count of bases 'T',... and deletion count in the end (character '`*`' in the mpileup)

## Calculate allele frequency differences ##

The exact allele frequency differences may be computed as follows:
```
perl <popoolation2-path>/snp-frequency-diff.pl --input p1_p2.sync --output-prefix p1_p2 --min-count 6 --min-coverage 50 --max-coverage 200
```

This script creates two output files having two different extensions:
  * `_`rc: this file contains the major and minor alleles for every SNP in a concise format
  * `_`pwc: this file contains the differences in allele frequencies for  every pairwise comparision of the populations present in the synchronized file
For details see the man pages of the script

The allele frequency differences can be found in the **`_`pwc** file, a small sample:
```
##chr   pos     rc      allele_count    allele_states   deletion_sum    snp_type        most_variable_allele    diff:1-2
2R      4459    N       2       C/T     0       pop     T       0.133
2R      9728    N       2       T/C     0       pop     T       0.116
```
The last column contains the obtained differences in allele frequencies for the allele provided in column _8_. Note that in this example the last column refers to a pairwise comparision between population _1 vs 2_, in case several populations are provided all pairwise comparisions will be appended in additional columns.
For a comparison of the observed allele frequency differences with the expected ones please see: [Validation](Validation.md)



## Fst-values: measure differentiation between populations ##

### Calculate Fst for every SNP ###

```
perl <popoolation2-path>/fst-sliding.pl --input p1_p2.sync --output p1_p2.fst --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 500
```

the option `--suppress-noninformative` suppresses output for windows not containing a SNP. When applied to windows of size one, this option suppresses output for bases that are no SNP.

For a validation of the Fst values please see [Validation](Validation.md)

### Calculate Fst values using a sliding window approach ###

```
perl <popoolation2-path>/fst-sliding.pl --input p1_p2.sync --output p1_p2_w500.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 500 --step-size 500 --pool-size 500
```

### Load Fst-values into IGV ###
First index the two bam files
```
samtools index map/pop1.bam
samtools index map/pop2.bam
```
Than convert the file containing the Fst-values into the `.igv` format
```
perl <popoolation2-path>/export/pwc2igv.pl --input p1_p2.fst --output p1_p2.igv
```

Finally load the results into the IGV

  * Start IGV
  * Import Genome -> Load the fasta file:  2R.fasta
  * Load the two `.bam` files
  * Load the `.igv` file containing the Fst values

Voila, the result should look something like this:
![http://popoolation2.googlecode.com/files/igv_fst_single.png](http://popoolation2.googlecode.com/files/igv_fst_single.png)

Note that the result shown above only contain the Fst for a single pairwise comparison (pop1 vs pop2). PoPoolation2 is per default calculating the Fst for all possible pairwise comparisons between populations (when several populations are present in the synchronized file), and these results may just as easily be converted into the `.igv` format. When loaded into IGV it may look like the following:

![http://popoolation2.googlecode.com/files/igv_fst_multi.png](http://popoolation2.googlecode.com/files/igv_fst_multi.png)


## Fisher's Exact Test: estimate the significance of allele frequency differences ##

The Fishers exact test can be used to test whether any differences in allele frequencies are statistically significant. At low coverages the absolute changes of allele frequencies or the Fst values may be strongly influenced by sampling effects, therefore the Fishers exact test may be used to identify significant changes in allele frequency.

```
perl <popoolation2-path>/fisher-test.pl --input p1_p2.sync --output p1_p2.fet --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative
```

Note that this script provides -log10(P-value) as output.

### Load the Fisher's exact test results into the IGV ###
First convert the Fisher's exact test results into the `.igv` file format:
```
perl <popoolation2-path>/export/pwc2igv.pl --input p1_p2.fet --output p1_p2_fet.igv
```
Than start the IGV, load the `.bam` alignments and the file: `p1_p2_fet.igv`. The result should be very similar as in the following
example:
![http://popoolation2.googlecode.com/files/fet_example.png](http://popoolation2.googlecode.com/files/fet_example.png)

## Cochran-Mantel-Haenszel test: detect consistent allele frequency changes in several biological replicates ##
The cmh-test (http://stat.ethz.ch/R-manual/R-devel/library/stats/html/mantelhaen.test.html) is used for detecting significant and consistent changes in allele frequencies when independent measurements of the allele frequencies have been obtained (e.g.: biological replicates). For example when you are interested in the SNPs responsible for plant height, you may obtain several measurements of the allele frequencies of tall and small plants from different ecological regions. Using the cmh-test you may want to detect consistent differences (consistent between small and tall) in allele frequency for the following comparisions: `small_alaska vs tall_alaska + small_norway vs tall_norway + small_sibiria vs tall_sibiria`

In this walkthrough it is necessary to create a data set having a biological replicate. We are just duplicating the last two columns of the synchronized file used above:
```
mkdir cmh
cat p1_p2.sync|awk 'BEGIN{OFS="\t"}{print $0,$4,$5}' > cmh/p1_p2_p1_p2.sync
```

Here is an example of the result:
```
2R	2298	N	3:0:0:0:0:0	3:0:0:0:0:0	3:0:0:0:0:0	3:0:0:0:0:0
2R	2299	N	0:4:0:0:0:0	0:4:0:0:0:0	0:4:0:0:0:0	0:4:0:0:0:0
2R	2300	N	5:0:0:0:0:0	5:0:0:0:0:0	5:0:0:0:0:0	5:0:0:0:0:0

```

Compare the allele frequencies of popualtion 1 with population 2 and of population 3 with population 4 using the cmh-test.

```
perl <popoolation2-path>/cmh-test.pl --input p1_p2_p1_p2.sync --output p1_p2_p1_p2.cmh --min-count 12 --min-coverage 50 --max-coverage 200 --population 1-2,3-4
```
Note that it is necessary to double the `--min-count` as we also doubled the number of populations. For a validation of the cmh-test please see [Validation](Validation.md)

### Display the cmh-test results in the IGV ###

First convert the results of the cmh-test into the `.gwas` format of IGV.
```
perl <popoolation2-path>/export/cmh2gwas.pl --input p1_p2_p1_p2.cmh --output p1_p2_p1_p2.gwas --min-pvalue 1.0e-20
```
Note that IGV has some problems displaying very low p-values, thus the option `--min-pvalue` was used introduced as threshold. P-values smaller than the threshold will be set to the value of the threshold.

Open the IGV as described before and load the file `p1_p2_p1_p2.gwas`. The results should be similar to the following example

![http://popoolation2.googlecode.com/files/IGV_gwas_example.png](http://popoolation2.googlecode.com/files/IGV_gwas_example.png)

## Calculate Fst for genes ##

Download the annotated exons for the first 1mio bp of chromosome 2R http://popoolation2.googlecode.com/files/2R_exons.gtf

Convert the synchronized file into a gene-based synchronized file

```
perl <popoolation2-path>/create-genewise-sync.pl --input p1_p2.sync --gtf 2R_exons.gtf --output p1_p2_genes.sync
```

Calculate the Fst for every gene:

```
perl <popoolation2-path>/fst-sliding.pl --min-count 6 --min-coverage 50 --max-coverage 200 --pool-size 500 --min-covered-fraction 0.0 --window-size 1000000 --step-size 1000000 --input p1_p2_genes.sync --output p1_p2_genewise.fst
```

The output should be similar to the following example:
```
CG17683-RB	500000	17	0.002	100.0	1:2=0.01349929
CG40127-RA	500000	5	0.001	99.7	1:2=0.01897121
CG17665-RB	500000	27	0.003	99.9	1:2=0.03149757
```

The output above indicates that gene `CG17683-RB` has an Fst of `0.01349929`, gene `CG40127-RA` an Fst of `0.01897121` and so on.
Of course this analysis may also be conducted with more than two populations, in which case again all pairwise comparisons will be calculated.

Similarly the fisher's exact test may be performed for every gene, or the CMH-test for the SNPs of every gene.

**Note**: this analysis may be performed for any feature (intron, cds etc) as long as it is present in the gtf file and has a unique `gene_id`


# After PoPoolation2 #
What to do now? You have successfully identified several hundreds SNPs being associated with a trait of interest (or environment) but what is their biological meaning? You may use Gowinda to test whether your SNPs show an enrichment for any GO category (or any gene set): http://code.google.com/p/gowinda/
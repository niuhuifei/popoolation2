

# Introduction #

In the [Tutorial](Tutorial.md) a small sample data set consisting of simulated fastq-reads was introduced, containing simulated data reads of two different populations. We introduced 10.000 SNPs into the sample and simulated varying allele frequencies for these SNPs. SNPs were only introduced in unambigous regions of the genome and a sequencing error rate of 1% was used. As the expected allele frequencies, for these SNPs are known, we can compare them with the observed ones allowing to validate the functionality of PoPoolation2.



# Data and study design #

A proper validation requires that the expected values are obtained by different scripts than the observed values.
For this reason we calculated all observed values (e.g: CMH P-values, Fst-values, allele frequency difference) with PoPoolation2 and all expected values with separate/independent scripts (see 'validation-scripts').

Furthermore all expected values are directly calculated from the targeted allele frequencies http://popoolation2.googlecode.com/files/expected-snp-frequencies.sync, whereas the observed values are calculated from the results obtained after simulating reads and mapping them to the reference genome http://popoolation2.googlecode.com/files/observed.sync.zip

  * the 'validation-scripts' can be found here: http://popoolation2.googlecode.com/files/validation_scripts.zip
  * The expected allele frequencies for the 10.000 SNPs can be found here http://popoolation2.googlecode.com/files/expected-snp-frequencies.sync
  * The observed allele frequencies of all bases, including the 10.000 SNPs, can be found here http://popoolation2.googlecode.com/files/observed.sync.zip




# Number of SNPs #

First we tested whether all 10.000 simulated SNPs were recovered by PoPoolation2.We found 9.999 were recovered and a single SNP (at position 1,120,507) was missing. This SNP is missing due to a low coverage in one population (15), causing the SNP to be ignored during filtering (minimum coverage 50). We furthermore found that for this SNP, only one allelic state was present ('C') and the other allelic state ('T') was entirely missing. We can only speculate as to what is causing this problem, we however suspect ambiguous mapping of reads with the missing allelic state ('T').

We furthermore identified 211 SNPs that were not in the expected set of SNPs, thus this 211 SNPs are entirely due to sequencing errors and false alignments. However these 'artefactual' SNPs show very small allele frequency differences (mean: 0.02386) and are thus not likely to cause signals of differentiation between the two populations. See also the following histogram of the allele frequency difference of these 211 SNPs.


![http://popoolation2.googlecode.com/files/aretefactual_snps.png](http://popoolation2.googlecode.com/files/aretefactual_snps.png)

# Observed and expected allele frequency differences #
We found a strong correlation between the expected and the observed allele frequency differences (R^2=0.9979; P < 2.2e-16; 9,999 tested SNPs), demonstrating that PoPoolation2 highly reliably recovers allele frequency differences. See also the following graph for the detailed correlation between observed and expected allele frequency differences:

![http://popoolation2.googlecode.com/files/correlation.png](http://popoolation2.googlecode.com/files/correlation.png)


The 95% quantiles of the error in the estimated allele frequency differences are -2.0% and 2.4%. That is, 95% of the estimated allele frequency differences are within -2% to 2.4% of the real allele frequency difference. A detailed distribution of the error in the estimated allele frequency difference can be found in the following graph:

![http://popoolation2.googlecode.com/files/errordistri.png](http://popoolation2.googlecode.com/files/errordistri.png)

# Observed and expected Fst values #
We found a very strong correlation between the expected and the observed Fst values (R^2 = 0.9967; P < 2.2e-16). The observed Fst was calculated for every SNP as shown in the [Tutorial](Tutorial.md). For details see the following graph:

![http://popoolation2.googlecode.com/files/fst-correlation.png](http://popoolation2.googlecode.com/files/fst-correlation.png)

When computing the difference between observed and expected Fst we found a small bias, see graph:
![http://popoolation2.googlecode.com/files/fst-error.png](http://popoolation2.googlecode.com/files/fst-error.png)

This distribution is shifted towards the right side, which means that on average the observed Fst values are smaller than the expected ones, thus the Fst estimate is slightly underestimated.

# Observed and expected CMH-test p-values #
Observed CMH values were created as shown in the [Tutorial](Tutorial.md). We used the -log10 of the p-values and found a strong correlation between the observed and the expected p-values (Spearman's rank correlation: Rho= 0.9990084, P < 2.2e-16; Linear model: R^2=0.9978, P < 2.2e-16). For details see the following graph.

![http://popoolation2.googlecode.com/files/cmh_exp_obs.png](http://popoolation2.googlecode.com/files/cmh_exp_obs.png)

The following graph shows the error distribution of the CMH-values

![http://popoolation2.googlecode.com/files/cmh_exp_obs_hist.png](http://popoolation2.googlecode.com/files/cmh_exp_obs_hist.png)

Again a small bias can be found, in general the observed data show slightly elevated log transformed p-values as compared to the expected ones.

# Observed and expected Fisher exact test p-values #

We calculated the significance of allele frequency differences using a Fisher's exact test as described in [Tutorial](Tutorial.md). We found a strong correlation (Spearman's rank correlation: Rho=0.9989923, P < 2.2e-16; Linear model: R^2=0.9974, P < 2.2e-16) between the observed and expected p-values obtained with Fisher's exact test, for details see:

![http://popoolation2.googlecode.com/files/obs_expected_fet.png](http://popoolation2.googlecode.com/files/obs_expected_fet.png)

Note that we used the -log10(p-value) for the correlation.

For a distribution of the errors obtained with the Fisher's exact test see the following graph:

![http://popoolation2.googlecode.com/files/fet_errordist.png](http://popoolation2.googlecode.com/files/fet_errordist.png)

Also the p-values calculated with the Fisher's exact test show the small bias mentioned above, in general the observed data show slightly elevated log transformed p-values as compared to the expected ones.


# Conclusion #
in the paragraphs above, we tested the main functionality of PoPoolation2
  * measure allele frequency differences between populations
  * calculate pairwise Fst-values between populations
  * use the Fisher's exact test to estimate the significance of allele frequency differences between populations
  * compute the cmh-test for estimating the significance of allele frequency differences when having several biological replicates
As we used separate scripts for calculate the observed and the expected values, we conclude that PoPoolation2 highly accurately reproduces differences in allele frequencies between populations. However, small differences are still found between observed and expected values (see above), which may be caused by the simulated sequencing errors (1% error rate) or inaccuracies during mapping of the reads.
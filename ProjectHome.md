

# Introduction #
PoPoolation2 allows to compare allele frequencies for SNPs between two or more populations and to identify significant differences. PoPoolation2 requires next generation sequencing data of pooled genomic DNA (Pool-Seq). It may be used for measuring differentiation between populations, for genome wide association studies and for experimental evolution.

# News #
As of version 1.201 Mpileup2sync has been implemented in Java multithreading. This novel implementation is about 78x faster!

# Usage #

  * Tutorial: [Tutorial](Tutorial.md)
  * Manual: [Manual](Manual.md)
  * Bugs and Questions: [BugsAndQuestions](BugsAndQuestions.md)
  * Slides from a course: http://drrobertkofler.wikispaces.com/PoPoolationGenomics

# How to cite PoPoolation2 #

Please cite the following:
  * Kofler,R.,Vinay Pandey, R. & Schloetterer, C PoPoolation2: Identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq); Bioinformatics; Vol. 27 no. 24 2011, pages 3435–3436; http://bioinformatics.oxfordjournals.org/content/27/24/3435

You may also be interested in our Pool-seq review (Nature Reviews Genetics) where we provide some recommendations for the analysis of Pool-seq data:

  * http://www.nature.com/nrg/journal/vaop/ncurrent/abs/nrg3803.html

# Partner Projects #

**Gowinda**: unbiased analysis of gene set enrichement (e.g: Gene Ontology) for Genome Wide Association Studies. Gowinda may thus be used for biological interpretation of the results of PoPoolation and PoPoolation2:  http://code.google.com/p/gowinda/

**PoPoolation**: A pipeline for analyzing pooled next generation sequencing data for single populations. Currently PoPoolation allows to calculate Tajima’s Pi, Watterson’s Theta and Tajima’s D for reference sequences using a sliding window approach http://code.google.com/p/popoolation/

**PoPoolation TE**: A quick and simple pipeline for the analysis of transposable element insertion frequencies in populations from pooled next generation sequencing data. PoPoolation TE identifies TE insertions that are present in the reference genome as well as novel TE insertions and estimates their population frequencies. This also allows for an comparison of TE insertion frequencies between different populations http://code.google.com/p/popoolationte.

**PoPoolation DB**: A user-friendly web-based database for the retrieval of natural variation in Drosophila melanogaster http://www.popoolation.at/pgt/

# Authors #

  * Robert Kofler http://drrobertkofler.wikispaces.com/
  * Christian Schlötterer http://i122server.vu-wien.ac.at/pop/lab_members/christian_schloetterer.html

# microbeGWAS
A fast, lightweight GWAS implementation for bacterial genomes.

This package implements a version of GWAS that is specific to bacterial genomes in order to produce an implementation that is as computationally efficient and streamlined as possible. Partial LD pruning for SNPs in complete LD with each other using the efficient [SNPrune](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0404-z) method described by Calus & Vandenplas is first performed, after which SNPs are fit to a logistic classifier. 

The data used is simulated data from the tool [BacGWASim](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000337). This tool simulates a phylogenetic tree of genomes in order to impose mutation and recombination upon a bacterial genome, resulting in a set of generated genomes and corresponding phenotypes. The pre-generated simulation dataset for [sample size 700](https://figshare.com/articles/bacterial_GWAS_benchmark_simulations_Sample_size_700/9956426), which has exactly 100000 SNPs was used to test the package.

Installation 

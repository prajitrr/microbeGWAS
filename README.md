# microbeGWAS
Prajit Rajkumar, Biology Undergraduate Student

A lightweight GWAS implementation for bacterial genomes.

This package implements a version of GWAS that is specific to bacterial genomes in order to produce an implementation that is as minimal and streamlined as possible. Partial LD pruning for SNPs in complete LD with each other using the efficient [SNPrune](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0404-z) method described by Calus & Vandenplas is first performed, after which SNPs are fit to a logistic classifier. 

The data used is simulated data from the tool [BacGWASim](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000337). This tool simulates a phylogenetic tree of genomes in order to impose mutation and recombination upon a bacterial genome, resulting in a set of generated genomes and corresponding phenotypes. The pre-generated simulation dataset for [sample size 700](https://figshare.com/articles/bacterial_GWAS_benchmark_simulations_Sample_size_700/9956426), which has exactly 100000 SNPs was used to test the package. Ten different phenotypes that are all valid inputs to this tool are also provided in this dataset.

# Installation 
1. Install `bcftools` by following the instructions [here](https://www.htslib.org/download/).
2. Use the following commands to install and setup the package.
   ```
   git clone https://github.com/prajitrr/microbeGWAS
   cd ./microbeGWAS
   pip install requirements.txt
   python setup.py install
   ```
   Additional Considerations: Note that if you do not have root permissions, instead of running `python setup.py install`, you may have to run the commands below.
   ```
   python setup.py install --user
   export PATH="$PATH:$PWD"
   ```
   Also, if you are using a virtual environment such as `anaconda`, the only dependencies other than `bcftools`, `numpy` and `scipy` may come pre-installed.

# Usage
This tool can be used with any bacterial VCF file and a set of phenotypes in the format found [here](https://figshare.com/articles/bacterial_GWAS_benchmark_simulations_Sample_size_700/9956426). Usage is limited to Unix based operating systems.

A sample dataset is provided within the `data` folder of the package. To run the prgoram with the sample dataset, run the following commands, while inside the `microbeGWAS` directory.
```
microbeGWAS ./data/sampleSize_700.vcf.gz ./data/phenRep0.tsv -o output.txt
```
The output produced is a tab-delimited list of all SNPs along with their effect sizes (logistic regression coefficients) and associated p-values. 
If no output file is specified, the program will write to `sys.stdout`.
An optional flag `--prune` or `-p` after the required file inputs can also be added to invoke minimal LD pruning.

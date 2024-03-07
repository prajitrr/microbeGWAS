# microbeGWAS
Prajit Rajkumar, Biology Undergraduate Student

A lightweight GWAS implementation for bacterial genomes.

This package implements a version of GWAS that is specific to bacterial genomes in order to produce an implementation that is as minimal and streamlined as possible. Partial LD pruning for SNPs in complete LD with each other using the efficient [SNPrune](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-018-0404-z) method described by Calus & Vandenplas is first performed, after which SNPs are fit to a logistic classifier. 

The data used is simulated data from the tool [BacGWASim](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000337). This tool simulates a phylogenetic tree of genomes in order to impose mutation and recombination upon a bacterial genome, resulting in a set of generated genomes and corresponding phenotypes. The pre-generated simulation dataset for [sample size 700](https://figshare.com/articles/bacterial_GWAS_benchmark_simulations_Sample_size_700/9956426), which has exactly 100000 SNPs was used to test the package. Ten different phenotypes that are all valid inputs to this tool are also provided in this dataset.

# Installation 
1. Install `bcftools` by following the instructions [here](https://www.htslib.org/download/).
2. Make sure that `numpy`, `pandas`, and `scikit-learn` are installed, using `pip` or another installer.
3. Use the following commands to install and setup the package.
   ```
   git clone https://github.com/prajitrr/microbeGWAS
   cd ./microbeGWAS
   python setup.py install
   ```

# Usage
This tool can be used with any bacterial VCF file and a set of phenotypes in the format found [here](https://figshare.com/articles/bacterial_GWAS_benchmark_simulations_Sample_size_700/9956426). Usage is limited to Unix based operating systems.

Begin by using bcftools to extraction information from the vcf file. Run the following commands.
```
bcftools query -f '[%GT]\n'" ~/PATH_TO_VCF/your_vcf.vcf > snp.txt
```
```
bcftools query -f '%POS\n'" ~/PATH_TO_VCF/your_vcf.vcf > pos.txt
```
The using the following syntax, input the files from the previous step along with the phenotypes file.
```
microbeGWAS ~/PATH_TO_POSITION_FILE/pos.txt ~/PATH_TO_SNP_FILE/snp.txt ~/PATH_TO_PHEN_FILE.tsv -o output.txt
```
The current output produced is a list of all SNPs after minimal LD pruning along with their effect sizes. 
If no output file is specified, the program will write to stdout.
An optional flag `--test` or `t` can also be added to invoke a train-test split of 80:20. A logistic model is first fit on the train data and then tested on the test set, and the F1 score that results is outputted. Otherwise, the entire dataset is used to fit the model.

# TODO
Implement a better way to calculate p-values for the effect sizes.

The current method, which is not implemented here is too RAM intensive due to having to invert a Hessian matrix of a significant size. This will likely involve increasing the severity of LD pruning beyond only complete LD.

Combine the bcftools step into the whole script to make the workflow more streamlined and faster.

Benchmark the tool against plink or another tool.

Port the logistic classifer and pandas file reading over to numpy so that the only dependencies are bcftools and numpy.\\

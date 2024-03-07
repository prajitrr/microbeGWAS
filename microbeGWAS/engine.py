#microbeGWAS/microbeGWAS/engine.py
import numpy as np
import pandas as pd

import argparse
import subprocess

def main():
    #Parse command line arguments
    parser = argparse.ArgumentParser(prog="microbeGWAS", 
                                     description="Command-line tool to perform GWAS on haploid bacterial data."
                                    )
    parser.add_argument("vcf", help="Input VCF file, zipped or unzipped", \
                       type=str)
  
    parser.add_argument("-p", "--phenotypes", help="Specify TSV file of phenotypes for each sample.", \
                       metavar="FILE", type=str, required=True)
    
    parser.add_argument("-o", "--out", help="Write output to specified file path.", \
                       "Default: stdout", metavar="FILE", type=str, required=False)
    
    args = parser.parse_args()

    #Call bcftools to process the VCF file and capture the output

if __name__ == "__main__":
    main()

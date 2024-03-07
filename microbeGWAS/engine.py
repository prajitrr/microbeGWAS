#microbeGWAS/microbeGWAS/engine.py
import numpy as np
import pandas as pd

import argparse
import subprocess
import sys

from sklearn import linear_model
from sklearn.metrics import f1_score

def read_np_matrix(path):
    with open(path, "rb") as f:
        return np.array(
            [np.fromiter((chr(c) for c in l.strip()), dtype=np.uint8) for l in f],
        )

def read_np_array(path):
    with open(path, "rb") as f:
        return np.array(
            [int(l.strip()) for l in f],
        ).astype(np.uint8)

def row_removal(chunk):
    comparison = chunk[:,1:]
    uniques, indices = np.unique(comparison, axis=0, return_index=True)
    uniques_pos = chunk[:, 0][indices]
    return np.concatenate((uniques_pos.reshape(-1,1), uniques), axis=1)

def save_to_file_with_header(file_path, data, header):
    np.savetxt(file_path, data, delimiter='\t', header=header, comments='')

def main():
    #Parse command line arguments
    parser = argparse.ArgumentParser(prog="microbeGWAS", 
                                     description="Command-line tool to perform GWAS on haploid bacterial data."
                                    )
    parser.add_argument("positions", help="Text file of SNP positions.", \
                       type=str)

    parser.add_argument("snps", help="Text file of SNPS.", \
                       type=str)
  
    parser.add_argument("phenotype", help="Specify TSV file of phenotypes for each sample.", \
                        type=str)
    
    parser.add_argument('--test', '-t', action='store_true', help='Set the tool to randomly split input data into train and test data (80/20 split).', required=False)

    parser.add_argument("-o", "--out", help="Write output to specified text file path. Default: stdout", \
                        metavar="FILE", type=str, required=False, default="stdout")
    
    args = parser.parse_args()

    #Call bcftools to process the VCF file and capture the output
    #snp_array = np.genfromtxt(subprocess.run("bcftools query -f '[%GT]\n'" + " " + args.vcf, shell=True, capture_output=True, text=True).stdout, delimiter=1, dtype=int)
    #pos_array = np.genfromtxt(subprocess.run("bcftools query -f '%POS\n'" + " " + args.vcf, shell=True, capture_output=True, text=True).stdout, delimiter=1, dtype=int)

    pos_array = read_np_array(args.positions)
    snp_array = read_np_matrix(args.snps)

    length = len(snp_array[0])
    half_len = np.floor((length)/2)
    
    #Read in phenotypes
    phen_array = (pd.read_csv(args.phenotype, sep='\t', nrows=length)["Phenotype"])
    phen_array[np.isnan(phen_array)] = 255
    phen_array = phen_array.astype(np.uint8)
    

    #Perform LD pruning
    sums = np.sum(snp_array, axis=1).reshape(-1,1)
    min_allele = (sums <= half_len).astype(int)
    sums = sums * np.power(-1, 1 - min_allele) + length * (1 - min_allele)
    
    ld_prune_arr = np.concatenate((pos_array.reshape(-1,1),snp_array, sums, min_allele), axis=1)
    sort_indices = np.argsort(ld_prune_arr[:, length + 1])
    ld_prune_arr = ld_prune_arr[sort_indices]
    ld_chunks = np.split(ld_prune_arr, np.where(np.diff(ld_prune_arr[:,length + 1]))[0]+1)
    
    pruned_chunks = [row_removal(chunk) for chunk in ld_chunks]
    pruned_array = np.vstack(pruned_chunks)
    
    #Run GWAS
    pruned_snp_array = pruned_array[:,1:length].transpose()
    pruned_snp_array =  np.concatenate((pruned_snp_array, phen_array.reshape(-1, 1)), axis=1)
    drops = pruned_snp_array[:, -1] == 255
    pruned_snp_array = pruned_snp_array[~drops]

    model = linear_model.LogisticRegression(C=1e40, solver='newton-cg')
    
    if not args.test:
        X = pruned_snp_array[:,:-1]
        y = pruned_snp_array[:,-1]
        fitted_model = model.fit(X, y)
        data = np.column_stack((pruned_array[:,0], fitted_model.coefs_[0]))
        header = "POS\tBETA"
    else:
        np.random.shuffle(pruned_snp_array)
        train = pruned_snp_array[:round(0.8*len(pruned_snp_array))]
        test = pruned_snp_array[:len(pruned_snp_array) - round(0.8*len(pruned_snp_array))]
        
        X=train[:,:-1]
        y=train[:,-1]
        fitted_model = model.fit(X, y)
        
        f1 = f1_score(test[:,-1],fitted_model.predict(test[:,:-1]))
        f1_col = np.full((len(pruned_array[:,0]),), f1)
        data = np.column_stack((pruned_array[:,0], fitted_model.coefs_[0], f1_col))
        header = "POS\tBETA\tF1"

    if args.out:
        save_to_file_with_header(args.out, data, header)
    else:
        save_to_file_with_header(sys.stdout, data, header)
        

if __name__ == "__main__":
    main()

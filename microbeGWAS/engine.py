#microbeGWAS/microbeGWAS/engine.py
import numpy as np
from scipy.stats import t

import csv
import argparse
import subprocess
import sys

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
    
    parser.add_argument("vcf", help="VCF file of SNPs.", \
                       type=str)
  
    parser.add_argument("phenotype", help="Specify TSV file of phenotypes for each sample.", \
                        type=str)
    
    parser.add_argument('--prune', '-p', action='store_true', help='Enable minimal LD-pruning.', required=False)

    parser.add_argument("-o", "--out", help="Write output to specified text file path. Default: stdout", \
                        metavar="FILE", type=str, required=False, default="stdout")
    
    args = parser.parse_args()

    #Call bcftools to process the VCF file and capture the output
    #snp_array = np.genfromtxt(subprocess.run("bcftools query -f '[%GT]\n'" + " " + args.vcf, shell=True, capture_output=True, text=True).stdout, delimiter=1, dtype=int)
    #pos_array = np.genfromtxt(subprocess.run("bcftools query -f '%POS\n'" + " " + args.vcf, shell=True, capture_output=True, text=True).stdout, delimiter=1, dtype=int)

    #Define bcftools command
    bash_command = "bcftools query -f '[%GT]\n' " + args.vcf
    bash_command_2 = "bcftools query -f '%POS\n' " + args.vcf

    #Call bcftools to extract data from VCF file and 
    process = subprocess.Popen(['bash', '-c', bash_command], stdout=subprocess.PIPE)
    output, _ = process.communicate()
    snps_string = output.decode()
    snps = np.array([np.fromiter(((c) for c in l.strip()), dtype=np.uint8) for l in snps_string.split()],)
    snps = snps.astype(int)
    x = snps

    process = subprocess.Popen(['bash', '-c', bash_command_2], stdout=subprocess.PIPE)
    output, _ = process.communicate()
    positions = output.decode()
    pos_array = np.fromstring(positions, dtype=int, sep='\n')
    
    data = np.genfromtxt(args.phenotype, delimiter='\t', skip_header=1, usecols=1, dtype=np.float16)
    data = data[:len(x.T)]
    nan_indices = np.where(np.isnan(data))[0]
    mask = np.ones(data.shape[0], dtype=np.bool_)
    mask[nan_indices] = False
    data = data[mask].astype(np.bool_)
    y = data
    
    mask = mask[:len(x.T)]
    x = x.T[mask].T

    length = len(x[0])
    half_len = np.floor((length)/2)
    
    #Perform LD pruning
    if args.prune:
        sums = np.sum(x, axis=1).reshape(-1,1)
        min_allele = (sums <= half_len).astype(int)
        sums = sums * np.power(-1, 1 - min_allele) + length * (1 - min_allele)
        
        ld_prune_arr = np.concatenate((pos_array.reshape(-1,1),x, sums, min_allele), axis=1)
        sort_indices = np.argsort(ld_prune_arr[:, length + 1])
        ld_prune_arr = ld_prune_arr[sort_indices]
        ld_chunks = np.split(ld_prune_arr, np.where(np.diff(ld_prune_arr[:,length + 1]))[0]+1)
        
        pruned_chunks = [row_removal(chunk) for chunk in ld_chunks]
        pruned_array = np.vstack(pruned_chunks)
    else:
        pruned_array = np.concatenate((pos_array.reshape(-1,1),x), axis=1)
    
    #Run GWAS
    snp_array = (pruned_array[:,1:]).astype(np.bool_)
    dof = len(y) - 2
    
    n00 = np.sum(~snp_array & ~y.T, axis=1)
    n01 = np.sum(snp_array & ~y.T, axis=1)
    n10 = np.sum(~snp_array & y.T, axis=1)
    n11 = np.sum(snp_array & y.T, axis=1)
    
    p0 = n01/(n00 + n01)
    p1 = n11/(n10 + n11)
    n0 = n00 + n01
    n1 = n10 + n11

    beta = np.log((n11*n00)/(n10*n01 + 1e-100)+1e-100)
    #beta = np.log(n11+1e-100) + np.log(n00+1e-100) - np.log(n10+1e-100) - np.log(n01+1e-100)
    #alpha = np.log(n01/n00)

    p_vals = 2*(1 - t.cdf(np.abs(np.sqrt((p0*(1-p0)*n0 * p1*(1-p1)*n1)/(p0*(1-p0)*n0 + p1*(1-p1)*n1 + 1e-100))*np.log((p1*(1-p0))/(p0*(1-p1)+1e-100)+1e-100)), dof))

    header = "POS\tBETA\tP_VAL"
    
    data = np.column_stack((pruned_array[:,0], beta, p_vals))

    if args.out:
        save_to_file_with_header(args.out, data, header)
    else:
        save_to_file_with_header(sys.stdout, data, header)
        

if __name__ == "__main__":
    main()

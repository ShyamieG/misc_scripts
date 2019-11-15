## SCRIPT: variant_lookup.py
## AUTHOR: SSG
## OVERVIEW: Returns a list of variant IDs (given genomic locations) or genomic positions (given variant IDs) based on a user-provided variant database in VCF format (eg. dbSNP)
## REQUIRED ARGUMENTS:
##         -l, --var_list (list of SNPs to look up, no header, whitespace delimited, 1 column for variant IDs, 2 columns for chromosome # and base pair position)
##         -d, --var_database (a database, like dnSNP, against which to query the SNPs
## OPTIONAL ARGUMENTS:
##         -r, --range (chromosome/base pair range to restrict database query)
##         -p, --print_alleles (whether or not to print the reference and alternate alleles, default 'False')
##         -o, --Output (desired output file name, default "var_lookup_output")
##
## USAGE EXAMPLES:
## $ python variant_lookup.py --var_list my_vars --var_database dbSNP.vcf --range "1:10000-20000"
## $ python variant_lookup.py --var_list my_vars --var_database dbSNP.vcf --range "22"

#! /usr/bin/env python
import argparse as ap
import time
import os
import sys
import logging as log

import numpy as np
import pandas as pd
from cyvcf2 import VCF

### Define main function
def main(args):
    ## Parse arguments
    argp = ap.ArgumentParser(description="")
    argp.add_argument("-l", "--var_list", required=True)
    argp.add_argument("-d", "--var_database", required=True)
    argp.add_argument("-r", "--range")
    argp.add_argument("-p", "--print_alleles", type=bool, default=False)
    argp.add_argument("-o", "--output", default="var_lookup_output")
    args = argp.parse_args(args)
    print(args)

    var_list = pd.read_table(args.var_list, header=None, delim_whitespace=True)
    var_db = VCF(args.var_database)
    outfile=args.output

    ## Check what type of variant list was provided
    if len(var_list.columns) == 1:
        list_type="varID"

    if len(var_list.columns) == 2:
        list_type="var_pos"

    ## Print list of variant positions
    with open(outfile, "a+") as f:
        if list_type == "varID":
            var_list = var_list.values
            for var in var_db(args.range):
                if var.ID in var_list:
                    if args.print_alleles == True:
                        f.write(var.ID + " " + str(var.CHROM) + " " + str(var.POS) + " " + str(var.REF) + " " + str(var.ALT) + "\n")
                    else:
                        f.write(var.ID + " " + str(var.CHROM) + " " + str(var.POS) + "\n")

        elif list_type == "var_pos":
            var_list.columns = ["chr", "pos"]
            var_list['chr'] = var_list.chr.astype('str')
            for var in var_db(args.range):
                if len(var_list[(var_list["pos"]==var.POS) & (var_list["chr"]==var.CHROM)])!=0:
                    if args.print_alleles == True:
                        f.write(var.ID + " " + str(var.CHROM) + " " + str(var.POS) + " " + str(var.REF) + " " + str(var.ALT) + "\n")
                    else:
                        f.write(var.ID + " " + str(var.CHROM) + " " + str(var.POS) + "\n")

    return 0

### Run script
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

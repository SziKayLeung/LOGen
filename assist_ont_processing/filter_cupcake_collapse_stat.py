#!/usr/bin/env python3

"""
Author: Szi Kay Leung (S.K.Leung@exeter.ac.uk)
Date: Oct 2024
Aim: Uses read_stat.txt from IsoSeq3-collapse and abundance to filter and regenerate a new read_stat.txt
"""

import os
import argparse
import pandas as pd

def filter_and_merge(args):
    print("Filtering read count", args.fl_read)
    filtered_abundance = args.abundance.loc[args.abundance["count_fl"] >= args.fl_read,]
    
    filtered_read_stat = args.read_stat.loc[args.read_stat["pbid"].isin(filtered_abundance["pbid"].values),]

    # assert errors
    assert set(filtered_read_stat["pbid"].values).issubset(filtered_abundance["pbid"].values)
    assert set(filtered_abundance["pbid"].values).issubset(filtered_read_stat["pbid"].values)

    print("Writing output:", args.output)
    filtered_read_stat.to_csv(args.output, index=False, header=True, sep="\t", encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description="ONT demultiplex and retrieve read counts from IsoSeq3 collapse output")
    parser.add_argument("read_stat", help='\t\tread_stat.txt generated from Cupcake/Isoseq3 collapse')
    parser.add_argument('abundance', help='\t\tabundance.txt generated from Cupcake/Isoseq3 collapse')
    parser.add_argument('-c','--fl_read', required=True, type=int, help='\t\tMinimum number of FL reads to filter')
  
    args = parser.parse_args()
    args.output = args.read_stat.replace('.txt', '') + ".min" + str(args.fl_read) + "FL.txt" 
    
        
    # read in files
    print("Reading in files...")
    args.read_stat = pd.read_csv(args.read_stat, sep = "\t")
    args.abundance = pd.read_csv(args.abundance, comment='#', sep = "\t")

    filter_and_merge(args)
    print("All Done")
    

if __name__ == "__main__":
    main()
    
